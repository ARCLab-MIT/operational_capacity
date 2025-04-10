import spacetrack.operators as op
from spacetrack import SpaceTrackClient
import datetime as dt
import numpy as np
from time import sleep
from pathlib import Path

# Delay between API calls to stay under 30 queries per minute
QUERY_DELAY = 2  # seconds


def main(): 
    # Prompt user for Space-Track log-in credentials
    print("Log in to Space-Track.org with your account credentials.\n")
    st_email = input("Space-Track Email: ")
    st_pass = input("Space-Track Password: ")

    # The goal is to pull 1 week of TLEs every year from 1970 to modern day. Use Jan 1 as the starting date of each year
    years = np.arange(2024, 2025, 1)
    months = np.arange(1, 2, 1)  # Use specific months or range: np.arange(1, 13, 1)
    day = 1
    for year in years: 
        for month in months:
            if year == 2025 and month > 2:
                break
            pull_1wk_tles(year, month, day, st_email, st_pass)
            print(str(month) + '/' + str(year))


def pull_1wk_tles(year, month, day, st_email, st_pass): 
    # Pulls all TLEs 1 week forward from the date of interest
    path = 'data/batch_tles_annual/'
    Path(path).mkdir(parents=True, exist_ok=True)

    # Log in to Space-Track using your email and password
    st = SpaceTrackClient(identity=st_email, password=st_pass)

    start_datetime = dt.datetime(year, month, day)
    one_week = dt.timedelta(weeks=1)
    end_datetime = start_datetime + one_week

    drange = op.inclusive_range(start_datetime, end_datetime)
    #raw_output_file = f'{path}{month:02d}{year}_raw.txt'  # Raw data file
    raw_output_file = f'{path}{year}_raw.txt'  # Raw data file
    #clean_output_file = f'{path}{month:02d}{year}.txt'  # Cleaned file
    clean_output_file = f'{path}{year}.txt'  # Cleaned file
    generated_tles_file = Path(__file__).parent / "generated_tles.txt"  # Path to generated TLEs


    try:
        # Attempt to pull one week of data
        with open(raw_output_file, 'w') as f: 
            response = st.gp_history(epoch=drange, mean_motion=op.greater_than(8), orderby="NORAD_CAT_ID desc", format="3le")
            f.write(response)
            print("Printed file successfully!")
            progress_bar()
            #sleep(QUERY_DELAY)  # Add delay after the query
            remove_duplicate_tles(raw_output_file, clean_output_file)
            print(f"Duplicate TLEs removed. Clean file saved as {clean_output_file}")

            # Append cleaned TLEs to generated_tles.txt
            with open(generated_tles_file, "r") as gen_f, open(clean_output_file, "a") as clean_f:
                clean_f.write(gen_f.read())  # Append generated TLEs to the cleaned TLEs
            print(f"Appended {generated_tles_file} onto {clean_output_file}")

    except:
        print(f"Failed to pull 1 week of data for {year}-{month:02d}-{day:02d}. Trying smaller chunks.")
        pull_smaller_chunks(st, start_datetime, end_datetime, raw_output_file)
        remove_duplicate_tles(raw_output_file, clean_output_file)

        # Append cleaned data to generated_tles.txt
        with open(generated_tles_file, "r") as gen_f, open(clean_output_file, "a") as clean_f:
            clean_f.write(gen_f.read())  # Append generated TLEs to the cleaned TLEs
        print(f"Appended {generated_tles_file} onto {clean_output_file}")


def pull_smaller_chunks(st, start_datetime, end_datetime, output_file):
    # Pulls data in smaller chunks (2 days at a time) and combines results
    temp_files = []
    chunk_size = dt.timedelta(days=1)
    current_start = start_datetime

    while current_start < end_datetime:
        current_end = min(current_start + chunk_size, end_datetime)
        drange = op.inclusive_range(current_start, current_end)
        temp_file = f'temp_{current_start.strftime("%Y%m%d")}_{current_end.strftime("%Y%m%d")}.txt'
        try:
            with open(temp_file, 'w') as f:
                response = st.gp_history(epoch=drange, mean_motion=op.greater_than(8), orderby="NORAD_CAT_ID desc", format="3le")
                f.write(response)
                temp_files.append(temp_file)
                print(f"Pulled data for {current_start} to {current_end}")
                progress_bar()
                #sleep(QUERY_DELAY)  # Add delay after each query
        except:
            print(f"Failed to pull data for {current_start} to {current_end}")
        current_start += chunk_size

    # Combine all temporary files into the final output
    with open(output_file, 'w') as outfile:
        for temp_file in temp_files:
            try:
                with open(temp_file, 'r') as infile:
                    outfile.write(infile.read())
                Path(temp_file).unlink()  # Delete temporary file after merging
            except:
                print(f"Could not read or delete temporary file: {temp_file}")
    print(f"Combined data into {output_file}")


def progress_bar():
    for i in range(26):
        if i < 13:
            print("*" * (i + 1))
        else:
            print("*" * (26 - i))
        sleep(0.1)


def remove_duplicate_tles(input_file, output_file):
    # Reads a file of TLEs and removes duplicate entries
    seen_satellites = set()
    unique_tles = []
    
    with open(input_file, "r") as file:
        lines = file.readlines()
    
    i = 0
    while i < len(lines):
        if lines[i].startswith("0 "):  # Satellite name line
            sat_name = lines[i].strip()

            #removes satellites with TBA - TO BE ASSIGNED
            if sat_name != "0 TBA - TO BE ASSIGNED" and sat_name not in seen_satellites:
                seen_satellites.add(sat_name)
                unique_tles.extend(lines[i:i+3])  # Keep this TLE set
            # Skip adding TLEs for "TBA - TO BE ASSIGNED"
            i += 3  # Move to the next TLE set
        else:
            i += 1  # Skip any unexpected lines
    
    with open(output_file, "w") as file:
        file.writelines(unique_tles)



if __name__ == "__main__":
    main()
