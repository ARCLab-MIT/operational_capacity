'''
import spacetrack.operators as op
from spacetrack import SpaceTrackClient
import datetime as dt
import numpy as np
from time import sleep
import csv
from pathlib import Path
from random import sample

def main(): 
    # The goal is to pull 1 week of tles every year from 1970 to modern day. Use jan 1 as the initiate date of each year
    years = np.arange(1990,1992,1)
    months = [1]#np.arange(1,12,1)
    day = 1
    # years = np.array([2022])
    for year in years: 
        for month in months:
            pull_1wk_tles(year, month, day)
            print(str(month)+'/'+str(year))


def pull_1wk_tles(year,month,day): 
    # pulls all tles 1 week forward from the date of interest 
    # Check whether ./Data/TLEs/ exists. If not, make it.
    path = Path(__file__).resolve().parent / "Data/Historical_TLE_Data"
    Path(path).mkdir(parents=True, exist_ok=True)

    # Prompt user for Space-Track log-in credentials
    print("Log in to using your Space-Track.org account credentials.\n")
    st_email = input("Space-Track Email: ")
    st_pass =  input("Space-Track Password: ")

    # Log in to Space-Track using your email and password
    st = SpaceTrackClient(identity=st_email, password=st_pass)

    print("Define a study period (using yyyy-mm-dd formats).\n")
    # studyperiod_start = start_date #input("Start date: ")
    # studyperiod_end = '1999-03-30'#input("End date: ")

    # Only pull TLEs from between the start date and the end date
    # start_datetime = dt.datetime(int(studyperiod_start[0:4]), int(studyperiod_start[5:7]), int(studyperiod_start[8:10]))
    start_datetime = dt.datetime(year, month, day)
    one_week = dt.timedelta(weeks = 1) 
    end_datetime = start_datetime + one_week
    drange = op.inclusive_range(start_datetime, end_datetime)

    try: 
        #file_name = f"{year}{month:02d}{day:02d}_tle.txt"
        #file_path = base_path / file_name
        with open('../Data/Historical_TLE_Data/'+str(year)+str(month)+str(day) +'_oneweek.txt','w') as f: 
            response = st.gp(epoch=drange, orderby='epoch desc', format='tle')
            # print(response)
            f.write(response)
            f.close()
            print("Printed file")
            for i in range(26):
                if i < 13:
                    print("*"*(i+1))
                    sleep(0.5)
                else:
                    print("*"*(26-i))
                    sleep(0.5)
    except:
        print('Couldnt pull data from ' + str(year) + str(month) + str(day))


if __name__ == "__main__":
    main()
'''

import spacetrack.operators as op
from spacetrack import SpaceTrackClient
import datetime as dt
import numpy as np
from time import sleep
import csv
from pathlib import Path
import os


def main():
    # The goal is to pull TLEs on the first day of every month from 1990 to modern day.
    years = np.arange(2023, 2024, 1)
    month = 1  # All 12 months
    day = 1

    # Prompt user for Space-Track log-in credentials once
    print("Log in using your Space-Track.org account credentials.\n")
    st_email = input("Space-Track Email: ")
    st_pass = input("Space-Track Password: ")
    
    # Log in to Space-Track
    st = SpaceTrackClient(identity=st_email, password=st_pass)
    
    for year in years:
        pull_tles(st, year, month, day)
        print(f"{month}/{year}")


def pull_tles(st, year, month, day):
    # Ensure files are saved inside the project folder
    base_path = Path(__file__).resolve().parent / "Data/Historical_TLE_Data"
    base_path.mkdir(parents=True, exist_ok=True)

    #tle_date = dt.datetime(year, month, day).strftime("%Y-%m-%d")
    start_datetime = dt.datetime(year, month, day)
    one_week = dt.timedelta(weeks = 1) 
    end_datetime = start_datetime + one_week
    drange = op.inclusive_range(start_datetime, end_datetime)

    try:
        file_name = f"{year}{month:02d}{day:02d}_tle.csv"
        file_path = base_path / file_name
        response = st.gp(epoch=drange, orderby='epoch desc', format='csv')
        #response = st.gp(epoch=tle_date, orderby='epoch desc', format='tle')
        
        if response.strip():  # Ensure response is not empty
            with open(file_path, 'w') as f:
                f.write(response)
                print(f"Printed file: {file_path}")
                for i in range(26):
                    if i < 13:
                        print("*"*(i+1))
                        sleep(0.5)
                    else:
                        print("*"*(26-i))
                        sleep(0.5)
        else:
            print(f"No data received for {year}-{month:02d}-{day:02d}, skipping file.")
    
    except Exception as e:
        print(f"Couldn't pull data from {year}-{month:02d}-{day:02d}: {e}")



if __name__ == "__main__":
    main()
