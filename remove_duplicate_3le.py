from pathlib import Path

def remove_duplicate_tles(input_file, output_file):
    seen_satellites = set()
    unique_tles = []
    
    with open(input_file, "r") as file:
        lines = file.readlines()
    
    i = 0
    while i < len(lines):
        if lines[i].startswith("0 "):  # Satellite name line
            sat_name = lines[i].strip()
            if sat_name not in seen_satellites:
                seen_satellites.add(sat_name)
                unique_tles.extend(lines[i:i+3])  # Keep this TLE set
            i += 3  # Move to the next TLE set
        else:
            i += 1  # Skip any unexpected lines
    
    with open(output_file, "w") as file:
        file.writelines(unique_tles)

# Define the input and output file paths
year = 2025
input_file = Path(__file__).resolve().parent / ("data/Historical_TLE_Data/" + str(year) + ".txt")
output_file = Path(__file__).resolve().parent / ("data/Historical_TLE_Data/" + str(year) + "_no_repeats.txt")

# Example usage
remove_duplicate_tles(input_file, output_file)

print(f"Processed file saved to {output_file}")
