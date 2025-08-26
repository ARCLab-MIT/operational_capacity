from pathlib import Path
import os

def remove_duplicate_tles(input_file, output_file):
    seen_ids = set()
    unique_tles = []
    
    with open(input_file, "r") as file:
        lines = file.readlines()
    
    i = 0
    while i + 2 < len(lines):  # Ensure full 3-line TLE
        if lines[i].startswith("0 ") and lines[i+1].startswith("1 ") and lines[i+2].startswith("2 "):
            line1 = lines[i+1]
            norad_id = line1[2:7].strip()
            
            # Skip TLEs where NORAD ID starts with 'T'
            if norad_id.upper().startswith('T'):
                i += 3
                continue
            
            if norad_id not in seen_ids:
                seen_ids.add(norad_id)
                unique_tles.extend(lines[i:i+3])  # Add the full TLE set
        i += 3  # Move to next TLE
    
    with open(output_file, "w") as file:
        file.writelines(unique_tles)

input_dir = "data/batch_tles_annual"
output_dir = "data/batch_tles_annual"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Loop from Jan 2019 to Dec 2024
for year in range(2019, 2026):
    for month in range(1, 2):
        #filename_in = f"{month:02d}{year}_raw.txt"
        #filename_out = f"{month:02d}{year}.txt"
        filename_in = f"{year}_raw.txt"
        filename_out = f"{year}.txt"
        input_path = os.path.join(input_dir, filename_in)
        output_path = os.path.join(output_dir, filename_out)
        remove_duplicate_tles(input_path, output_path)