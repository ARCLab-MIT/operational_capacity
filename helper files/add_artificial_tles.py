import os

input_dir = "data/batch_tles_monthly"
output_dir = "data/batch_tles_monthly_random"
starlink_file = "generated_tles_random.txt"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load TLEs from generated_starlink.txt
with open(starlink_file, 'r') as f:
    starlink_tles = f.read()

# Loop from Jan 2019 to Dec 2024
for year in range(2019, 2025):
    for month in range(1, 13):
        filename = f"{month:02d}{year}.txt"
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)

        if not os.path.exists(input_path):
            print(f"File not found: {input_path}")
            continue

        with open(input_path, 'r') as infile:
            lines = infile.readlines()

        filtered_tles = []
        i = 0
        while i < len(lines):
            if lines[i].startswith("Alt_"):
                i += 3  # Skip this TLE block
            else:
                filtered_tles.extend(lines[i:i+3])
                i += 3

        # Combine and write to output
        with open(output_path, 'w') as outfile:
            outfile.writelines(filtered_tles)
            outfile.write(starlink_tles)

        print(f"Saved: {output_path}")

# def filter_starlink_tles(input_file, output_file):
#     with open(input_file, 'r') as infile:
#         lines = infile.readlines()

#     filtered_tles = []

#     # Iterate in chunks of 3 lines (name, line1, line2)
#     for i in range(0, len(lines), 3):
#         if i + 2 >= len(lines):
#             break  # Skip incomplete TLE sets

#         name_line = lines[i].strip()
#         line1 = lines[i + 1]
#         line2 = lines[i + 2]

#         if "STARLINK" in name_line:
#             filtered_tles.extend([name_line + "\n", line1, line2])

#     with open(output_file, 'w') as outfile:
#         outfile.writelines(filtered_tles)

# # Example usage
# filter_starlink_tles('data/batch_tles_monthly/122024.txt', 'data/batch_tles_monthly/filtered_starlink_tles.txt')

