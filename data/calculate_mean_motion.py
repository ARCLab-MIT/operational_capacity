import math
earth_radius = 6378.137  # km
mu = 398600.435507  # km^3/s^2
import random

def calculate_mean_motion(semi_major_axis):
    # Calculate the orbital period (T) in seconds
    period = 2 * math.pi * math.sqrt(semi_major_axis**3 / mu)

    # Mean motion (n) in radians per second
    mean_motion = 2 * math.pi / period  # radians per second

    # Convert mean motion to revolutions per day
    mean_motion_per_day = mean_motion * 86400 / (2 * math.pi)  # revolutions per day
    
    return mean_motion_per_day

print(calculate_mean_motion(1275 + earth_radius))

def count_tles_in_mean_motion_range(tle_file_path, lower_bound=12.967, upper_bound=13.096):
    count = 0

    with open(tle_file_path, 'r') as file:
        lines = file.readlines()

    # Process TLEs in 3-line sets (name + line 1 + line 2)
    for i in range(0, len(lines) - 2, 3):
        line2 = lines[i + 2].strip()

        if line2.startswith('2 '):
            try:
                mean_motion = float(line2[52:63])
                if lower_bound <= mean_motion <= upper_bound:
                    count += 1
            except ValueError:
                continue  # Skip lines with malformed float values

    return count

# Example usage
tle_file = 'data/batch_tles_annual/2023.txt'  # Replace with the actual path to your file
matching_tle_count = count_tles_in_mean_motion_range(tle_file)
print(f"Number of TLEs with mean motion between 15.136 and 15.303: {matching_tle_count}")


# def count_matching_tles(file_path):
#     count = 0

#     with open(file_path, 'r') as file:
#         lines = file.readlines()

#     for i in range(0, len(lines), 3):
#         if i + 2 >= len(lines):
#             break  # Avoid index error if file isn't perfectly divisible by 3

#         name_line = lines[i].strip()
#         line1 = lines[i+1].strip()
#         line2 = lines[i+2].strip()

#         try:

#             norad_number = int(line1[2:7])

#             inclination = float(line2[8:16])

#             if norad_number >= 90000:
#                 mean_motion = float(line2[46:54])
#             else:
#                 mean_motion = float(line2[52:63])

#             if 49.99 <= inclination <= 55.01 and 15.022 <= mean_motion <= 15.088:
#                 print(f"Match: {name_line}")  # Print satellite name
#                 count += 1

#         except ValueError:
#             continue  # Skip malformed lines

#     return count

# # Example usage:
# file_path = "data/batch_tles_annual/2000.txt"
# result = count_matching_tles(file_path)
# print(f"Number of matching TLEs: {result}")

# monthly_occupation_550 = [0, 0, 0, 0, 0.5097190711919235, 0, 0, 0, 0, 0, 1.8705277024930076, 0, 0, 0, 1.0054648029715363, 0.47683275231726846, 0.5291426166626773, 0.4532509741286715, 0, 0, 0, 1.5325679571017419, 1.080748823036795, 0, 2.750945597085564, 5.239196835263681, 5.158337414445081, 9.83008112925603, 13.044093884427276, 10.545077005669363, 4.126338335987661, 9.913526514722198, 11.920536909309817, 6.877229847066196, 11.003621875851685, 8.991967494058654, 12.379042150733554, 11.161112824945631, 8.71120420948541, 9.628180251292015, 9.169597144813487, 8.71118484164103, 11.038033278036014, 6.178727905320794, 14.107450620722988, 8.711207386123846, 11.462057946100272, 11.462099214835503, 8.710489076545764, 9.881082662161367, 13.75385777277527, 9.169749761087633, 5.9598034031770135, 14.872875211570035, 10.283181366835558, 12.883207872503666, 8.331493284291703, 7.794263135014882, 14.743050062221176, 9.628271441371016, 15.48662094910458, 12.995244938016649, 9.839848312414038, 7.335751711577553, 11.462182688473696, 9.2658610158474, 11.462077161338055, 9.169708814404412, 12.516482899847203, 10.54527487103691, 5.960323689655398, 7.794251497499646]
# monthly_population_550 = [6, 6, 6, 6, 6, 6, 39, 7, 58, 57, 57, 12, 26, 26, 63, 138, 175, 216, 291, 369, 425, 476, 529, 598, 671, 767, 825, 880, 927, 1043, 1213, 1406, 1456, 1463, 1497, 1502, 1498, 1489, 1486, 1492, 1491, 1480, 1480, 1491, 1484, 1484, 1499, 1474, 1471, 1479, 1593, 1476, 1478, 1455, 1454, 1444, 1524, 1464, 1436, 1433, 1429, 1439, 1435, 1374, 1360, 1355, 1342, 1331, 1308, 1288, 1186, 1021]
# annual_occupation_550 = [0, 0, 0, 0, 0, 0, 0, 0, 0.9288803993673741, 0, 0, 0, 0, 1.1206697298186084, 0, 0, 0, 0, 0.5019032005216694]
# annual_population_550 = [7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]

# import pandas as pd
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates

# # Define your start and end
# start_year = 2000
# end_year = 2025  # because we want to include Dec 2024

# # Create x-axis time vectors
# plot_years = pd.date_range(start=f'{start_year}-01-01', end='2018-01-01', freq='YS')  # Annual data: 2000–2018
# plot_months = pd.date_range(start='2019-01-01', end='2024-12-01', freq='MS')          # Monthly data: Jan 2019–Dec 2024

# # Combine occupation and population data
# full_dates = list(plot_years) + list(plot_months)
# full_occupation = annual_occupation_550 + monthly_occupation_550
# full_population = annual_population_550 + monthly_population_550

# # Set up plot
# fig, ax1 = plt.subplots(figsize=(14, 6))

# # Left y-axis: Conjunction Frequency
# color1 = 'tab:blue'
# ax1.set_xlabel('Date')
# ax1.set_ylabel('Conjunction Frequency at 550 km', color=color1)
# ax1.plot(full_dates, full_occupation, color=color1, label='Conjunction Frequency', linestyle='-')
# ax1.tick_params(axis='y', labelcolor=color1)

# # Right y-axis: Number of Objects
# ax2 = ax1.twinx()
# color2 = 'tab:red'
# ax2.set_ylabel('Number of Objects at 550 km', color=color2)
# ax2.plot(full_dates, full_population, color=color2, label='Number of Objects', linestyle='-')
# ax2.tick_params(axis='y', labelcolor=color2)

# # Format x-axis
# ax1.set_xlim([pd.Timestamp(f'{start_year}-01-01'), pd.Timestamp('2024-12-31')])
# ax1.xaxis.set_major_locator(mdates.YearLocator())
# ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

# # Optional: rotate labels
# for label in ax1.get_xticklabels():
#     label.set_rotation(45)

# # Title, grid, layout
# plt.title('Conjunction Frequency vs Number of Objects at 550 km (2000–2024)')
# ax1.grid(True)
# fig.tight_layout()
# #plt.show()
