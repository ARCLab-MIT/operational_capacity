import math
import numpy as np

# Constants
earth_radius = 6378.137  # km
mu = 398600.435507  # km^3/s^2

# Example function to compute mean motion (you can replace this with your calculation)
def calculate_mean_motion(semi_major_axis):
    # Calculate the orbital period (T) in seconds
    period = 2 * math.pi * math.sqrt(semi_major_axis**3 / mu)

    # Mean motion (n) in radians per second
    mean_motion = 2 * math.pi / period  # radians per second

    # Convert mean motion to revolutions per day
    mean_motion_per_day = mean_motion * 86400 / (2 * math.pi)  # revolutions per day
    
    return mean_motion_per_day

# TLE Generator
def generate_tle(name, inclination, semi_major_axis, norad_id):
    # Format inclination (the second line has inclination in degrees)
    inclination_str = f"{inclination:8.4f}"  # Format inclination to 4 decimal places
    
    # Format mean motion (the second line has mean motion)
    mean_motion = calculate_mean_motion(semi_major_axis)
    mean_motion_str = f"{mean_motion:8.4f}"  # Format mean motion to 4 decimal places
    
    # Create TLE lines, adjusting for the new norad_id
    tle_line_1 = f"1 {norad_id:05d}U 59001A   24268.51499735  .00001072  00000-0  57112-3 0  9991"
    tle_line_2 = f"2 {norad_id:05d} {inclination_str} 309.4174 0  71.7704 303.5363 {mean_motion_str}"
    
    return name, tle_line_1, tle_line_2

# Generate TLEs with varying inclination and mean motion
def generate_and_write_tles(file_name, inclination_range, semi_major_axis_range):
    norad_id = 90000  # Start with norad_id = 90000
    with open(file_name, 'w') as f:
        for inclination in inclination_range:
            for semi_major_axis in semi_major_axis_range:
                # Generate TLE for each combination of inclination and mean motion
                name = f"Alt_{semi_major_axis}_Inc_{inclination}"
                _, tle_line_1, tle_line_2 = generate_tle(name, inclination, semi_major_axis, norad_id)
                
                # Write the TLE to the file
                f.write(f"{name}\n")
                f.write(f"{tle_line_1}\n")
                f.write(f"{tle_line_2}\n")
                
                # Increment norad_id for the next TLE
                norad_id += 1

# Define ranges for inclination and semi-major axis
inclination_range = np.arange(0, 111, 5)  # Whole number degrees from 0 to 110 degrees
semi_major_axis_range = np.arange(300, 2001, 10) + earth_radius  # Semi-major axis from 300 km to 2000 km at 50 km intervals
print(semi_major_axis_range)

# Generate TLEs and write to a file
generate_and_write_tles("generated_tles.txt", inclination_range, semi_major_axis_range)
