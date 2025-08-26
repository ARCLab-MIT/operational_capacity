import numpy as np

# Replace this with your actual file path
tle_file_path = 'data/batch_tles_monthly_random/122024.txt'

mean_anomalies = []
raans = []
args_of_perigee = []

with open(tle_file_path, 'r') as file:
    lines = file.readlines()

for i in range(len(lines)):
    if lines[i].startswith('Alt_'):
        continue  # Skip metadata lines

    if not lines[i].startswith('Alt_') and lines[i].startswith('1 ') and (i + 1) < len(lines):
        line2 = lines[i + 1].strip()
        if line2.startswith('2 '):
            try:
                # Extract fixed-width substrings (1-based positions in TLE format)
                raan = float(line2[17:25].strip())
                argp = float(line2[34:42].strip())
                mean_anom = float(line2[43:51].strip())
                mean_motion = float(line2[52:63].strip())
                
                if mean_motion < 15.1203 and mean_motion > 15.0549:
                    raans.append(raan)
                    args_of_perigee.append(argp)
                    mean_anomalies.append(mean_anom)

            except ValueError:
                    continue  # Skip bad lines

# Convert to NumPy arrays
raans = np.array(raans)
args_of_perigee = np.array(args_of_perigee)
mean_anomalies = np.array(mean_anomalies)

# Helper to format stats
def print_stats(name, values):
    print(f"\n{name}:")
    print(f"  Mean: {np.mean(values):.3f}")
    print(f"  Min : {np.min(values):.3f}")
    print(f"  Max : {np.max(values):.3f}")
    print(f"  Std : {np.std(values):.3f}")

# Print results
print_stats("RAAN", raans)
print_stats("Argument of Perigee", args_of_perigee)
print_stats("Mean Anomaly", mean_anomalies)
