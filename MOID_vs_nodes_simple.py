import numpy as np
from skyfield.api import load, EarthSatellite
from scipy.spatial.distance import cdist
from datetime import timedelta
from tle_conjunctions_wp_v2 import find_min_conj_dist

def approximate_moid_sgp4(sat1, sat2, t0, duration_seconds=2*60*60, step_seconds=10):
    ts = load.timescale()

    # Create list of Skyfield Time objects spaced by step_seconds
    dt_list = [t0.utc_datetime() + timedelta(seconds=i) for i in range(0, duration_seconds, step_seconds)]
    times = ts.utc([dt.year for dt in dt_list],
                   [dt.month for dt in dt_list],
                   [dt.day for dt in dt_list],
                   [dt.hour for dt in dt_list],
                   [dt.minute for dt in dt_list],
                   [dt.second for dt in dt_list])
    
    pos1 = np.array([sat1.at(t).position.km for t in times])
    pos2 = np.array([sat2.at(t).position.km for t in times])

    dist_matrix = cdist(pos1, pos2)  # shape (N, N)
    min_dist = np.min(dist_matrix)
    return min_dist

import random
from datetime import datetime
from skyfield.api import EarthSatellite, load

def generate_random_tles(epoch=None):
    """
    Generate two synthetic but valid TLEs for use with Skyfield or SGP4.
    The TLEs will correspond to satellites in nearby orbits with randomized parameters.

    Returns:
        (tle1_lines, tle2_lines): Tuple of two TLE 2-line string pairs
    """
    if epoch is None:
        now = datetime.utcnow()
        epoch_year = now.year % 100
        epoch_day = now.timetuple().tm_yday + now.hour / 24.0
    else:
        epoch_year = epoch.year % 100
        epoch_day = epoch.timetuple().tm_yday + epoch.hour / 24.0

    def make_tle(satnum):
        inclination = round(random.uniform(45, 100), 4)
        raan = round(random.uniform(0, 360), 4)
        ecc = round(random.uniform(0.01, 0.05), 7)  # 7-digit eccentricity
        # ecc = round(random.uniform(0.0001, 0.01), 7)  # 7-digit eccentricity
        argp = round(random.uniform(0, 360), 4)
        mean_anomaly = round(random.uniform(0, 360), 4)
        mean_motion = round(random.uniform(14.0, 15.8), 8)  # ~500–900 km orbit
        bstar = f"{random.randint(10000, 99999):05d}-3"

        line1 = f"1 {satnum:05d}U 24001A   {epoch_year:02d}{epoch_day:012.8f}  .00000000  00000-0  {bstar} 0  9990"
        line2 = f"2 {satnum:05d} {inclination:8.4f} {raan:8.4f} {int(ecc*1e7):07d} {argp:8.4f} {mean_anomaly:8.4f} {mean_motion:11.8f}"
        return line1, line2

    tle1 = make_tle(random.randint(10000, 99999))
    tle2 = make_tle(random.randint(10000, 99999))

    tle1_line1, tle1_line2 = tle1
    tle2_line1, tle2_line2 = tle2

    return tle1, tle2, tle1_line1, tle1_line2, tle2_line1, tle2_line2




n = 500
moid_km = np.zeros(n)
min_nodes = np.zeros(n)
for i in range(n):
    tle1, tle2, tle1_line1, tle1_line2, tle2_line1, tle2_line2 = generate_random_tles()

    ts = load.timescale()
    t0 = ts.utc(2025, 6, 18)

    sat1 = EarthSatellite(*tle1)
    sat2 = EarthSatellite(*tle2)

    moid_km[i] = approximate_moid_sgp4(sat1, sat2, t0)
    print(f"Approximate MOID ≈ {moid_km[i]:.2f} km")

    # find analytical nodal conjunction distance
    d_a, d_b, c,d,e = find_min_conj_dist(tle1_line1, tle1_line2, tle2_line1, tle2_line2)
    min_nodes[i] = min(d_a, d_b)
    print(f"Analytical Nodal Conjunction Distance ≈ {min_nodes[i]:.2f} km")

#plot
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 5))
plt.scatter(moid_km, moid_km-min_nodes, marker='o', color='b')
plt.title('Approximate MOID between Randomly Generated Satellites')
plt.xlabel('MOID (km)')
plt.ylabel('diff MOID vs Nodal Conjunction Distance (km)')
plt.tight_layout()
plt.show()

