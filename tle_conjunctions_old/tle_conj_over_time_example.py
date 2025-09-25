import numpy as np
from skyfield.api import load, EarthSatellite
from scipy.spatial.distance import cdist
from datetime import timedelta
from tle_conjunctions_wp_v2 import find_min_conj_dist_v2, time_until_collision_v2, get_satellite
import csv
import datetime as dt
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'

def compute_sat_params(tle_l1, tle_l2):

    # get all the satellite orbit parameters from the TLEs. Store in a dictionary with the satellite norad id as the key
    earth_radius = 6378.15
    mu = 398600
    satellite = get_satellite(tle_l1[i], tle_l2[i])
    a = satellite.model.a * earth_radius  # Semi-major axis in km
    e = satellite.model.ecco  # Eccentricity
    a_range = [a*(1-e), a*(1+e)]
    incl = satellite.model.inclo  # Inclination in radians
    raan = satellite.model.nodeo  # Right ascension of ascending node in radians
    arg_per = satellite.model.argpo  # Argument of perigee in radians
    meanan = satellite.model.mo  # Mean anomaly at epoch in radians
    epoch = satellite.epoch.utc_datetime()  # Epoch time
    # compute orbit period using a and mu
    period = 2 * np.pi * np.sqrt(a ** 3/mu)

    sat_params = {
        'a': a,
        'a_range': a_range,
        'e': e,
        'incl': incl,
        'raan': raan,
        'arg_per': arg_per,
        'meanan': meanan,
        'epoch': epoch,
        'period': period,
    }

    return sat_params

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

# load and read csv file
def read_tle_csv(file_path):
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        lines = list(reader)
        line1 = [l[23] for l in lines]
        line2 = [l[24] for l in lines]
        epoch = [l[7] for l in lines]
        # convert epoch to datetime from format: '2024-01-01 09:02:47'
        if file_path == 'data/COSMOS2221_20250624_1571652595.csv':
            epoch = [dt.datetime.strptime(e, '%Y-%m-%d %H:%M:%S') for e in epoch[1:]]
        elif file_path == 'data/TIMED_20250624_1809607376.csv':
            epoch = [dt.datetime.strptime(e, '%m/%d/%y %H:%M') for e in epoch[1:]]

    return line1[1:], line2[1:], epoch


s1l1, s1l2, s1_epoch = read_tle_csv('data/TIMED_20250624_1809607376.csv')
s2l1, s2l2, s2_epoch = read_tle_csv('data/COSMOS2221_20250624_1571652595.csv')
sat_params_s1 = []
sat_params_s2 = []
print('Computing sat params from TLEs!')
for i in range(len(s1l1)):
    sat_params_s1.append(compute_sat_params(s1l1, s1l2))
    print(i/len(s1l1))
for i in range(len(s2l1)):
    sat_params_s2.append(compute_sat_params(s2l1, s2l2))
    print(i/len(s2l1))


s1_epoch_ts = [(s1_epoch[i]).timestamp() for i in range(len(s1_epoch))]
s2_epoch_ts = [(s2_epoch[i]).timestamp() for i in range(len(s2_epoch))]
s1_epoch_ts = np.array(s1_epoch_ts)
s2_epoch_ts = np.array(s2_epoch_ts)

n = 1000 # 1000 timesteps throughout the year
tvec = np.linspace(0, 365*24*60*60, n) # 1 year in seconds
t = [dt.datetime(2024,1,1) + dt.timedelta(seconds=tvec[i]) for i in range(len(tvec))]
t_ts = [(t[i]).timestamp() for i in range(n)]
t_ts = np.array(t_ts)

# for each timestep, compute the analytical nodal conjunction properties
# find analytical nodal conjunction distance
d_a = np.zeros((len(t),))
d_b = np.zeros((len(t),))
alt_a = np.zeros((len(t),))
alt_b = np.zeros((len(t),))
T_lcm = np.zeros((len(t),))
ttc_a = np.zeros((len(t),))
ttc_b = np.zeros((len(t),))
int_pt_a = np.zeros((len(t), 3))  # 3D point of intersection
int_pt_b = np.zeros((len(t), 3))  # 3D point of intersection
col_a = np.zeros((len(t),), dtype=bool)  # collision flag for sat1
col_b = np.zeros((len(t),), dtype=bool)  # collision flag for sat2
for i in range(len(t)):
    # get index for sat1 and sat2 that is the most recent before time t
    diff1 = t_ts[i] - s1_epoch_ts
    # diff1[diff1 < 0] == 1000000 # fix this later!
    idx1 = np.argmin(np.abs(diff1))

    diff2 = t_ts[i]-s2_epoch_ts
    # diff2[diff2 < 0] == 1000000 # fix this later!
    idx2 = np.argmin(np.abs(diff2))

    # d_a[i], d_b[i], alt_a[i], alt_b[i], T_lcm[i] = find_min_conj_dist(s1l1[idx1], s1l2[idx1], s2l1[idx2], s2l2[idx2])
    d_a[i], d_b[i], alt_a[i], alt_b[i], c, d, T_lcm[i] = find_min_conj_dist_v2(sat_params_s1[idx1], sat_params_s2[idx2])

    t_tol = 10
    ttc_a[i], int_pt_a[i], col_a[i] = time_until_collision_v2(sat_params_s1[idx1], sat_params_s2[idx2], T_lcm[i], 0, t_tol)
    ttc_b[i], int_pt_b[i], col_b[i] = time_until_collision_v2(sat_params_s1[idx1], sat_params_s2[idx2], T_lcm[i], 1, t_tol)

    # min_nodes[i] = min(d_a, d_b)
    # print(f"Analytical Nodal Conjunction Distance ≈ {min_nodes[i]:.2f} km")

check = np.sum(col_a) + np.sum(col_b)
print(f"Total number of collision timesteps: {check}")

t = np.array(t)
plt.figure(figsize = (5,3))
plt.plot(t, d_a, color = 'darkcyan', alpha = 0.6, label = 'Node A')
plt.plot(t, d_b, color = 'firebrick', alpha = 0.6, label = 'Node B')
# plt.axvline(x=t[ttc_a > 0][0], color='darkcyan', linestyle='-', alpha=0.5, label='Collision Node A')
# plt.axvline(x=t[ttc_b > 0][0], color='firebrick', linestyle='-', alpha=0.5, label='Collision Node B')
plt.grid(axis = 'both', color = 'lightgray', linewidth = 0.25)
# plt.xlabel('Date')
# plt.xticks([])
plt.axvline(x = dt.datetime(2024,2,28), color = 'goldenrod', alpha = 0.8)
plt.xlabel('Date')
plt.ylabel('Nodal offset distance [km]')
plt.tight_layout()
plt.show()


plt.figure()
plt.plot(col_a, color='darkcyan', label='Node A Collision')
plt.plot(col_b, color='firebrick', label='Node B Collision')
plt.show()

# plot d_a and d_b (two curves) over 
# t = np.array(t)
plt.figure(figsize = (8,5))
plt.subplot(3,1,1)
plt.plot(t, d_a, color = 'darkcyan', alpha = 0.6, label = 'Node A')
plt.plot(t, d_b, color = 'firebrick', alpha = 0.6, label = 'Node B')
# if col_a[i] == True, put a darkcyan vertical line
for i in range(len(t)):
    if col_a[i]:
        plt.axvline(x=t[i], color='darkcyan', linestyle='-', alpha=0.5)
    if col_b[i]:
        plt.axvline(x=t[i], color='firebrick', linestyle='-', alpha=0.5)
# plt.axvline(x=t[ttc_a > 0][0], color='darkcyan', linestyle='-', alpha=0.5, label='Collision Node A')
# plt.axvline(x=t[ttc_b > 0][0], color='firebrick', linestyle='-', alpha=0.5, label='Collision Node B')
plt.grid(axis = 'both', color = 'gray', linewidth = 0.25)
# plt.xlabel('Date')
# plt.xticks([])
plt.ylabel('Nodal offset distance [km]')
plt.subplot(3,1,2)
plt.plot(t, alt_a-6378.15, color = 'darkcyan', alpha = 0.6, label = 'Node A')
plt.plot(t, alt_b-6378.15, color = 'firebrick', alpha = 0.6, label = 'Node B')
plt.grid(axis = 'both', color = 'gray', linewidth = 0.25)
# plt.xlabel('Date')
# plt.xticks([])
plt.ylabel('Nodal altitude[km]')
plt.subplot(3,1,3)
plt.plot(t, T_lcm, color = 'k')
plt.grid(axis = 'both', color = 'gray', linewidth = 0.25)
plt.xlabel('Date')
plt.ylabel(r'$T_c$')
plt.subplots_adjust(hspace=0.2, wspace=0.1)
# plt.tight_layout()
plt.show()



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

