from skyfield.api import EarthSatellite, load
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import timedelta
import time
from concurrent.futures import ThreadPoolExecutor
import pickle
from concurrent.futures import ProcessPoolExecutor
import os
from tqdm import tqdm
from scipy.optimize import fsolve
from scipy.signal import argrelextrema
import csv
import pandas as pd
from fractions import Fraction
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
# from memory_profiler import profile

earth_radius = 6378.137  # km
mu = 398600.435507  # km^3/s^2

def main():
    # Example TLEs
    # tle1_line1 = "1 25544U 98067A   20264.89222821  .00001264  00000-0  29661-4 0  9991"
    # tle1_line2 = "2 25544  51.6456  21.0985 0005571  41.2232 318.8816 15.49114663246565"

    # tle2_line1 = "1 43205U 18015A   20264.87073778  .00000647  00000-0  18304-4 0  9994"
    # tle2_line2 = "2 43205  53.0551  55.6069 0012501  66.8434 293.3607 15.08856543157085"

    # load txt file with TLEs and read them
    # check to see if sat_params.pkl exists, if not, create it
    years = np.arange(2000, 2001, 1)
    months = np.arange(1, 2, 1)  # Use specific months or range: np.arange(1, 13, 1)
    for year in years: 
        for month in months:
            print(f"Month: {month}, Year: {year}")
            try:
                sat_params = []
                sat_params = pickle.load(open('data/sat_params_annual_prop/only_real/sat_params_' + str(year) + '.pkl', 'rb'))
                print('Loading satellite params from file!')

            except FileNotFoundError:
                tle_file = open('data/batch_tles_annual/' + str(year)+'.txt', 'r')
                tle_lines = tle_file.readlines()

                print('No sat param file found. Loading TLE parameters...')

                tle_l1 = np.zeros(int(len(tle_lines)/3), dtype=object)
                tle_l2 = np.zeros(int(len(tle_lines)/3), dtype=object)
                obj_name = np.zeros(int(len(tle_lines)/3), dtype=object)
                a_range = np.zeros(int(len(tle_lines)/3), dtype=object)
                obj_norad = np.zeros(int(len(tle_lines)/3), dtype=object)
                sat_params = {}
                j = 0
                for i in range(int(len(tle_lines)/3)):
                    obj_name[i] = tle_lines[3*i][2:-1]
                    obj_norad[i] = int(tle_lines[3*i+1][2:7])
                    tle_l1[i] = tle_lines[3*i+1][:-1]
                    tle_l2[i] = tle_lines[3*i+2][:-1]
                

                    # get all the satellite orbit parameters from the TLEs. Store in a dictionary with the satellite norad id as the key
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

                    # filter to only include objects within a certain altitude range
                    if a_range[0] < 2000 + earth_radius and a_range[1] > 300 + earth_radius:
                        epoch_end = datetime(year, month, 14, 0, 0, 0, tzinfo=timezone.utc)
                        new_elements = propagate_keplerian(a, e, incl, raan, arg_per, meanan, epoch, epoch_end)
                        a_prop, e_prop, incl_prop, raan_prop, arg_per_prop, meanan_prop = new_elements
                        a_range_prop = [a_prop * (1-e_prop), a_prop * (1+e_prop)]
                        period_prop = 2 * np.pi * np.sqrt(a_prop ** 3/mu)
                        sat_params[j] = {
                            'a': a_prop,
                            'a_range': a_range_prop,
                            'e': e_prop,
                            'incl': incl_prop,
                            'raan': raan_prop,
                            'arg_per': arg_per_prop,
                            'meanan': meanan_prop,
                            'epoch': epoch_end,
                            'period': period_prop,
                            'name': obj_name[i],
                            'norad': obj_norad[i]
                        }
                        j += 1

                with open('data/sat_params_annual_prop/only_real/sat_params_' + str(year) + '.pkl', 'wb') as f:
                    pickle.dump(sat_params, f)

            # Initialize variables
            n = len(sat_params)
            parallel_compute = True

            # Check to see if data file already exists
            try:
                with open('data/conj_params_annual_prop/only_real/conj_params_' + str(year) + '.pkl', 'rb') as f:
                    dist_a, dist_b, alt_a, alt_b, T_lcm = pickle.load(f)
                    print('Loading conjunction parameters from file!')

            except FileNotFoundError:

                # Initialize variables
                dist_a = np.full((n, len(sat_params)), np.nan, dtype=np.float32)
                dist_b = np.full((n, len(sat_params)), np.nan, dtype=np.float32)
                alt_a = np.full((n, len(sat_params)), np.nan, dtype=np.float32)
                alt_b = np.full((n, len(sat_params)), np.nan, dtype=np.float32)
                T_lcm = np.full((n, len(sat_params)), np.nan, dtype=np.float32)
                t1 = time.time()

                # Ensure the directory exists
                os.makedirs('data', exist_ok=True)

                print('Computing Conjunction Geometries!')

                # Conjunction geometries file can be very large, use a temporary file to store intermediate results
                if parallel_compute == True:
                    max_workers = min(32, os.cpu_count() or 1)
                    chunk_size = 200  # Adjust this for performance vs. memory tradeoff

                    for start in range(0, n, chunk_size):
                        end = min(start + chunk_size, n)
                        print(f"Processing chunk {start} to {end}...")

                        with ProcessPoolExecutor(max_workers=max_workers) as executor:
                            futures = {executor.submit(compute_conjunction_distance_v2, i, sat_params): i for i in range(start, end)}

                            for future in as_completed(futures):
                                idx = futures[future]
                                try:
                                    dist_a[idx, :], dist_b[idx, :], alt_a[idx, :], alt_b[idx, :], T_lcm[idx, :] = future.result()
                                except Exception as e:
                                    print(f"Error at index {idx}: {e}")
                                if idx % 20 == 0:
                                    print(f"{idx}/{n} completed")
                else:
                    results = []
                    for i in range(n):
                        #dist_a[i,:], dist_b[i,:], alt_a[i,:], alt_b[i,:], T_lcm[i,:] = compute_conjunction_distance_v2(i, sat_params)
                        results.append(compute_conjunction_distance_v2(i, sat_params))

                        # Write results every batch_size iterations
                        if (i + 1) % batch_size == 0 or i == n - 1:
                            with open(temp_path, 'ab') as f:
                                pickle.dump(results, f)
                            results = []  # Clear memory
                            print(f"Saved batch {i + 1}/{n} to disk")

                t2 = time.time()
                print(f"Time taken: {t2 - t1} seconds")

                # Save the results to a file
                try:
                    with open('data/conj_params_annual_prop/only_real/conj_params_' + str(year) + '.pkl', 'wb') as f:
                        pickle.dump((dist_a, dist_b, alt_a, alt_b, T_lcm), f)
                        print(f"Conjunction parameters file created!")
                except OSError as e:
                    if e.errno == 28:
                        print("Error: No space left on device. Please free up some space and try again.")
                    else:
                        print(f"An error occurred while saving the pickle file: {e}")
            
            ###############################################################################

            # Define the threshold distance for conjunction
            threshold = 0.2  # kilometers of miss distance
            t_tol = 1 # seconds of miss distance
            LEO_cutoff = 2000 + earth_radius # km

            print(f"Average distance between intersection points 1a and 2a: {np.nanmean(dist_a)}")
            print(f"Maximum distance between intersection points 1a and 2a: {np.nanmax(dist_a)}")
            print(f"Minimum distance between intersection points 1a and 2a: {np.nanmin(dist_a)}")
            print(f"Average distance between intersection points 1b and 2b: {np.nanmean(dist_b)}")
            print(f"Maximum distance between intersection points 1b and 2b: {np.nanmax(dist_b)}")
            print(f"Minimum distance between intersection points 1b and 2b: {np.nanmin(dist_b)}")
            count_dist_a = np.sum(dist_a < threshold)
            count_dist_b = np.sum(dist_b < threshold)
            print(f"Number of items in dist_a less than 0.1: {count_dist_a}")
            print(f"Number of items in dist_b less than 0.1: {count_dist_b}")

            # Find pairs of satellites with distances below the threshold
            print("Number of conjunctions before filtering: ", str(len(sat_params)**2))
            below_threshold_indices = np.where((dist_a < threshold) & (0.2*24*3600 < T_lcm) & (alt_a < LEO_cutoff)) 
            below_threshold_indices_b = np.where((dist_b < threshold) &  (0.2*24*3600 < T_lcm) & (alt_b < LEO_cutoff))
            print("Number of conjunctions after filtering: ", str(len(below_threshold_indices[0]) + len(below_threshold_indices_b[0])))


            # Filter unique pairs (i < j)
            below_threshold_pairs = [(i, j) for i, j in zip(*below_threshold_indices) if i < j]
            print(f"Number of collision pairs at collision point a {len(below_threshold_pairs)}")
            below_threshold_pairs_b = [(i, j) for i, j in zip(*below_threshold_indices_b) if i < j]
            print(f"Number of collision pairs at collision point b {len(below_threshold_pairs_b)}")
                
            # Process each pair to calculate time until collision
            collision_count = 0
            none_count = 0
            intersection_points = []
            collision_times = []
            altitudes = []
            a_or_b = 0
            out_of_range_1 = []
            out_of_range_2 = []
            saved_vars = []

            # Define output file
            output_file = f"data/collision_results_annual_prop/skip_same_name/collision_results_{year}_clean.csv"

            # Determine whether conjunctions occur
            for i, j in tqdm(below_threshold_pairs, desc="Processing Pairs", unit="pair"):
                name1, norad_id1 = sat_params[i]['name'], sat_params[i]['norad']
                name2, norad_id2 = sat_params[j]['name'], sat_params[j]['norad']
                a_range1 = sat_params[i]['a_range']
                a_range2 = sat_params[j]['a_range']
                incl1 = sat_params[i]['incl'] * 180/np.pi
                incl2 = sat_params[j]['incl'] * 180/np.pi

                name1_upper = name1.upper()
                name2_upper = name2.upper()

                # If either name has DEB, R/B, TBA, or OBJECT, don't skip
                if any(keyword in name1_upper or keyword in name2_upper for keyword in ['DEB', 'R/B', 'TBA', 'OBJECT']):
                    pass  # Continue processing this pair
                else:
                    # Extract only words made of letters (exclude any containing digits or symbols)
                    words1 = set(re.findall(r'\b[A-Z]+\b', name1_upper))
                    words2 = set(re.findall(r'\b[A-Z]+\b', name2_upper))

                    # two objects in the same constellation are assumed to be in resonance
                    if words1 & words2:
                        #print(f"Skipping pair {name1} and {name2} due to common words in names {words1 & words2}.")
                        continue  # Skip this pair

                # Calculate expected time until collision
                time_to_collision_a, intersection_point_a, collision = time_until_collision_v2(sat_params[i], sat_params[j], T_lcm[i,j], a_or_b, t_tol)

                if time_to_collision_a is not None:
                    time_to_collision_days = time_to_collision_a/86400
                    collision_times.append(time_to_collision_a)
                    alt = np.linalg.norm(intersection_point_a) - earth_radius
                    intersection_points.append(intersection_point_a)
                    altitudes.append(alt)
                    collision_count = collision_count + 1
                else:
                    # if only collecting results for objects with collisions
                    time_to_collision_days = None
                    alt = None
                    none_count = none_count + 1

                    # # if collecting results for all nodes
                    # time_to_collision_days = None
                    # collision_times.append(time_to_collision_a)
                    # alt = np.linalg.norm(intersection_point_a) - earth_radius
                    # intersection_points.append(intersection_point_a)
                    # altitudes.append(alt)
                    # collision_count = collision_count + 1

                T_lcm_days = T_lcm[i,j]/86400

                # comment out if statement if collecting results for all nodes
                if collision == True: 
                    saved_vars.append([name1, norad_id1, a_range1, name2, norad_id2, a_range2, intersection_point_a, a_or_b, alt, collision, time_to_collision_days, T_lcm_days, 1/T_lcm_days, incl1, incl2])
                
            # repeat for collision point b
            a_or_b = 1
            for i, j in tqdm(below_threshold_pairs_b, desc="Processing Pairs", unit="pair"):
                name1, norad_id1 = sat_params[i]['name'], sat_params[i]['norad']
                name2, norad_id2 = sat_params[j]['name'], sat_params[j]['norad']
                a_range1 = sat_params[i]['a_range']
                a_range2 = sat_params[j]['a_range']
                incl1 = sat_params[i]['incl'] * 180/np.pi
                incl2 = sat_params[j]['incl'] * 180/np.pi

                name1_upper = name1.upper()
                name2_upper = name2.upper()

                # If either name has DEB, R/B, TBA, or OBJECT, don't skip
                if any(keyword in name1_upper or keyword in name2_upper for keyword in ['DEB', 'R/B', 'TBA', 'OBJECT']):
                    pass  # Continue processing this pair
                else:
                    # Extract only words made of letters (exclude any containing digits or symbols)
                    words1 = set(re.findall(r'\b[A-Z]+\b', name1_upper))
                    words2 = set(re.findall(r'\b[A-Z]+\b', name2_upper))

                    if words1 & words2:
                        #print(f"Skipping pair {name1} and {name2} due to common words in names {words1 & words2}.")
                        continue  # Skip this pair

                # Calculate expected time until collision
                time_to_collision_b, intersection_point_b, collision = time_until_collision_v2(sat_params[i], sat_params[j], T_lcm[i,j], a_or_b, t_tol)

                if time_to_collision_b is not None:
                    time_to_collision_days = time_to_collision_b/86400
                    collision_times.append(time_to_collision_b)
                    alt = np.linalg.norm(intersection_point_b) - earth_radius
                    intersection_points.append(intersection_point_b)
                    altitudes.append(alt)
                    collision_count = collision_count + 1
                else:
                    # if only collecting results for objects with collisions
                    time_to_collision_days = None
                    alt = None
                    none_count = none_count + 1

                    # # if collecting results for all nodes
                    # time_to_collision_days = None
                    # collision_times.append(time_to_collision_b)
                    # alt = np.linalg.norm(intersection_point_b) - earth_radius
                    # intersection_points.append(intersection_point_b)
                    # altitudes.append(alt)
                    # collision_count = collision_count + 1

                T_lcm_days = T_lcm[i,j]/86400
                
                # comment out if statement if collecting results for all nodes
                if collision == True: 
                    saved_vars.append([name1, norad_id1, a_range1, name2, norad_id2, a_range2, intersection_point_b, a_or_b, alt, collision, time_to_collision_days, T_lcm_days, 1/T_lcm_days, incl1, incl2])

            # write saved_vars to a csv file
            with open(output_file, mode="w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(["Satellite 1","NORAD ID 1", "Altitude Range [km]", "Satellite 2", "NORAD ID 2", "Altitude Range [km]", "Conjunction Node Position [eci, km]", "Collision Point A or B", "Conjunction Point Altitude [km]", "Will they Collide?", "Time Until Collision [days]", "Synodic Period [days]", "Collision Frequency [1/days]", "Inclination of Sat 1", "Inclination of Sat 2"])
                writer.writerows(saved_vars)

################################################################################################################################
# HELPER FUNCTIONS
################################################################################################################################

# calculate the time until collision for two objects that share a conjunction node
def time_until_collision_v2(sat_params1, sat_params2, T_lcm, a_or_b, t_tol):
        
        m1, ecc1, a1, period1 = sat_params1['meanan'], sat_params1['e'], sat_params1['a'], sat_params1['period']
        m2, ecc2, a2, period2 = sat_params2['meanan'], sat_params2['e'], sat_params2['a'], sat_params2['period']

        #calculate eccentric anomalies of starting points
        e1 = kepler_eqn(m1, ecc1)
        e2 = kepler_eqn(m2, ecc2)

        #use true anomalies to calculate eccentric anomalies of all possible collision points
        intersection_point_1a, intersection_point_1b, intersection_point_2a, intersection_point_2b, f1a, f1b, f2a, f2b = find_orbital_intersection_v2(sat_params1, sat_params2)

        if a_or_b == 0:
            E_col_1 = 2 * np.arctan(np.sqrt((1 - ecc1) / (1 + ecc1)) * np.tan(f1a / 2))
            E_col_2 = 2 * np.arctan(np.sqrt((1 - ecc2) / (1 + ecc2)) * np.tan(f2a / 2))
            intersection_point = intersection_point_1a
        elif a_or_b == 1:
            E_col_1 = 2 * np.arctan(np.sqrt((1 - ecc1) / (1 + ecc1)) * np.tan(f1b / 2))
            E_col_2 = 2 * np.arctan(np.sqrt((1 - ecc2) / (1 + ecc2)) * np.tan(f2b / 2))
            intersection_point = intersection_point_1b

        while e1 > E_col_1:
            E_col_1 = E_col_1 + 2*np.pi

        while e2 > E_col_2:
            E_col_2 = E_col_2 + 2*np.pi

        #time from starting points to collision points
        tau1 = (period1/(2 * np.pi)) * (E_col_1 - e1 - ecc1 * (np.sin(E_col_1) - np.sin(e1)))
        tau2 = (period2/(2 * np.pi)) * (E_col_2 - e2 - ecc2 * (np.sin(E_col_2) - np.sin(e2)))

        max_mult_1 = T_lcm//period1
        max_mult_2 = T_lcm//period2

        t1a = np.arange(max_mult_1) * period1 + tau1
        t2a = np.arange(max_mult_2) * period2 + tau2

        # if any of the differences between any t1a and t2a are less than the tolerance, then there is a collision
        # Compute the pairwise differences
        diff_matrix = np.abs(t1a[:, np.newaxis] - t2a)

        # find instances where the difference is less than the tolerance
        col_idx = np.where(diff_matrix < t_tol)

        if col_idx[0].size == 0:
            time_to_collision = None
            collision = False

        else:
            time_to_collision = t1a[col_idx[0][0]]
            collision = True

        return time_to_collision, intersection_point, collision


# Helper function to parse TLEs and return satellite object
def get_satellite(tle_line1, tle_line2):
    satellite = EarthSatellite(tle_line1, tle_line2, 'Satellite', load.timescale())
    return satellite

# Find the closest rational approximation of T1/T2 within a given tolerance
def closest_rational_approx(T1, T2, tol=1e-3):
    ratio = T1 / T2  # Compute ratio first
    frac = Fraction.from_float(ratio).limit_denominator(int(1/tol))  # Use from_float for floating points
    return frac.denominator * T1  # Compute T_lcm

def kepler_eqn(M, ecc, tol=1e-6, max_iter=100):
    """
    Solve Kepler's equation for the eccentric anomaly E using a for loop.
    
    Parameters:
    M (float): Mean anomaly (radians)
    ecc (float): Eccentricity
    tol (float): Tolerance for the solution
    max_iter (int): Maximum number of iterations
    
    Returns:
    float: Eccentric anomaly (radians)
    """
    E = M  # Initial guess
    for _ in range(max_iter):
        E_new = E - (E - ecc * np.sin(E) - M) / (1 - ecc * np.cos(E))
        if abs(E_new - E) < tol:
            return E_new
        E = E_new
    
    return E  # Return the closest approximation if tolerance is not met

def orbital_plane_normal(inclination, raan):
    """
    Compute the normal to the orbital plane using inclination and RAAN (Right Ascension of Ascending Node).
    """

    # Normal vector in the inertial frame
    normal = np.array([
        np.sin(raan) * np.sin(inclination),
        -np.cos(raan) * np.sin(inclination),
        np.cos(inclination)
    ])
    
    return normal / np.linalg.norm(normal)  # Normalize the vector

def compute_plane_intersection(normal1, normal2):
    """
    Compute the direction vector of the line of intersection between two planes.
    """
    line_direction = np.cross(normal1, normal2)
    cross_magnitude = np.linalg.norm(line_direction)
    #if cross_magnitude < 1e-4:
        #print("satellites are coplanar")

    return line_direction / np.linalg.norm(line_direction)  # Normalize the vector


def find_intersection_point_analytically_v2(sat_params_i, line_direction):
    # find the angle between the line of intersection and the satellite's apse line (line in the direction of eccentricity)
    # convert apse line from perifocal to ECI frame
    apse_line_pf = np.array([1,0,0])
    apse_line = pf_to_eci_v2(apse_line_pf, sat_params_i)

    angle = np.arccos(np.dot(line_direction, apse_line))
    true_anomaly_a = angle
    true_anomaly_b = -angle

    # compute the distance of the orbit when the true anomaly is the angle between the line of intersection and the apse line
    d = sat_params_i['a']*(1-sat_params_i['e']**2)/(1+sat_params_i['e']*np.cos(angle))

    # compute the position vector of the satellite at d distance from the focus of the orbit along the line_direction
    pos_point_a = d * line_direction

    # repeat the same process for the other intersection point (switch the line direction)
    line_direction = -line_direction
    angle = np.arccos(np.dot(line_direction, apse_line))
    d = sat_params_i['a']*(1-sat_params_i['e']**2)/(1+sat_params_i['e']*np.cos(angle))
    pos_point_b = d * line_direction

    return pos_point_a, pos_point_b, true_anomaly_a, true_anomaly_b

def pf_to_eci_v2(pos_pf, sat_params_i):
    # convert from perifocal frame to ECI frame position vector
    i, raan, arg_pe = sat_params_i['incl'], sat_params_i['raan'], sat_params_i['arg_per'] # Arguments in radians
    pos_eci = np.zeros(3)
    pos_eci[0] = (np.cos(raan)*np.cos(arg_pe) - np.sin(raan)*np.sin(arg_pe)*np.cos(i))*pos_pf[0] + (-np.cos(raan)*np.sin(arg_pe) - np.sin(raan)*np.cos(arg_pe)*np.cos(i))*pos_pf[1]
    pos_eci[1] = (np.sin(raan)*np.cos(arg_pe) + np.cos(raan)*np.sin(arg_pe)*np.cos(i))*pos_pf[0] + (-np.sin(raan)*np.sin(arg_pe) + np.cos(raan)*np.cos(arg_pe)*np.cos(i))*pos_pf[1]
    pos_eci[2] = (np.sin(i)*np.sin(arg_pe))*pos_pf[0] + (np.sin(i)*np.cos(arg_pe))*pos_pf[1]
    return pos_eci

def compute_positions_keplerian(satellite, deg_step=1):
    """
    Compute the positions of the satellite over an entire revolution using Keplerian dynamics.

    Parameters:
    satellite (EarthSatellite): The satellite object.
    deg_step (int): The true anomaly step in degrees for position computation.

    Returns:
    np.ndarray: Array of positions in ECI space (km).
    """
    # Extract Keplerian elements from the satellite object
    a = satellite.model.a*earth_radius  # Semi-major axis in km
    e = satellite.model.ecco  # Eccentricity
    i = satellite.model.inclo  # Inclination 
    raan = satellite.model.nodeo  # Right ascension of ascending node
    arg_pe = satellite.model.argpo  # Argument of perigee 
    M0 = satellite.model.mo  # Mean anomaly at epoch 
    epoch = satellite.epoch.utc_datetime()  # Epoch time

    # theta range for one orbital period
    theta_step = [np.radians(deg_step * i) for i in range(int(360 / deg_step))]

    # Compute ranges as a function of theta
    r = a * (1 - e**2) / (1 + e * np.cos(theta_step))
    x = r * np.cos(theta_step)
    y = r * np.sin(theta_step)
    z = np.zeros_like(x)

    # Rotate the position vector from the perifocal to the ECI frame
    # Rotation matrices
    R3_W = np.array([[np.cos(raan), -np.sin(raan), 0],
                    [np.sin(raan), np.cos(raan), 0],
                    [0, 0, 1]])
    R1_i = np.array([[1, 0, 0],
                    [0, np.cos(i), -np.sin(i)],
                    [0, np.sin(i), np.cos(i)]])
    R3_w = np.array([[np.cos(arg_pe), -np.sin(arg_pe), 0],
                    [np.sin(arg_pe), np.cos(arg_pe), 0],
                    [0, 0, 1]])

    # Combined rotation matrix
    R = R3_W @ R1_i @ R3_w

    # Rotate the position vectors
    positions = np.array([R @ np.array([x[i], y[i], z[i]]) for i in range(len(x))])
    return positions

# Function to plot the orbits and intersection points
def plot_orbits_and_intersection(satellite1, satellite2, intersection_line, intersection_point_1a, intersection_point_1b, intersection_point_2a, intersection_point_2b):
    ts = load.timescale()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Get orbit positions for each satellite
    positions1 = compute_positions_keplerian(satellite1)
    positions2 = compute_positions_keplerian(satellite2)

    # Plot the orbits
    ax.plot(positions1[:, 0], positions1[:, 1], positions1[:, 2], label='Satellite 1 Orbit')
    ax.plot(positions2[:, 0], positions2[:, 1], positions2[:, 2], label='Satellite 2 Orbit')

    # Plot the intersection line
    line_scale = 10000  # Scale the line to visualize it
    ax.plot([0, intersection_line[0] * line_scale], 
            [0, intersection_line[1] * line_scale], 
            [0, intersection_line[2] * line_scale], 
            color='black', linestyle = '--', label='Intersection Line')
    
    # Plot in both positive and negative directions
    ax.plot([0, -intersection_line[0] * line_scale],
            [0, -intersection_line[1] * line_scale],
            [0, -intersection_line[2] * line_scale],
            color='black', linestyle='--')

    # Plot the intersection points
    if intersection_point_1a is not None:
        ax.scatter(intersection_point_1a[0], intersection_point_1a[1], intersection_point_1a[2], color='tab:blue', s=30, label='Intersection Point 1a')
    if intersection_point_1b is not None:
        ax.scatter(intersection_point_1b[0], intersection_point_1b[1], intersection_point_1b[2], color='tab:blue', s=30, label='Intersection Point 1b')

    if intersection_point_2a is not None:
        ax.scatter(intersection_point_2a[0], intersection_point_2a[1], intersection_point_2a[2], color='tab:orange', s=30, label='Intersection Point 2a')
    if intersection_point_2b is not None:
        ax.scatter(intersection_point_2b[0], intersection_point_2b[1], intersection_point_2b[2], color='tab:orange', s=30, label='Intersection Point 2b')

    # Set labels and legend
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    # make axes equal
    ax.set_box_aspect([1,1,1])
    ax.legend()
    ax.set_title('Orbital Intersection Visualization')

    plt.show()


def find_orbital_intersection_v2(sat_params1, sat_params2):
    inclination1, raan1, period1 = sat_params1['incl'], sat_params1['raan'], sat_params1['period']
    inclination2, raan2, period2 = sat_params2['incl'], sat_params2['raan'], sat_params2['period']
    
    # Compute normal vectors for each orbital plane
    normal1 = orbital_plane_normal(inclination1, raan1)
    normal2 = orbital_plane_normal(inclination2, raan2)

    # Compute the line of intersection between the two planes
    intersection_line = compute_plane_intersection(normal1, normal2)

    # Find the intersection points on both orbits using the new function
    intersection_point_1a, intersection_point_1b, true_anomaly_1a, true_anomaly_1b = find_intersection_point_analytically_v2(sat_params1, intersection_line)
    intersection_point_2a, intersection_point_2b, true_anomaly_2a, true_anomaly_2b = find_intersection_point_analytically_v2(sat_params2, intersection_line)

    return intersection_point_1a, intersection_point_1b, intersection_point_2a, intersection_point_2b, true_anomaly_1a, true_anomaly_1b, true_anomaly_2a, true_anomaly_2b

def find_min_conj_dist_v2(sat_params1, sat_params2):

    intersection_point_1a, intersection_point_1b, intersection_point_2a, intersection_point_2b, true_anomaly_1a, true_anomaly_1b, true_anomaly_2a, true_anomaly_2b = find_orbital_intersection_v2(sat_params1, sat_params2)
    period1 = sat_params1['period']
    period2 = sat_params2['period']

    dist_a = np.linalg.norm(intersection_point_1a - intersection_point_2a)
    dist_b = np.linalg.norm(intersection_point_1b - intersection_point_2b)
    alt_a = np.mean([np.linalg.norm(intersection_point_1a), np.linalg.norm(intersection_point_2a)])
    alt_b = np.mean([np.linalg.norm(intersection_point_1b), np.linalg.norm(intersection_point_2b)])

    T_lcm = closest_rational_approx(period1, period2)

    return dist_a, dist_b, alt_a, alt_b, intersection_point_1a, intersection_point_1b, T_lcm


# Define the function to compute the minimum conjunction distance for a single value of j
def compute_conjunction_distance_v2(i, sat_params):
    dist_a = np.full(len(sat_params), np.nan)
    dist_b = np.full(len(sat_params), np.nan)
    alt_a = np.full(len(sat_params), np.nan)
    alt_b = np.full(len(sat_params), np.nan)
    pos_a = np.full((len(sat_params), 3),  np.nan)
    pos_b = np.full((len(sat_params), 3), np.nan)
    T_lcm = np.full(len(sat_params), np.nan)

    for j in range(len(sat_params)):
        if i >= j:
            continue
        if sat_params[j]['a_range'][0] > (sat_params[i]['a_range'][1] + 0.1) or sat_params[j]['a_range'][1] < (sat_params[i]['a_range'][0] - 0.1):
            continue # Skip any cases where the orbits do not have any likelihood of conjunction
        # Find the intersection points and plot the result
        dist_a[j], dist_b[j], alt_a[j], alt_b[j], pos_a[j,:], pos_b[j,:], T_lcm[j] = find_min_conj_dist_v2(sat_params[i], sat_params[j])

    return dist_a, dist_b, alt_a, alt_b, T_lcm

def propagate_keplerian(a, e, i, RAAN, omega, M0, t_0, t):
    """
    Propagates the orbital elements using Kepler's equation.

    Args:
        input_oe (numpy.ndarray): Initial orbital elements [a, e, i, RAAN, omega, M].
        param (dict): Dictionary containing parameters like mu, t, and t_0.

    Returns:
        numpy.ndarray: Propagated orbital elements [a, e, i, RAAN, omega, M].
    """

    # Mean motion
    n = np.sqrt(mu) * a**(-3/2)

    dt = (t - t_0).total_seconds()  # Time difference in seconds

    M_no_wrap = M0 + n * (dt)  # Mean anomaly at time t
    M = np.mod(M_no_wrap, 2 * np.pi)  # Wrap to [0, 2*pi]

    return a, e, i, RAAN, omega, M


if __name__ == '__main__':
    main()