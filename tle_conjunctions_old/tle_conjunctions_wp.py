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


def main():
    # Example TLEs
    # tle1_line1 = "1 25544U 98067A   20264.89222821  .00001264  00000-0  29661-4 0  9991"
    # tle1_line2 = "2 25544  51.6456  21.0985 0005571  41.2232 318.8816 15.49114663246565"

    # tle2_line1 = "1 43205U 18015A   20264.87073778  .00000647  00000-0  18304-4 0  9994"
    # tle2_line2 = "2 43205  53.0551  55.6069 0012501  66.8434 293.3607 15.08856543157085"

    # load txt file with TLEs and read them
    tle_file = open('3le_leo_092424.txt', 'r')
    tle_lines = tle_file.readlines()

    print('Loading TLE parameters...')

    tle_l1 = np.zeros(int(len(tle_lines)/3), dtype=object)
    tle_l2 = np.zeros(int(len(tle_lines)/3), dtype=object)
    obj_name = np.zeros(int(len(tle_lines)/3), dtype=object)
    a_range = np.zeros(int(len(tle_lines)/3), dtype=object)
    for i in range(int(len(tle_lines)/3)):
        obj_name[i] = tle_lines[3*i][2:-1]
        tle_l1[i] = tle_lines[3*i+1][:-1]
        tle_l2[i] = tle_lines[3*i+2][:-1]

    # Initialize variables
    n = 10 #len(tle_l1)
    parallel_compute = False

    # Check to see if data file already exists
    try:
        with open('data/conj_params_' + str(n) + '.pkl', 'rb') as f:
            dist_a, dist_b, alt_a, alt_b, T_lcm = pickle.load(f)
    except FileNotFoundError:
        a_range = np.zeros((len(tle_l1), 2))  # Assuming a_range is a 2D array

        # Compute the altitude range
        for i in range(len(tle_l1)):
            satellite = get_satellite(tle_l1[i], tle_l2[i])
            a = satellite.model.a  # Semi-major axis in Earth radii
            e = satellite.model.ecco  # Eccentricity
            a_range[i] = [a * (1 - e), a * (1 + e)]
        
        # a_range = a_range*6378.15 # convert to km
        # Initialize variables
        dist_a = np.full((n, len(tle_l1)), np.nan)
        dist_b = np.full((n, len(tle_l1)), np.nan)
        alt_a = np.full((n, len(tle_l1)), np.nan)
        alt_b = np.full((n, len(tle_l1)), np.nan)
        T_lcm = np.full((n, len(tle_l1)), np.nan)
        t1 = time.time()

        # Ensure the directory exists
        os.makedirs('data', exist_ok=True)

        # Use ProcessPoolExecutor to parallelize the loop over i
        if parallel_compute == True:
            with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
                futures = [executor.submit(compute_conjunction_distance, i, tle_l1, tle_l2, a_range) for i in range(n)]
                for idx, future in enumerate(futures):
                    dist_a[idx, :], dist_b[idx, :], alt_a[idx, :], alt_b[idx, :], T_lcm[idx, :] = future.result()
                    print(idx / n)
        else:
            for i in range(n):
                dist_a[i,:], dist_b[i,:], alt_a[i,:], alt_b[i,:], T_lcm[i,:] = compute_conjunction_distance(i, tle_l1, tle_l2, a_range)
                dist_a_close = dist_a[i,:][dist_a[i,:]<1]
                dist_b_close = dist_b[i,:][dist_b[i,:]<1]
                alt_a_close = alt_a[i,:][dist_a[i,:]<1]
                alt_b_close = alt_b[i,:][dist_b[i,:]<1]

                # plot alt_a and alt_b with horizontal lines marking the a_range
                # plt.figure()
                # plt.plot(alt_a_close-6378.15, 'b.')
                # plt.plot(alt_b_close-6378.15, 'g.')
                # plt.axhline(a_range[i,0]*6378.15-6378.15)
                # plt.axhline(a_range[i,1]*6378.15-6378.15)
                # plt.show()

                # plt.figure()
                # plt.plot(T_lcm[i,:][dist_a[i,:]<1]/(24*3600), '.')
                # plt.ylabel('T_lcm (days)')
                # plt.show()
                
                print(i/n)

        plt.figure()
        plt.imshow(T_lcm/(3600*24), aspect = 'auto')
        plt.title('T_lcm, days')
        plt.colorbar()
        plt.show()

        t2 = time.time()
        print(f"Time taken: {t2 - t1} seconds")

        # Save the results to a file
        try:
            with open('data/conj_params_' + str(n) + '.pkl', 'wb') as f:
                pickle.dump((dist_a, dist_b, alt_a, alt_b, T_lcm), f)
                print(f"File 'data/conj_params_{n}.pkl' created!")
        except OSError as e:
            if e.errno == 28:
                print("Error: No space left on device. Please free up some space and try again.")
            else:
                print(f"An error occurred while saving the pickle file: {e}")


    # plot a heatmap of alt_a vs dist_a and alt_b vs dist_b (all in one plot)
    # plt.figure()
    # plt.scatter(alt_a.flatten()-6378.15, dist_a.flatten(), color='tab:blue', label='a')
    # plt.scatter(alt_b.flatten()-6378.15, dist_b.flatten(), color='tab:blue', label='b')
    # plt.xlabel('Conjunction altitude (km)')
    # plt.ylabel('Closest approach distance (km)')
    # # plt.legend()
    # plt.show()

    # repeat the plot above, but make it a heatmap with number of points within each grid cell
    # make a 2D histogram of alt_a vs dist_a and alt_b vs dist_b
        # Filter out NaN values
    # valid_indices_a = ~np.isnan(alt_a.flatten()) & ~np.isnan(dist_a.flatten())
    # valid_indices_b = ~np.isnan(alt_b.flatten()) & ~np.isnan(dist_b.flatten())

    # alt_a_valid = alt_a.flatten()[valid_indices_a] - 6378.15
    # dist_a_valid = dist_a.flatten()[valid_indices_a]
    # alt_b_valid = alt_b.flatten()[valid_indices_b] - 6378.15
    # dist_b_valid = dist_b.flatten()[valid_indices_b]
    # T_lcm_valid = T_lcm.flatten()[valid_indices_a] #use a or b?

    # # make a combined dist and alt array
    # alt_valid = np.concatenate((alt_a_valid, alt_b_valid))
    # dist_valid = np.concatenate((dist_a_valid, dist_b_valid))

    ###############################################################################
    # After we've identified satellite conjunctions and T_lcm, we need to know if the satellites actually conjoin within T_lcm. 
    # Define the threshold distance for conjunction
    threshold = 0.1  # kilometers

    # Find pairs of satellites with distances below the threshold
    below_threshold_indices = np.where((dist_a < threshold) & (T_lcm < 7776000))
    below_threshold_indices_b = np.where((dist_b < threshold) & (T_lcm < 7776000))

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

    # Define output file
    output_file = "collision_results.csv"

    # Open file for writing
    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)
    
    # Write header
        writer.writerow(["Satellite 1","NORAD ID 1", "Satellite 2", "NORAD ID 2", "Conjunction Node Position [eci, km]", "Collision Point A or B", "Conjunction Point Altitude [km]", "Will they Collide?", "Time Until Collision [days]", "Synodic Period [days]", "Collision Frequency [1/days]"])

        for i, j in tqdm(below_threshold_pairs, desc="Processing Pairs", unit="pair"):
            tle1_line1, tle1_line2 = tle_l1[i], tle_l2[i]
            tle2_line1, tle2_line2 = tle_l1[j], tle_l2[j]
            name1 = obj_name[i]
            name2 = obj_name[j]
            norad_id1 = tle1_line2.split()[1]
            norad_id2 = tle2_line2.split()[1]

            # Calculate expected time until collision
            time_to_collision_a, intersection_point, T_lcm, collision = time_until_collision(tle1_line1, tle1_line2, tle2_line1, tle2_line2, a_or_b)
            alt = np.linalg.norm(intersection_point) - 6378.15

            if time_to_collision_a is not None:
                time_to_collision_days = time_to_collision_a/86400
            else:
                time_to_collision_days = time_to_collision_a
            T_lcm_days = T_lcm/86400

            if time_to_collision_a is not None:
                collision_times.append(time_to_collision_a)
                intersection_points.append(intersection_point)
                altitudes.append(alt)

            if time_to_collision_a != None:
                collision_count = collision_count + 1
            else:
                none_count = none_count + 1

            writer.writerow([name1, norad_id1, name2, norad_id2, intersection_point, a_or_b, alt, collision, time_to_collision_days, T_lcm_days, 1/T_lcm_days])

        for i, j in tqdm(below_threshold_pairs_b, desc="Processing Pairs", unit="pair"):
            tle1_line1, tle1_line2 = tle_l1[i], tle_l2[i]
            tle2_line1, tle2_line2 = tle_l1[j], tle_l2[j]
            name1 = obj_name[i]
            name2 = obj_name[j]
            norad_id1 = tle1_line2.split()[1]
            norad_id2 = tle2_line2.split()[1]
            a_or_b = 1

                # Calculate expected time until collision
            time_to_collision_a, intersection_point, T_lcm, collision = time_until_collision(tle1_line1, tle1_line2, tle2_line1, tle2_line2, a_or_b)
            alt = np.linalg.norm(intersection_point) - 6378.15

            if time_to_collision_a is not None:
                time_to_collision_days = time_to_collision_a/86400
            else:
                time_to_collision_days = time_to_collision_a
            T_lcm_days = T_lcm/86400

            if time_to_collision_a is not None:
                collision_times.append(time_to_collision_a)
                intersection_points.append(intersection_point)
                altitudes.append(alt)

                # Print or log the results
            if time_to_collision_a != None:
                    #print(f"Satellite pair ({i}, {j}):")
                    #print(time_to_collision_a)
                collision_count = collision_count + 1
            else:
                none_count = none_count + 1

            writer.writerow([name1, norad_id1, name2, norad_id2, intersection_point, a_or_b, alt, collision, time_to_collision_days, T_lcm_days, 1/T_lcm_days])

        print(f"none count {none_count}")
        print(f"collision count {collision_count}")

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # Extract x, y, z coordinates from the points
        x_coords = [point[0] for point in intersection_points]
        y_coords = [point[1] for point in intersection_points]
        z_coords = [point[2] for point in intersection_points]

        # Plot the points, with colors corresponding to collision times
        collision_times_days = [time / 86400 for time in collision_times]
        scatter = ax.scatter(x_coords, y_coords, z_coords, c=collision_times_days, cmap='viridis', marker='o')
        
        # Add color bar to indicate time-to-collision values
        cbar = plt.colorbar(scatter, ax=ax, label='Expected Time to Collision (days)')
        
        # Add axis labels
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_zlabel('Z Coordinate')
        ax.set_title('Conjunction Node Positions')
        
        # Show the plot
        plt.show()
        plt.scatter(altitudes, collision_times_days)
        plt.xlabel('Altitude of Conjunction Node')
        plt.ylabel('Expected Time To Collision')
        plt.show()
    #print(f"  Time until collision at point A: {time_to_collision_a:.2f} seconds")



# Function to calculate time until collision
def time_until_collision(tle1_line1, tle1_line2, tle2_line1, tle2_line2, a_or_b):

        satellite1 = get_satellite(tle1_line1, tle1_line2)
        satellite2 = get_satellite(tle2_line1, tle2_line2)
        M1 = np.deg2rad(satellite1.model.mo)
        M2 = np.deg2rad(satellite2.model.mo)
        ecc1 = satellite1.model.ecco
        ecc2 = satellite2.model.ecco
        a1 = satellite1.model.a * 6378
        a2 = satellite2.model.a * 6378
        mu = 398600.435507 #[km^3/s^2]
        T1 = 2 * np.pi * np.sqrt(a1 ** 3/mu)
        T2 = 2 * np.pi * np.sqrt(a2 ** 3/mu)

        #calculate eccentric anomalies of starting points
        E1 = M1 + ecc1 * np.sin(M1)
        E2 = M2 + ecc2 * np.sin(M2)

        #use true anomalies to calculate eccentric anomalies of all possible collision points
        intersection_point, f1a, f1b, f2a, f2b = find_orbital_intersection(tle1_line1, tle1_line2, tle2_line1, tle2_line2)

        if a_or_b == 0:
            E_col_1 = 2 * np.arctan(np.sqrt((1 - ecc1) / (1 + ecc1)) * np.tan(f1a / 2))
            E_col_2 = 2 * np.arctan(np.sqrt((1 - ecc2) / (1 + ecc2)) * np.tan(f2a / 2))
        elif a_or_b == 1:
            E_col_1 = 2 * np.arctan(np.sqrt((1 - ecc1) / (1 + ecc1)) * np.tan(f1b / 2))
            E_col_2 = 2 * np.arctan(np.sqrt((1 - ecc2) / (1 + ecc2)) * np.tan(f2b / 2))

        while E1 > E_col_1:
            E_col_1 = E_col_1 + 2*np.pi

        while E2 > E_col_2:
            E_col_2 = E_col_2 + 2*np.pi

        #time from starting points to collision points
        #if a_or_b == 0:
        tau1 = (T1/(2 * np.pi)) * (E_col_1 - E1 - ecc1 * (np.sin(E_col_1) - np.sin(E1)))
        tau2 = (T2/(2 * np.pi)) * (E_col_2 - E2 - ecc2 * (np.sin(E_col_2) - np.sin(E2)))
        #if a_or_b == 1:
        #    tau1 = (T1/(2 * np.pi)) * (E_col_1 - E1 - ecc1 * (np.sin(E_col_1) - np.sin(E1)))
        #    tau2 = (T2/(2 * np.pi)) * (E_col_2 - E2 - ecc2 * (np.sin(E_col_2) - np.sin(E2)))

        #ETTC for collision point a
        T_lcm = 1 / abs((1 / T1) - (1 / T2))

        #if T_lcm > 7776000:
        #    time_to_collision = None
        #    collision = False
        #    return time_to_collision, intersection_point, T_lcm, collision

        tvec = np.arange(2000) 
        t1a = tau1 + T1*tvec
        t2a = tau2 + T2*tvec
        t1a = t1a[t1a <= 7786000]
        t2a = t2a[t2a <= 7786000]
        #print(f"T_lcm:, {T_lcm}, Length of t1a: {len(t1a)}, Last values: {t1a[-4]}, {t1a[-3]}, {t1a[-2]}, {t1a[-1]}")
        #print(f"T_lcm:, {T_lcm}, Length of t2a: {len(t2a)}, Last values: {t2a[-4]}, {t2a[-3]}, {t2a[-2]}, {t2a[-1]}")

        intersection_within_threshold = []
        i, j = 0, 0
        time_to_collision = 0

        ''' Entropy and Variance section:
        distance_between_sats = []
        vec = tau1a + T1*(np.arange(1))
        for i in range (len(vec)):
            sat_2_position = position_of_sat(satellite2, vec[i])
            distance_between_sats_temp = np.linalg.norm(sat_2_position - intersection_point)
            distance_between_sats.append(distance_between_sats_temp)

        distance_between_sats = np.array(distance_between_sats)
        min_indices = argrelextrema(distance_between_sats, np.less)[0]
        # Extract minima y-values
        y_mins = distance_between_sats[min_indices]
        # Compute variance of local minima
        variance_minima = np.var(y_mins)
        if y_mins.any():
            norm_variance = variance_minima/(min(y_mins))
        else:
            norm_variance = variance_minima
        '''

        while i < len(t1a) and j < len(t2a):
            if np.abs(t1a[i] - t2a[j]) <= 0.1: #0.1 seconds, collision threshold
                # If within threshold, add to the result
                break
            elif t1a[i] < t2a[j]:
                i += 1
            else:
                j += 1

        if i < len(t1a) and t1a[i] < 7776000:
            time_to_collision = t1a[i]
            collision = True
        else:
            time_to_collision = None
            collision = False

        #plt.scatter(vec, distance_between_sats)
        #plt.xlabel("Time from Epoch (seconds)")
        #plt.ylabel("Distance between satellites when Sat A is at conjunction point [km]")
        #plt.show()
        #print(f"time to collision {time_to_collision}")
        return time_to_collision, intersection_point, T_lcm, collision

# Helper function to parse TLEs and return satellite object
def get_satellite(tle_line1, tle_line2):
    satellite = EarthSatellite(tle_line1, tle_line2, 'Satellite', load.timescale())
    return satellite

# Function to compute the normal vector to an orbital plane
def orbital_plane_normal(inclination, raan):
    """
    Compute the normal to the orbital plane using inclination and RAAN (Right Ascension of Ascending Node).
    """
    # inclination_rad = np.radians(inclination)
    # raan_rad = np.radians(raan)

    # Normal vector in the inertial frame
    normal = np.array([
        np.sin(raan) * np.sin(inclination),
        -np.cos(raan) * np.sin(inclination),
        np.cos(inclination)
    ])
    
    return normal / np.linalg.norm(normal)  # Normalize the vector

# Function to compute the intersection line of two orbital planes
def compute_plane_intersection(normal1, normal2):
    """
    Compute the direction vector of the line of intersection between two planes.
    """
    line_direction = np.cross(normal1, normal2)
    cross_magnitude = np.linalg.norm(line_direction)
    #if cross_magnitude < 1e-4:
        #print("satellites are coplanar")

    return line_direction / np.linalg.norm(line_direction)  # Normalize the vector

def find_intersection_point_analytically(satellite, line_direction):
    # find the angle between the line of intersection and the satellite's apse line (line in the direction of eccentricity)
    # apse_line_pf = np.array([np.cos(satellite.model.argpo), np.sin(satellite.model.argpo), 0])
    # convert apse line from perifocal to ECI frame
    apse_line_pf = np.array([1,0,0])
    apse_line = pf_to_eci(apse_line_pf, satellite)

    angle = np.arccos(np.dot(line_direction, apse_line))
    true_anomaly_a = angle
    true_anomaly_b = -angle
    # compute the distance of the orbit when the true anomaly is the angle between the line of intersection and the apse line
    d = (satellite.model.a)*6378.15 * (1 - satellite.model.ecco ** 2) / (1 + satellite.model.ecco * np.cos(angle))
    # compute the position vector of the satellite at d distance from the focus of the orbit along the line_direction
    pos_point_a = d * line_direction

    # repeat the same process for the other intersection point (switch the line direction)
    line_direction = -line_direction
    angle = np.arccos(np.dot(line_direction, apse_line))
    d = (satellite.model.a)*6378.15 * (1 - satellite.model.ecco ** 2) / (1 + satellite.model.ecco * np.cos(angle))
    pos_point_b = d * line_direction

    return pos_point_a, pos_point_b, true_anomaly_a, true_anomaly_b

def pf_to_eci(pos_pf, satellite):
    # convert from perifocal frame to ECI frame position vector
    # Extract Keplerian elements
    i = satellite.model.inclo  # Inclination in radians
    raan = satellite.model.nodeo  # Right ascension of ascending node in radians
    arg_pe = satellite.model.argpo  # Argument of perigee in radians
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
    a = satellite.model.a*6378.15  # Semi-major axis in km
    e = satellite.model.ecco  # Eccentricity
    i = satellite.model.inclo  # Inclination 
    raan = satellite.model.nodeo  # Right ascension of ascending node
    arg_pe = satellite.model.argpo  # Argument of perigee 
    M0 = satellite.model.mo  # Mean anomaly at epoch 
    epoch = satellite.epoch.utc_datetime()  # Epoch time

    # Orbital period (seconds)
    # mu = 398600.4418  # Gravitational parameter for Earth (km^3/s^2)
    # orbital_period = 2 * np.pi * np.sqrt(a**3 / mu)

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


def find_min_conj_dist(tle1_line1, tle1_line2, tle2_line1, tle2_line2):
    # Load the satellites from TLEs
    satellite1 = get_satellite(tle1_line1, tle1_line2)
    satellite2 = get_satellite(tle2_line1, tle2_line2)

    # Get orbital parameters (inclination and RAAN) from TLEs
    inclination1, raan1 = satellite1.model.inclo, satellite1.model.nodeo
    inclination2, raan2 = satellite2.model.inclo, satellite2.model.nodeo

    # Compute normal vectors for each orbital plane
    normal1 = orbital_plane_normal(inclination1, raan1)
    normal2 = orbital_plane_normal(inclination2, raan2)

    # Compute the line of intersection between the two planes
    intersection_line = compute_plane_intersection(normal1, normal2)
    # print(f"Intersection line direction: {intersection_line}")

    # Find the intersection points on both orbits
    # intersection_point_1 = find_intersection_point(satellite1, intersection_line)
    # intersection_point_2 = find_intersection_point(satellite2, intersection_line)

    # Find the intersection points on both orbits using the new function
    intersection_point_1a, intersection_point_1b, true_anomaly_1a, true_anomaly_1b = find_intersection_point_analytically(satellite1, intersection_line)
    intersection_point_2a, intersection_point_2b, true_anomaly_2a, true_anomaly_2b = find_intersection_point_analytically(satellite2, intersection_line)

    # min_dist = min(np.linalg.norm(intersection_point_1a - intersection_point_2a), np.linalg.norm(intersection_point_1b - intersection_point_2b))
    # arg = np.argmin([np.linalg.norm(intersection_point_1a - intersection_point_2a), np.linalg.norm(intersection_point_1b - intersection_point_2b)])
    # conj_alt = np.linalg.norm(intersection_point_1a) if arg == 0 else np.linalg.norm(intersection_point_1b)

    dist_a = np.linalg.norm(intersection_point_1a - intersection_point_2a)
    dist_b = np.linalg.norm(intersection_point_1b - intersection_point_2b)
    alt_a = np.mean([np.linalg.norm(intersection_point_1a), np.linalg.norm(intersection_point_2a)])
    alt_b = np.mean([np.linalg.norm(intersection_point_1b), np.linalg.norm(intersection_point_2b)])
    # print the distance between intersection points 1a and 2a, 1b and 2b

    a1 = satellite1.model.a * 6378
    a2 = satellite2.model.a * 6378
    mu = 398600.435507 #[km^3/s^2]
    T1 = 2 * np.pi * np.sqrt(a1 ** 3/mu)
    T2 = 2 * np.pi * np.sqrt(a2 ** 3/mu)
    # T_lcm = 1 / abs((1 / T1) - (1 / T2))
    T_lcm = compute_T_lcm(T1,T2)



    return dist_a, dist_b, alt_a, alt_b, T_lcm

def compute_T_lcm(T1, T2, tol=1e-3):
    """Find the closest rational approximation of T1/T2 within a given tolerance.
    
    Args:
        T1 (float): First period.
        T2 (float): Second period.
        tol (float): Allowed tolerance in the approximation.
    
    Returns:
        tuple: (p, q) such that T1/T2 ≈ p/q within tolerance.
    """
    ratio = T1 / T2
    frac = Fraction.from_float(ratio).limit_denominator(int(1/tol))  # Finds best rational approximation
    T_lcm = frac.denominator*T1 # CHECK TO MAKE SURE THAT DENOM ALWAYS PAIRS WITH T1
    return T_lcm

# Define the function to compute the minimum conjunction distance for a single value of j
def compute_conjunction_distance(i, tle_l1, tle_l2, a_range):
    dist_a = np.full(len(tle_l1), np.nan)
    dist_b = np.full(len(tle_l1), np.nan)
    alt_a = np.full(len(tle_l1), np.nan)
    alt_b = np.full(len(tle_l1), np.nan)
    T_lcm = np.full(len(tle_l1), np.nan)

    for j in range(len(tle_l1)):
        if i >= j:
            continue
        if a_range[j][0] > a_range[i][1] or a_range[j][1] < a_range[i][0]: 
            continue # Skip any cases where the orbits do not have any likelihood of conjunction
        tle1_line1 = tle_l1[i]
        tle1_line2 = tle_l2[i]
        tle2_line1 = tle_l1[j]
        tle2_line2 = tle_l2[j]
        # Find the intersection points and plot the result
        dist_a[j], dist_b[j], alt_a[j], alt_b[j], T_lcm[j] = find_min_conj_dist(tle1_line1, tle1_line2, tle2_line1, tle2_line2)
    print(i)
    return dist_a, dist_b, alt_a, alt_b, T_lcm

def position_of_sat(satellite, delta_t):
    distance = 0
    n = satellite.model.no_kozai/60 #mean motion in radians/min
    e = satellite.model.ecco
    a = satellite.model.a*6378.15
    i = satellite.model.inclo
    arg_pe = satellite.model.argpo
    raan = satellite.model.nodeo 
    M0 = satellite.model.mo
    M = M0 + n*delta_t
    def kepler_eq(E, M, e):
        return E - e * np.sin(E) - M
    E_t = fsolve(kepler_eq, M, args=(M, e))[0]  # Solve for Eccentric Anomaly
    nu_t = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E_t / 2),
                      np.sqrt(1 - e) * np.cos(E_t / 2))
    r_t = a * (1 - e ** 2) / (1 + e * np.cos(nu_t))
    r_perifocal = np.array([r_t * np.cos(nu_t), r_t * np.sin(nu_t), 0])

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
    r_eci = R3_W @ R1_i @ R3_w @ r_perifocal
    return r_eci

if __name__ == '__main__':
    main()