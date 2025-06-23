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
    years = np.arange(2019, 2025, 1)
    months = np.arange(1, 13, 1)  # Use specific months or range: np.arange(1, 13, 1)
    for year in years: 
        for month in months:
            print(f"Month: {month}, Year: {year}")
            try:
                sat_params = []
                #sat_params = pickle.load(open('data/sat_params_annual_prop/sat_params_' + str(year) + '_no_artifical_tles.pkl', 'rb'))
                #sat_params = pickle.load(open(f"data/sat_params_monthly_prop/starlink_artificial/sat_params_{month:02d}{year}.pkl", "rb"))
                #sat_params = pickle.load(open(f"data/sat_params_monthly_prop/starlink_only.pkl", "rb"))
                sat_params = pickle.load(open(f"data/sat_params_monthly_prop/only_real/sat_params_{month:02d}{year}.pkl", "rb"))
                print('Loading satellite params from file!')

            except FileNotFoundError:
                #tle_file = open('data/batch_tles_annual/' + str(year)+'_no_artificial_tles.txt', 'r')
                #tle_file = open(f"data/batch_tles_monthly/starlink_artificial/{month:02d}{year}.txt", "r")
                #tle_file = open(f"data/batch_tles_monthly/filtered_starlink_tles.txt", "r")
                tle_file = open(f"data/batch_tles_monthly/only_real/{month:02d}{year}.txt", "r")
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

                    if a_range[0] < 2000 + earth_radius and a_range[1] > 300 + earth_radius:
                        if obj_norad[i] < 90000:
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
                        else:
                            sat_params[j] = {
                                'a': a,
                                'a_range': a_range,
                                'e': e,
                                'incl': incl,
                                'raan': raan,
                                'arg_per': arg_per,
                                'meanan': meanan,
                                'epoch': epoch,
                                'period': period,
                                'name': obj_name[i],
                                'norad': obj_norad[i]
                            }
                        j += 1

                    #print(i/len(tle_lines)*3)

                '''
                # save sat_params to a file
                print(type(sat_params))
                print(len(sat_params))
                print(sat_params[0])
                filtered_sat_params = {}

                filtered_sat_params = {}

                for key, value in sat_params.items():  # Iterate through dictionary items
                    if value['a'] >= (2000 + earth_radius):  
                        print(value)
                    else:
                        filtered_sat_params[key] = value  # Keep valid satellites

                sat_params = filtered_sat_params  # Replace with filtered dictionary
                '''

                print(len(tle_lines)/3)
                print(f"sat params length {len(sat_params)}")

                #with open('data/sat_params_annual_prop/sat_params_' + str(year) + '_no_artificial_tles.pkl', 'wb') as f:
                #with open(f"data/sat_params_monthly_prop/starlink_artificial/sat_params_{month:02d}{year}.pkl", "wb") as f:
                with open(f"data/sat_params_monthly_prop/only_real/sat_params_{month:02d}{year}.pkl", "wb") as f:
                    pickle.dump(sat_params, f)


            # Initialize variables
            n = len(sat_params)
            parallel_compute = True

            # Check to see if data file already exists
            try:
                #with open('data/conj_params_annual_prop/conj_params_' + str(year) + 'no_artificial_tles.pkl', 'rb') as f:
                #with open(f"data/conj_params_monthly_prop/starlink_artificial/conj_params_{month:02d}{year}.pkl", "rb") as f:
                with open(f"data/conj_params_monthly_prop/only_real/conj_params_{month:02d}{year}.pkl", "rb") as f:
                #with open(f"data/conj_params_monthly_prop/starlink_only.pkl", "rb") as f:
                    dist_a, dist_b, alt_a, alt_b, T_lcm = pickle.load(f)
                    print('Loading conjunction parameters from file!')

            except FileNotFoundError:

                # Initialize variables
                dist_a = np.full((n, len(sat_params)), np.nan)
                dist_b = np.full((n, len(sat_params)), np.nan)
                alt_a = np.full((n, len(sat_params)), np.nan)
                alt_b = np.full((n, len(sat_params)), np.nan)
                # pos_a = np.full((n, len(sat_params), 3), np.nan)
                # pos_b = np.full((n, len(sat_params), 3), np.nan)
                T_lcm = np.full((n, len(sat_params)), np.nan)
                t1 = time.time()

                # Ensure the directory exists
                os.makedirs('data', exist_ok=True)

                print('Computing Conjunction Geometries!')

                #batch_size = 1000
                #temp_path = f'data/conj_params_with_debris_{n}_temp.pkl'

            #with open(temp_path, 'wb') as f:
                #    pickle.dump([], f)  # Initialize with an empty list

                # compute close approach distances for every object pair (don't repeat object pairs in reverse order)
                if parallel_compute == True:
                    max_workers = os.cpu_count()
                    with ProcessPoolExecutor(max_workers=max_workers) as executor:
                        futures = [executor.submit(compute_conjunction_distance_v2, i, sat_params) for i in range(n)]

                        results = []
                        for idx, future in enumerate(futures):
                            #results.append(future.result())
                            """
                            if (idx + 1) % batch_size == 0 or idx == n - 1:
                                with open(temp_path, 'ab') as f:
                                    pickle.dump(results, f)
                                results = []  # Clear memory
                                print(f"Saved batch {idx + 1}/{n} to disk")
                            """
                            dist_a[idx, :], dist_b[idx, :], alt_a[idx, :], alt_b[idx, :], T_lcm[idx, :] = future.result()
                            print(idx / n)
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
                        
                        # dist_a_close = dist_a[i,:][dist_a[i,:]<1]
                        # dist_b_close = dist_b[i,:][dist_b[i,:]<1]
                        # alt_a_close = alt_a[i,:][dist_a[i,:]<1]
                        # alt_b_close = alt_b[i,:][dist_b[i,:]<1]

                        # # plot alt_a and alt_b with horizontal lines marking the a_range
                        # plt.figure()
                        # plt.plot(alt_a_close-6378.15, 'b.')
                        # plt.plot(alt_b_close-6378.15, 'g.')
                        # plt.axhline(sat_params[i]['a_range'][0]-6378.15)
                        # plt.axhline(sat_params[i]['a_range'][1]-6378.15)
                        # plt.show()

                        # plt.figure()
                        # plt.plot(T_lcm[i,:][dist_a[i,:]<1]/(24*3600), '.')
                        # plt.ylabel('T_lcm (days)')
                        # plt.show()
                        
                        print(i/n)

                t2 = time.time()
                print(f"Time taken: {t2 - t1} seconds")

                # Save the results to a file
                try:
                    #with open('data/conj_params_annual_prop/conj_params_' + str(year) + '_no_artificial_tles.pkl', 'wb') as f:
                    with open(f"data/conj_params_monthly_prop/only_real/conj_params_{month:02d}{year}.pkl", "wb") as f:
                    #with open(f"data/conj_params_monthly_prop/starlink_only.pkl", "wb") as f:
                        pickle.dump((dist_a, dist_b, alt_a, alt_b, T_lcm), f)
                        print(f"File 'data/conj_params_{year}_no_artificial_tles.pkl' created!")
                except OSError as e:
                    if e.errno == 28:
                        print("Error: No space left on device. Please free up some space and try again.")
                    else:
                        print(f"An error occurred while saving the pickle file: {e}")


            # plot T_lcm as imshow
            # plt.figure()
            # plt.imshow(T_lcm/(24*3600), cmap='viridis', aspect='auto')
            # plt.colorbar()
            # plt.xlabel('Satellite Index')
            # plt.ylabel('Satellite Index')
            # plt.title('Synodic Period')
            # plt.show()

            # plot a histogram of the T_lcm values that are not nan
            # plt.figure()
            # plt.hist(T_lcm[~np.isnan(T_lcm)]/(24*3600), bins=1000)
            # # log scale
            # plt.yscale('log')
            # plt.xlabel('LCM Period (days)')
            # plt.ylabel('Frequency')
            # plt.title('Histogram of LCM Periods')
            # plt.show()

            # s_on_s = 0
            # total = 0
            # for i in range(len(sat_params)):
            #     for j in range(len(sat_params)):
            #         if T_lcm[i,j] < 0.2*(3600*24):
            #             # check to see if sat_params[i]['name'] and sat_params[j]['name'] include 'STARLINK'
            #             if 'STARLINK' in sat_params[i]['name'] and 'STARLINK' in sat_params[j]['name']:
            #                 s_on_s += 1
            #                 total += 1
            #             else: 
            #                 total += 1

            #     # print(i)

            # print(f"Frac STARLINK on STARLINK conjunctions: {s_on_s/total}")
                        
            # quit()

            # repeat the plot above, but make it a heatmap with number of points within each grid cell
            # make a 2D histogram of alt_a vs dist_a and alt_b vs dist_b
                # Filter out NaN values
            # valid_indices_a = ~np.isnan(alt_a.flatten()) & ~np.isnan(dist_a.flatten())
            # valid_indices_b = ~np.isnan(alt_b.flatten()) & ~np.isnan(dist_b.flatten())

            # alt_a_valid = alt_a.flatten()[valid_indices_a] - earth_radius
            # dist_a_valid = dist_a.flatten()[valid_indices_a]
            # alt_b_valid = alt_b.flatten()[valid_indices_b] - earth_radius
            # dist_b_valid = dist_b.flatten()[valid_indices_b]
            # T_lcm_valid = T_lcm.flatten()[valid_indices_a or valid_indices_b] 

            # make a combined dist and alt array
            # alt_valid = np.concatenate((alt_a_valid, alt_b_valid))
            # dist_valid = np.concatenate((dist_a_valid, dist_b_valid))
            
            ###############################################################################

            # Define the threshold distance for conjunction
            threshold = 0.1  # kilometers of miss distance
            t_tol = 1 # seconds of miss distance
            LEO_cutoff = 2000 + earth_radius # km

            print(f"Average distance between intersection points 1a and 2a: {np.nanmean(dist_a)}")
            print(f"Maximum distance between intersection points 1a and 2a: {np.nanmax(dist_a)}")
            print(f"Minimum distance between intersection points 1a and 2a: {np.nanmin(dist_a)}")
            print(f"Average distance between intersection points 1b and 2b: {np.nanmean(dist_b)}")
            print(f"Maximum distance between intersection points 1b and 2b: {np.nanmax(dist_b)}")
            print(f"Minimum distance between intersection points 1b and 2b: {np.nanmin(dist_b)}")
            count_dist_a = np.sum(dist_a < 0.1)
            count_dist_b = np.sum(dist_b < 0.1)
            print(f"Number of items in dist_a less than 0.1: {count_dist_a}")
            print(f"Number of items in dist_b less than 0.1: {count_dist_b}")

            # Find pairs of satellites with distances below the threshold
            print("Number of conjunctions before filtering: ", str(len(sat_params)**2))
            below_threshold_indices = np.where((dist_a < threshold) & (0.2*24*3600 < T_lcm) & (alt_a < LEO_cutoff))#(T_lcm < 90*24*3600)& 
            below_threshold_indices_b = np.where((dist_b < threshold) &  (0.2*24*3600 < T_lcm) & (alt_b < LEO_cutoff)) #(T_lcm < 90*24*3600) &
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
            #output_file = f"data/collision_results_monthly_prop/collision_results_{month:02d}{year}.csv"
            #output_file = f"data/collision_results_annual_prop/collision_results_{year}_no_artificial_tles.csv"
            output_file = f"data/collision_results_monthly_prop/only_real/collision_results_{month:02d}{year}.csv"
            #output_file = f"data/collision_results_monthly_prop/starlink_only.csv"

            # Determine whether conjunctions occur
            for i, j in tqdm(below_threshold_pairs, desc="Processing Pairs", unit="pair"):
                name1, norad_id1 = sat_params[i]['name'], sat_params[i]['norad']
                name2, norad_id2 = sat_params[j]['name'], sat_params[j]['norad']
                a_range1 = sat_params[i]['a_range']
                a_range2 = sat_params[j]['a_range']
                incl1 = sat_params[i]['incl'] * 180/np.pi
                incl2 = sat_params[j]['incl'] * 180/np.pi

                # Calculate expected time until collision
                # time_to_collision_a, intersection_point_a, T_lcm, collision = time_until_collision(tle1_line1, tle1_line2, tle2_line1, tle2_line2, a_or_b)
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

                # Calculate expected time until collision
                # time_to_collision_a, intersection_point_a, T_lcm, collision = time_until_collision(tle1_line1, tle1_line2, tle2_line1, tle2_line2, a_or_b)
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


                # writer.writerow([name1, norad_id1, a_range1, name2, norad_id2, a_range2, intersection_point_a, a_or_b, alt, collision, time_to_collision_days, T_lcm_days, 1/T_lcm_days])
                
                # for i, j in tqdm(below_threshold_pairs_b, desc="Processing Pairs", unit="pair"):
                #     tle1_line1, tle1_line2 = tle_l1[i], tle_l2[i]
                #     tle2_line1, tle2_line2 = tle_l1[j], tle_l2[j]
                #     name1 = obj_name[i]
                #     name2 = obj_name[j]
                #     norad_id1 = tle1_line2.split()[1]
                #     norad_id2 = tle2_line2.split()[1]
                #     a_or_b = 1

                #     satellite1 = get_satellite(tle1_line1, tle1_line2)
                #     satellite2 = get_satellite(tle2_line1, tle2_line2)
                #     a1 = satellite1.model.a * earth_radius
                #     a2 = satellite2.model.a * earth_radius
                #     ecc1 = satellite1.model.ecco
                #     ecc2 = satellite2.model.ecco

                #     rp1 = a1 * (1-ecc1) - earth_radius
                #     ra1 = a1 * (1+ecc1) - earth_radius
                #     rp2 = a2 * (1-ecc2) - earth_radius
                #     ra2 = a2 * (1+ecc2) - earth_radius
                #     a_range1 = [rp1, ra1]
                #     a_range2 = [rp2, ra2]

                #         # Calculate expected time until collision
                #     time_to_collision_b, intersection_point_b, T_lcm, collision = time_until_collision(tle1_line1, tle1_line2, tle2_line1, tle2_line2, a_or_b)
                    
                #     alt = np.linalg.norm(intersection_point_b) - earth_radius

                #     if time_to_collision_b is not None:
                #         time_to_collision_days = time_to_collision_b/86400
                #     else:
                #         time_to_collision_days = time_to_collision_b
                #     T_lcm_days = T_lcm/86400

                #     if time_to_collision_b is not None:
                #         collision_times.append(time_to_collision_b)
                #         intersection_points.append(intersection_point_b)
                #         altitudes.append(alt)

                #         if alt < rp1 or alt > ra1:
                #             out_of_range_1.append(alt - ra1)
                #         if alt < rp2 or alt > ra2:
                #             out_of_range_2.append(alt - ra2)
                #         collision_count = collision_count + 1
                #     else:
                #         none_count = none_count + 1

                #     saved_vars.append([name1, norad_id1, a_range1, name2, norad_id2, a_range2, intersection_point_b, a_or_b, alt, collision, time_to_collision_days, T_lcm_days, 1/T_lcm_days])


            # write saved_vars to a csv file
            with open(output_file, mode="w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(["Satellite 1","NORAD ID 1", "Altitude Range [km]", "Satellite 2", "NORAD ID 2", "Altitude Range [km]", "Conjunction Node Position [eci, km]", "Collision Point A or B", "Conjunction Point Altitude [km]", "Will they Collide?", "Time Until Collision [days]", "Synodic Period [days]", "Collision Frequency [1/days]", "Inclination of Sat 1", "Inclination of Sat 2"])
                writer.writerows(saved_vars)
                #exit()

            # with open(output_file, mode="w", newline="") as file:
            #     writer.writerow(["Satellite 1","NORAD ID 1", "Altitude Range [km]", "Satellite 2", "NORAD ID 2", "Altitude Range [km]", "Conjunction Node Position [eci, km]", "Collision Point A or B", "Conjunction Point Altitude [km]", "Will they Collide?", "Time Until Collision [days]", "Synodic Period [days]", "Collision Frequency [1/days]"])
            
            # writer.writerow([name1, norad_id1, a_range1, name2, norad_id2, a_range2, intersection_point_b, a_or_b, alt, collision, time_to_collision_days, T_lcm_days, 1/T_lcm_days])
            # repeat line above but do it for 
                # print(f"none count {none_count}")
                # print(f"collision count {collision_count}")
                # print(f"Out of Range 1 {len(out_of_range_1), max(out_of_range_1), min(out_of_range_1)}")
                # print(f"Out of Range 2 {len(out_of_range_2), max(out_of_range_2), min(out_of_range_2)}")

            """
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            
            # Extract x, y, z coordinates from the points
            x_coords = [point[0] for point in intersection_points]
            y_coords = [point[1] for point in intersection_points]
            z_coords = [point[2] for point in intersection_points]
            
            # Normalize collision times for colormap (optional, for better visual scaling)
            #norm_collision_times = (collision_times - np.min(collision_times)) / (np.max(collision_times) - np.min(collision_times))

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
            """


def time_until_collision_v2(sat_params1, sat_params2, T_lcm, a_or_b, t_tol):
        
        m1, ecc1, a1, period1 = sat_params1['meanan'], sat_params1['e'], sat_params1['a'], sat_params1['period']
        m2, ecc2, a2, period2 = sat_params2['meanan'], sat_params2['e'], sat_params2['a'], sat_params2['period']

        # satellite1 = get_satellite(tle1_line1, tle1_line2)
        # satellite2 = get_satellite(tle2_line1, tle2_line2)
        # M1 = np.deg2rad(satellite1.model.mo)
        # M2 = np.deg2rad(satellite2.model.mo)
        # ecc1 = satellite1.model.ecco
        # ecc2 = satellite2.model.ecco
        # a1 = satellite1.model.a * earth_radius
        # a2 = satellite2.model.a * earth_radius
        # mu = 398600.435507 #[km^3/s^2]
        # T1 = 2 * np.pi * np.sqrt(a1 ** 3/mu)
        # T2 = 2 * np.pi * np.sqrt(a2 ** 3/mu)

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
        #if a_or_b == 0:
        tau1 = (period1/(2 * np.pi)) * (E_col_1 - e1 - ecc1 * (np.sin(E_col_1) - np.sin(e1)))
        tau2 = (period2/(2 * np.pi)) * (E_col_2 - e2 - ecc2 * (np.sin(E_col_2) - np.sin(e2)))
        #if a_or_b == 1:
        #    tau1 = (T1/(2 * np.pi)) * (E_col_1 - E1 - ecc1 * (np.sin(E_col_1) - np.sin(E1)))
        #    tau2 = (T2/(2 * np.pi)) * (E_col_2 - E2 - ecc2 * (np.sin(E_col_2) - np.sin(E2)))

        #may not even need to calculate this here if just using 90 days
        # T_lcm = closest_rational_approx(period1, period2)

        #if T_lcm > 7776000:
        #    time_to_collision = None
        #    collision = False
        #    return time_to_collision, intersection_point, T_lcm, collision

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

        # intersection_within_threshold = []


        # tvec = np.arange(2000) 
        # t1a = tau1 + period1*tvec
        # t2a = tau2 + period2*tvec
        # if T_lcm < 7776000:
        #     t1a = t1a[t1a <= T_lcm + 10000]
        #     t2a = t2a[t2a <= T_lcm + 10000]
        # else:
        #     t1a = t1a[t1a <= 7786000]
        #     t2a = t2a[t2a <= 7786000]
        # #print(f"T_lcm:, {T_lcm}, Length of t1a: {len(t1a)}, Last values: {t1a[-4]}, {t1a[-3]}, {t1a[-2]}, {t1a[-1]}")
        # #print(f"T_lcm:, {T_lcm}, Length of t2a: {len(t2a)}, Last values: {t2a[-4]}, {t2a[-3]}, {t2a[-2]}, {t2a[-1]}")

        # intersection_within_threshold = []
        # i, j = 0, 0
        # time_to_collision = 0


        # while i < len(t1a) and j < len(t2a):
        #     if np.abs(t1a[i] - t2a[j]) <= 0.1: #0.1 seconds, collision threshold
        #         # If within threshold, add to the result
        #         break
        #     elif t1a[i] < t2a[j]:
        #         i += 1
        #     else:
        #         j += 1

        # if i < len(t1a) and t1a[i] < 7776000:
        #     time_to_collision = t1a[i]
        #     collision = True
        # else:
        #     time_to_collision = None
        #     collision = False

        #plt.scatter(vec, distance_between_sats)
        #plt.xlabel("Time from Epoch (seconds)")
        #plt.ylabel("Distance between satellites when Sat A is at conjunction point [km]")
        #plt.show()
        #print(f"time to collision {time_to_collision}")
        return time_to_collision, intersection_point, collision


    # Function to calculate time until collision
def time_until_collision(tle1_line1, tle1_line2, tle2_line1, tle2_line2, a_or_b):

        satellite1 = get_satellite(tle1_line1, tle1_line2)
        satellite2 = get_satellite(tle2_line1, tle2_line2)
        #M1 = np.deg2rad(satellite1.model.mo)
        #M2 = np.deg2rad(satellite2.model.mo)
        M1 = satellite1.model.mo
        M2 = satellite2.model.mo
        ecc1 = satellite1.model.ecco
        ecc2 = satellite2.model.ecco
        a1 = satellite1.model.a * earth_radius
        a2 = satellite2.model.a * earth_radius
        mu = 398600.435507 #[km^3/s^2]
        T1 = 2 * np.pi * np.sqrt(a1 ** 3/mu)
        T2 = 2 * np.pi * np.sqrt(a2 ** 3/mu)

        #calculate eccentric anomalies of starting points
        E1 = kepler_eqn(M1, ecc1)
        E2 = kepler_eqn(M2, ecc2)

        #use true anomalies to calculate eccentric anomalies of all possible collision points
        intersection_point_1a, intersection_point_1b, f1a, f1b, f2a, f2b = find_orbital_intersection(tle1_line1, tle1_line2, tle2_line1, tle2_line2)

        if a_or_b == 0:
            E_col_1 = 2 * np.arctan(np.sqrt((1 - ecc1) / (1 + ecc1)) * np.tan(f1a / 2))
            E_col_2 = 2 * np.arctan(np.sqrt((1 - ecc2) / (1 + ecc2)) * np.tan(f2a / 2))
            intersection_point = intersection_point_1a
        elif a_or_b == 1:
            E_col_1 = 2 * np.arctan(np.sqrt((1 - ecc1) / (1 + ecc1)) * np.tan(f1b / 2))
            E_col_2 = 2 * np.arctan(np.sqrt((1 - ecc2) / (1 + ecc2)) * np.tan(f2b / 2))
            intersection_point = intersection_point_1b

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

        #may not even need to calculate this here if just using 90 days
        T_lcm = closest_rational_approx(T1, T2)

        #if T_lcm > 7776000:
        #    time_to_collision = None
        #    collision = False
        #    return time_to_collision, intersection_point, T_lcm, collision

        tvec = np.arange(2000) 
        t1a = tau1 + T1*tvec
        t2a = tau2 + T2*tvec
        if T_lcm < 7776000:
            t1a = t1a[t1a <= T_lcm + 10000]
            t2a = t2a[t2a <= T_lcm + 10000]
        else:
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
            if np.abs(t1a[i] - t2a[j]) <= 1: #0.1 seconds, collision threshold
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

def closest_rational_approx(T1, T2, tol=1e-3):
    """Find the closest rational approximation of T1/T2 within a given tolerance."""
    ratio = T1 / T2  # Compute ratio first
    frac = Fraction.from_float(ratio).limit_denominator(int(1/tol))  # Use from_float for floating points
    return frac.denominator * T1  # Compute T_lcm

import numpy as np

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

# Function to compute the point on an orbit that lies on the intersection line
# def find_intersection_point(satellite, line_direction, step_size=1.0):
#     """
#     Find the point on the orbit of the satellite that intersects with the plane's intersection line.
#     """
#     ts = load.timescale()

#     # Search over one orbital period
#     t0 = satellite.epoch
#     orbital_period = satellite.model.no_kozai ** -1  # Approx orbital period in minutes

#     time_range = [t0.utc_datetime() + timedelta(minutes=step_size * i) for i in np.arange(0, orbital_period, step_size)]

#     for time in time_range:
#         t = ts.utc(time)
#         geocentric = satellite.at(t).position.km  # Get position of the satellite in km
        
#         # Check the distance from the satellite position to the intersection line
#         distance_to_line = np.linalg.norm(np.cross(geocentric, line_direction))
#         if distance_to_line < 1e-6:  # Close to the line
#             return geocentric

#     return None

# def find_intersection_point_analytically(satellite, line_direction):
#     # Extract Keplerian elements
#     a = satellite.model.a  # Semi-major axis in Earth radii
#     e = satellite.model.ecco  # Eccentricity
#     i = satellite.model.inclo  # Inclination in radians
#     raan = satellite.model.nodeo  # Right ascension of ascending node in radians
#     arg_pe = satellite.model.argpo  # Argument of perigee in radians

#     # Compute the apse line in the perifocal frame
#     apse_line_perifocal = np.array([np.cos(arg_pe), np.sin(arg_pe), 0])

#     # Rotation matrices to transform from perifocal to ECI frame
#     R3_W = np.array([[np.cos(raan), -np.sin(raan), 0],
#                      [np.sin(raan), np.cos(raan), 0],
#                      [0, 0, 1]])
#     R1_i = np.array([[1, 0, 0],
#                      [0, np.cos(i), -np.sin(i)],
#                      [0, np.sin(i), np.cos(i)]])
#     R3_w = np.array([[np.cos(arg_pe), -np.sin(arg_pe), 0],
#                      [np.sin(arg_pe), np.cos(arg_pe), 0],
#                      [0, 0, 1]])

#     # Combined rotation matrix
#     R = R3_W @ R1_i @ R3_w

#     # Transform the apse line to the ECI frame
#     apse_line_eci = R @ apse_line_perifocal

#     # Normalize line_direction and apse_line_eci
#     line_direction = line_direction / np.linalg.norm(line_direction)
#     apse_line_eci = apse_line_eci / np.linalg.norm(apse_line_eci)

#     # Compute the true anomaly for the intersection points
#     angle = np.arccos(np.dot(line_direction, apse_line_eci))
#     true_anomaly_a = angle
#     true_anomaly_b = -angle

#     # Compute the distance of the orbit at the true anomalies
#     d_a = a * (1 - e ** 2) / (1 + e * np.cos(true_anomaly_a))
#     d_b = a * (1 - e ** 2) / (1 + e * np.cos(true_anomaly_b))

#     # Compute the position vectors in the perifocal frame
#     pos_perifocal_a = d_a * np.array([np.cos(true_anomaly_a), np.sin(true_anomaly_a), 0])
#     pos_perifocal_b = d_b * np.array([np.cos(true_anomaly_b), np.sin(true_anomaly_b), 0])

#     # Transform the position vectors to the ECI frame
#     pos_eci_a = R @ pos_perifocal_a
#     pos_eci_b = R @ pos_perifocal_b

#     return pos_eci_a, pos_eci_b

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
    d = (satellite.model.a)*earth_radius * (1 - satellite.model.ecco ** 2) / (1 + satellite.model.ecco * np.cos(angle))
    # compute the position vector of the satellite at d distance from the focus of the orbit along the line_direction
    pos_point_a = d * line_direction

    # repeat the same process for the other intersection point (switch the line direction)
    line_direction = -line_direction
    angle = np.arccos(np.dot(line_direction, apse_line))
    d = (satellite.model.a)*earth_radius * (1 - satellite.model.ecco ** 2) / (1 + satellite.model.ecco * np.cos(angle))
    pos_point_b = d * line_direction

    return pos_point_a, pos_point_b, true_anomaly_a, true_anomaly_b

def find_intersection_point_analytically_v2(sat_params_i, line_direction):
    # find the angle between the line of intersection and the satellite's apse line (line in the direction of eccentricity)
    # apse_line_pf = np.array([np.cos(satellite.model.argpo), np.sin(satellite.model.argpo), 0])
    # convert apse line from perifocal to ECI frame
    apse_line_pf = np.array([1,0,0])
    apse_line = pf_to_eci_v2(apse_line_pf, sat_params_i)

    angle = np.arccos(np.dot(line_direction, apse_line))
    true_anomaly_a = angle
    true_anomaly_b = -angle
    # compute the distance of the orbit when the true anomaly is the angle between the line of intersection and the apse line
    # d = (satellite.model.a)*earth_radius * (1 - satellite.model.ecco ** 2) / (1 + satellite.model.ecco * np.cos(angle))
    d = sat_params_i['a']*(1-sat_params_i['e']**2)/(1+sat_params_i['e']*np.cos(angle))
    # compute the position vector of the satellite at d distance from the focus of the orbit along the line_direction
    pos_point_a = d * line_direction

    # repeat the same process for the other intersection point (switch the line direction)
    line_direction = -line_direction
    angle = np.arccos(np.dot(line_direction, apse_line))
    # d = (satellite.model.a)*earth_radius * (1 - satellite.model.ecco ** 2) / (1 + satellite.model.ecco * np.cos(angle))
    d = sat_params_i['a']*(1-sat_params_i['e']**2)/(1+sat_params_i['e']*np.cos(angle))
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

def pf_to_eci_v2(pos_pf, sat_params_i):
    # convert from perifocal frame to ECI frame position vector
    # Extract Keplerian elements
    # i = satellite.model.inclo  # Inclination in radians
    # raan = satellite.model.nodeo  # Right ascension of ascending node in radians
    # arg_pe = satellite.model.argpo  
    i, raan, arg_pe = sat_params_i['incl'], sat_params_i['raan'], sat_params_i['arg_per'] # Arguments in radians
    pos_eci = np.zeros(3)
    pos_eci[0] = (np.cos(raan)*np.cos(arg_pe) - np.sin(raan)*np.sin(arg_pe)*np.cos(i))*pos_pf[0] + (-np.cos(raan)*np.sin(arg_pe) - np.sin(raan)*np.cos(arg_pe)*np.cos(i))*pos_pf[1]
    pos_eci[1] = (np.sin(raan)*np.cos(arg_pe) + np.cos(raan)*np.sin(arg_pe)*np.cos(i))*pos_pf[0] + (-np.sin(raan)*np.sin(arg_pe) + np.cos(raan)*np.cos(arg_pe)*np.cos(i))*pos_pf[1]
    pos_eci[2] = (np.sin(i)*np.sin(arg_pe))*pos_pf[0] + (np.sin(i)*np.cos(arg_pe))*pos_pf[1]
    return pos_eci

# def find_intersection_point_analytically(satellite, line_direction):
#     """
#     Find the point on the orbit of the satellite that intersects with the plane's intersection line 
#     using orbital elements instead of propagation.
#     """
#     # Orbital elements
#     a = satellite.model.a * 6378.137  # Semi-major axis (Earth radii to km conversion)
#     e = satellite.model.ecco  # Eccentricity
#     i = np.radians(satellite.model.inclo)  # Inclination (degrees to radians)
#     omega = np.radians(satellite.model.argpo)  # Argument of perigee (degrees to radians)
#     raan = np.radians(satellite.model.nodeo)  # RAAN (degrees to radians)

#     # Mean anomaly (assume a circular orbit for simplicity)
#     M = np.linspace(0, 2 * np.pi, num=100)  # Mean anomaly array

#     intersection_point = None
#     min_distance = float('inf')

#     for mean_anomaly in M:
#         # Eccentric anomaly using Kepler's equation
#         E = mean_anomaly + e * np.sin(mean_anomaly)  # First-order approximation

#         # True anomaly
#         true_anomaly = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))

#         # Radius at true anomaly
#         r = a * (1 - e * np.cos(E))

#         # Position in orbital plane (x', y', z')
#         x_prime = r * np.cos(true_anomaly)
#         y_prime = r * np.sin(true_anomaly)
#         z_prime = 0  # Assuming the orbital plane is XY

#         # Rotate to the inertial frame
#         # Applying rotations for RAAN and inclination
#         x = (np.cos(raan) * np.cos(omega) - np.sin(raan) * np.sin(omega) * np.cos(i)) * x_prime + \
#             (-np.cos(raan) * np.sin(omega) - np.sin(raan) * np.cos(omega) * np.cos(i)) * y_prime
#         y = (np.sin(raan) * np.cos(omega) + np.cos(raan) * np.sin(omega) * np.cos(i)) * x_prime + \
#             (-np.sin(raan) * np.sin(omega) + np.cos(raan) * np.cos(omega) * np.cos(i)) * y_prime
#         z = (np.sin(i) * np.sin(omega)) * x_prime + (np.sin(i) * np.cos(omega)) * y_prime

#         position = np.array([x, y, z])

#         # Calculate distance from the intersection line
#         distance_to_line = np.linalg.norm(np.cross(position, line_direction))

#         # Check if this position is the closest to the line direction
#         if distance_to_line < min_distance:
#             min_distance = distance_to_line
#             intersection_point = position

#     return intersection_point

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

    # Create time range for the orbits
    # t0 = satellite1.epoch
    # orbital_period1 = 2*np.pi*np.sqrt((satellite1.model.a*6378.137)**3/398600.4418)  # Orbital period in seconds
    # orbital_period2 = 2*np.pi*np.sqrt((satellite2.model.a*6378.137)**3/398600.4418)
    # orbital_period1 = satellite1.model.no_kozai ** -1 * 60  # Orbital period in seconds
    # orbital_period2 = satellite2.model.no_kozai ** -1 * 60

    # time_range1 = [t0.utc_datetime() + timedelta(seconds=orbital_period1/1000 * i) for i in np.linspace(0, orbital_period1/2, 1000)]
    # time_range2 = [t0.utc_datetime() + timedelta(seconds=orbital_period2/1000 * i) for i in np.linspace(0, orbital_period2/2, 1000)]

    # Get orbit positions for each satellite
    positions1 = compute_positions_keplerian(satellite1)
    positions2 = compute_positions_keplerian(satellite2)
    # positions1 = np.array([satellite1.at(ts.utc(time)).position.km for time in time_range1])
    # positions2 = np.array([satellite2.at(ts.utc(time)).position.km for time in time_range2])

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

# Main function to find intersection line and points
def find_orbital_intersection(tle1_line1, tle1_line2, tle2_line1, tle2_line2):
    t1 = time.time()
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
    #print(f"Intersection line direction: {intersection_line}")

    # Find the intersection points on both orbits
    # intersection_point_1 = find_intersection_point(satellite1, intersection_line)
    # intersection_point_2 = find_intersection_point(satellite2, intersection_line)

    # Find the intersection points on both orbits using the new function
    intersection_point_1a, intersection_point_1b, true_anomaly_1a, true_anomaly_1b = find_intersection_point_analytically(satellite1, intersection_line)
    intersection_point_2a, intersection_point_2b, true_anomaly_2a, true_anomaly_2b = find_intersection_point_analytically(satellite2, intersection_line)

    t2 = time.time()
    #print(f"Time taken: {t2 - t1} seconds")
    '''
    print(f"Intersection point on orbit 1a: {intersection_point_1a}")
    print(f"Intersection point on orbit 1b: {intersection_point_1b}")
    print(f"Intersection point on orbit 2a: {intersection_point_2a}")
    print(f"Intersection point on orbit 2b: {intersection_point_2b}")

    # print the distance between intersection points 1a and 2a, 1b and 2b
    print('Closest approach distance between points a: ' + str(np.linalg.norm(intersection_point_1a - intersection_point_2a)))
    print('Closest approach distance between points b: ' + str(np.linalg.norm(intersection_point_1b - intersection_point_2b)))
    '''
    # Plot the orbits and intersection
    #plot_orbits_and_intersection(satellite1, satellite2, intersection_line, intersection_point_1a, intersection_point_1b, intersection_point_2a, intersection_point_2b)
    return intersection_point_1a, intersection_point_1b, true_anomaly_1a, true_anomaly_1b, true_anomaly_2a, true_anomaly_2b


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

    # inclination1, raan1, period1 = sat_params1['incl'], sat_params1['raan'], sat_params1['period']
    # inclination2, raan2, period2 = sat_params2['incl'], sat_params2['raan'], sat_params2['period']
    
    # # Compute normal vectors for each orbital plane
    # normal1 = orbital_plane_normal(inclination1, raan1)
    # normal2 = orbital_plane_normal(inclination2, raan2)

    # # Compute the line of intersection between the two planes
    # intersection_line = compute_plane_intersection(normal1, normal2)

    # # Find the intersection points on both orbits using the new function
    # intersection_point_1a, intersection_point_1b, true_anomaly_1a, true_anomaly_1b = find_intersection_point_analytically_v2(sat_params1, intersection_line)
    # intersection_point_2a, intersection_point_2b, true_anomaly_2a, true_anomaly_2b = find_intersection_point_analytically_v2(sat_params2, intersection_line)

    dist_a = np.linalg.norm(intersection_point_1a - intersection_point_2a)
    dist_b = np.linalg.norm(intersection_point_1b - intersection_point_2b)
    alt_a = np.mean([np.linalg.norm(intersection_point_1a), np.linalg.norm(intersection_point_2a)])
    alt_b = np.mean([np.linalg.norm(intersection_point_1b), np.linalg.norm(intersection_point_2b)])

    T_lcm = closest_rational_approx(period1, period2)

    return dist_a, dist_b, alt_a, alt_b, intersection_point_1a, intersection_point_1b, T_lcm


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

    min_dist = min(np.linalg.norm(intersection_point_1a - intersection_point_2a), np.linalg.norm(intersection_point_1b - intersection_point_2b))
    arg = np.argmin([np.linalg.norm(intersection_point_1a - intersection_point_2a), np.linalg.norm(intersection_point_1b - intersection_point_2b)])
    conj_alt = np.linalg.norm(intersection_point_1a) if arg == 0 else np.linalg.norm(intersection_point_1b)

    dist_a = np.linalg.norm(intersection_point_1a - intersection_point_2a)
    dist_b = np.linalg.norm(intersection_point_1b - intersection_point_2b)
    alt_a = np.mean([np.linalg.norm(intersection_point_1a), np.linalg.norm(intersection_point_2a)])
    alt_b = np.mean([np.linalg.norm(intersection_point_1b), np.linalg.norm(intersection_point_2b)])
    # print the distance between intersection points 1a and 2a, 1b and 2b

    a1 = satellite1.model.a * earth_radius
    a2 = satellite2.model.a * earth_radius
    mu = 398600.435507 #[km^3/s^2]
    T1 = 2 * np.pi * np.sqrt(a1 ** 3/mu)
    T2 = 2 * np.pi * np.sqrt(a2 ** 3/mu)
    T_lcm = closest_rational_approx(T1, T2)

    return dist_a, dist_b, alt_a, alt_b, T_lcm

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
    
    # print(i)
    # if np.any(np.isnan(dist_a)):
    #     print(f"Error at index {i}")

    return dist_a, dist_b, alt_a, alt_b, T_lcm


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
        if a_range[j][0] > (a_range[i][1] + 100) or a_range[j][1] < (a_range[i][0] - 100):
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
    a = satellite.model.a*earth_radius
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