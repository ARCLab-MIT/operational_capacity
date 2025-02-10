import time
from math import gcd, ceil
from sympy import mod_inverse  # SymPy's modular inverse function
import numpy as np
import matplotlib.pyplot as plt 
import math
# make font arial and text, plot edges, and labels dimgray
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['text.color'] = 'dimgray'
plt.rcParams['axes.labelcolor'] = 'dimgray'
plt.rcParams['xtick.color'] = 'dimgray'
plt.rcParams['ytick.color'] = 'dimgray'
plt.rcParams['axes.edgecolor'] = 'dimgray'

def main():
    T1 = 5400  # Period of Satellite 1 (seconds)
    print(f"Period of Satellite 1: {T1/60} mins.")
    T2 = 5500  # Period of Satellite 2 (seconds)
    print(f"Period of Satellite 2: {T2/60} mins.")
    f = np.linspace(0,1,10000)#0.9  # Initial offset as a fraction of T2
    print(f"Number of initial states considered: {len(f)}")
    threshold = 100  # Time difference threshold for collision (seconds)
    print(f"Collision time difference threshold: {threshold} seconds.")

    t1 = time.time()
    collision_time = np.zeros(len(f))
    for i in range(len(f)):
        collision_time[i], n, m, lcm_yrs = find_collision_time_v3(T1, T2, f[i], threshold)
        # print(f"The first collision occurs at t = {collision_time[i]:.2f} seconds.")
        # print(f"n = {n:.0f}")
        # print(f"m = {m:.0f}")
    t2 = time.time()
    print(f"LCM days is {lcm_yrs*365.25}")
    
    print(f"Elapsed time: {t2 - t1:.2f} seconds.")

    # Assuming collision_time is a numpy array or list
    collision_time_years = collision_time / (3600 * 24 * 365.25)
    collision_time_years_nonzeros = collision_time_years[collision_time_years > 0]
    perc_nonzeros = len(collision_time_years_nonzeros) / len(collision_time_years) * 100
    print(f"Percentage of non-infinite collision times: {perc_nonzeros:.2f}%")

    plt.figure(figsize = (4,3))
    plt.hist(collision_time_years_nonzeros, bins=100, color='black')

    # Convert counts to percentages
    # percentages = counts * 100

    # Update the heights of the patches to reflect percentages
    # for patch, percentage in zip(patches, percentages):
    #     patch.set_height(percentage)

    plt.xlabel('Time to Collision [years]')
    plt.grid(True,axis = 'y', which = 'major', linestyle='-', lw=0.5, color = 'dimgray')
    plt.ylabel('Frequency')
    # plt.yscale('log')
    plt.tight_layout()
    plt.show()

    # compute the weighted average of the collision time
    weighted_avg = np.sum(collision_time_years_nonzeros) / len(collision_time_years_nonzeros)

    print(f"Weighted average of collision time: {weighted_avg:.2f} years")


def find_collision_time_v3(T1, T2, f, threshold):
    """
    Finds the first collision time for two satellites.

    Parameters:
    - T1: Orbital period of Satellite 1 (seconds).
    - T2: Orbital period of Satellite 2 (seconds).
    - f: Fraction of T2 away from the conjunction point at t=0 (0 <= f < 1).
    - threshold: Time difference threshold for collision (seconds).

    Returns:
    - t_collision: Time of the first collision (seconds).
    """
    # Calculate the least common multiple (LCM) of T1 and T2
    lcm_T1_T2 = math.lcm(T1, T2)
    lcm_years = lcm_T1_T2 / (3600 * 24 * 365.25)
    # print(f"LCM of T1 and T2: {lcm_T1_T2} seconds = {lcm_years} years")

    # Initialize arrival times for each satellite
    t_sat1 = 0  # Satellite 1 starts at conjunction point
    t_sat2 = f * T2  # Satellite 2 starts offset by f * T2

    # Indices for each satellite's arrival
    n, m = 0, 0

    # Loop to find the first collision
    while t_sat1 <= lcm_T1_T2:# and t_sat1 <= t_cutoff:
        # Check the absolute time difference
        if abs(t_sat1 - t_sat2) < threshold:
            return min(t_sat1, t_sat2), n, m, lcm_years

        # Advance the satellite with the earlier arrival time
        if t_sat1 < t_sat2:
            n += 1
            t_sat1 = n * T1
        else:
            m += 1
            t_sat2 = m * T2 + f * T2

    return 0, 0, 0, 0

if __name__ == "__main__":
    main()
    