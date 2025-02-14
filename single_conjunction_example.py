import numpy as np
from skyfield.api import EarthSatellite, load
import time
from fractions import Fraction


def main():
    # Example TLEs
    # 0 STARLINK-30164
    # tle1_line1 = "1 57500U 23113F   24268.60933512  .00002345  00000-0  18858-3 0  9995"
    # tle1_line2 = "2 57500  43.0055  22.9320 0001551 273.9909  86.0756 15.02545041 63777"

    # # 0 EAGLESAT 1
    # tle2_line1 = "1 43018U 17073F   24268.62229613  .00041108  00000-0  16974-2 0  9990"
    # tle2_line2 = "2 43018  97.5772 152.2690 0153034 285.0896  73.3447 15.13765642371092"

    # 0 EXPLORER 7
    tle1_line1 = '00022U 59009A   24268.53169352  .00011671  00000-0  62956-3 0  9990'
    tle1_line2 = '2 00022  50.2766 296.9541 0106313 136.0730 224.8736 15.10389732704408'

    # 0 CHAMRAN-1
    tle2_line1 = '1 61072U 24165B   24268.57911644  .00010710  00000-0  68133-3 0  9996'
    tle2_line2 = '2 61072  64.0910 306.1383 0020336 301.9004  58.0132 15.09043265  1558'



    # tle1_line1 = "1 25544U 98067A   20264.89222821  .00001264  00000-0  29661-4 0  9991"
    # tle1_line2 = "2 25544  51.6456  21.0985 0005571  41.2232 318.8816 15.49114663246565"

    # tle2_line1 = "1 43205U 18015A   20264.87073778  .00000647  00000-0  18304-4 0  9994"
    # tle2_line2 = "2 43205  53.0551  55.6069 0012501  66.8434 293.3607 15.08856543157085"
    intersection_point, f1a, f1b, f2a, f2b = find_orbital_intersection(tle1_line1, tle1_line2, tle2_line1, tle2_line2)

    alt = np.linalg.norm(intersection_point) - 6378.15

    print(alt)

# Helper function to parse TLEs and return satellite object
def get_satellite(tle_line1, tle_line2):
    satellite = EarthSatellite(tle_line1, tle_line2, 'Satellite', load.timescale())
    return satellite

# Main function to find intersection line and points
def find_orbital_intersection(tle1_line1, tle1_line2, tle2_line1, tle2_line2):
    t1 = time.time()
    # Load the satellites from TLEs
    satellite1 = get_satellite(tle1_line1, tle1_line2)
    satellite2 = get_satellite(tle2_line1, tle2_line2)

    # Get orbital parameters (inclination and RAAN) from TLEs
    inclination1, raan1 = satellite1.model.inclo, satellite1.model.nodeo
    inclination2, raan2 = satellite2.model.inclo, satellite2.model.nodeo
    a1 = satellite1.model.a*6378.15
    a2 = satellite2.model.a*6378.15

    mu = 398600.15
    # Compute the LCM period
    T1 = 2 * np.pi * np.sqrt(a1 ** 3 / mu)
    T2 = 2 * np.pi * np.sqrt(a2 ** 3 / mu)

    # Compute the closest rational approximation of T1/T2
    T_lcm = closest_rational_approx(T1, T2)


    print('T_lcm (d): ', T_lcm/(3600*24))

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
    return intersection_point_1a, true_anomaly_1a, true_anomaly_1b, true_anomaly_2a, true_anomaly_2b

def closest_rational_approx(T1, T2, tol=1e-3):
        """Find the closest rational approximation of T1/T2 within a given tolerance."""
        ratio = T1 / T2  # Compute ratio first
        frac = Fraction.from_float(ratio).limit_denominator(int(1/tol))  # Use from_float for floating points
        t1 = frac.numerator * T2  # Compute T_syn
        t2 = frac.denominator * T1  # Compute T_syn
        print("T_lcm error (s): ", abs(t1-t2))
        return frac.denominator * T1  # Compute T_lcm

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

main()