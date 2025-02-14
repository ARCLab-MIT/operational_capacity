import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime, timedelta, timezone
from sgp4.api import Satrec, WGS72
from skyfield.api import load
from astropy import units as u
from poliastro.bodies import Earth

def plot_3d_orbit(tle_line1, tle_line2, num_points=1000):
    # Load satellite from TLE
    satellite = Satrec.twoline2rv(tle_line1, tle_line2, WGS72)
    
    # Define a timespan for propagation
    ts = load.timescale()
    now = datetime.now(timezone.utc)  # Ensure datetime is explicitly in UTC
    times = [ts.utc((now + timedelta(minutes=i)).year,
                    (now + timedelta(minutes=i)).month,
                    (now + timedelta(minutes=i)).day,
                    (now + timedelta(minutes=i)).hour,
                    (now + timedelta(minutes=i)).minute,
                    (now + timedelta(minutes=i)).second)
             for i in np.linspace(0, 60*3, num_points)]  # 1 day coverage

    # Extract positions
    x_vals, y_vals, z_vals = [], [], []
    for time in times:
        error_code, pos, _ = satellite.sgp4(time.tt, 0.0)
        if error_code == 0:
            x_vals.append(pos[0] * 1000)  # Convert km to meters
            y_vals.append(pos[1] * 1000)
            z_vals.append(pos[2] * 1000)

    # plt.show()
    return(x_vals, y_vals, z_vals)

# get list of object names to pull TLEs for

name_list = ['COSMOS 2058', 'STARLINK-3200', 'STARLINK-3632', 'STARLINK-3557', 'STARLINK-3539', 'STARLINK-3568', 'STARLINK-3567', 'STARLINK-3647', 'STARLINK-3538', 'STARLINK-3964', 'STARLINK-3948', 'STARLINK-3953', 'STARLINK-3924', 'STARLINK-3895', 'STARLINK-3922', 'STARLINK-3906', 'STARLINK-3867', 'STARLINK-3936', 'STARLINK-4078', 'STARLINK-4207', 'STARLINK-4298', 'STARLINK-4534', 'STARLINK-4515', 'STARLINK-4512', 'STARLINK-4525', 'STARLINK-4491', 'STARLINK-5290', 'STARLINK-5122', 'STARLINK-5250']

# open the TLE file and get the indices of the lines that contain these names
tle_file = open('3le_leo_092424.txt', 'r')
tle_lines = tle_file.readlines()

line1 = []
line2 = []
ln_idx = np.zeros(len(name_list))
for i in range(len(name_list)):
    # find the line of the TLE that contains the name of the object
    idx = np.where([name_list[i] in line for line in tle_lines])[0][0]
    ln_idx[i] = idx
    line1.append(tle_lines[idx+1][:-1])
    line2.append(tle_lines[idx+2][:-1])

# for i in range(int(len(tle_lines)/3)):
#     obj_name[i] = tle_lines[3*i][2:-1]
#     obj_norad[i] = int(tle_lines[3*i+1][2:7])
#     tle_l1[i] = tle_lines[3*i+1][:-1]
#     tle_l2[i] = tle_lines[3*i+2][:-1]

# # Example usage with ISS TLE
# tle_line1 = "1 25544U 98067A   24042.54828704  .00016717  00000+0  30542-3 0  9992"
# tle_line2 = "2 25544  51.6424  52.6872 0002719  42.8774  21.0742 15.50233897444043"

# line1 = []
# line2 = []

# # 0 TIROS 5
# line1.append("1 00309U 62025A   24268.51244858  .00001918  00000-0  34782-3 0  9992")
# line2.append("2 00309  58.0916  15.9403 0192770 116.2893 245.8082 14.59068971289728")

# # 0 HORYU 2
# line1.append("1 38340U 12025D   24268.57500579  .00008102  00000-0  99411-3 0  9996")
# line2.append("2 38340  98.1670  10.4269 0010346  80.1604 280.0777 14.83334739665368")

# # 0 STARLINK-5276
# line1.append("1 55399U 23014J   24268.52097751  .00001474  00000-0  13129-3 0  9994")
# line2.append("2 55399  69.9997 323.6351 0003070 269.9764  90.1043 14.98338058 92498")


# Create 3D plot
r_e = 6378.137 # Earth's radius in km
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")
for i in range(len(line1)):
    x, y, z = plot_3d_orbit(line1[i], line2[i])
    if i == 0:
        ax.plot(x, y, z, "k", label="Satellite Orbit", linewidth=5)
    else:
        ax.plot(x, y, z, "lightgray", label="Satellite Orbit", linewidth=1)
    print(i / len(line1) * 100, '% done')

# Plot Earth
r_e = 6378.137*1e3  # Earth's radius in km
u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
xx = r_e * np.cos(u) * np.sin(v)
yy = r_e * np.sin(u) * np.sin(v)
zz = r_e * np.cos(v)
ax.plot_surface(xx, yy, zz, color='tab:blue', alpha=0.2)

# Labels and title
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")
ax.set_title("3D Satellite Orbit from TLE")

# Remove background grid and color
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.grid(False)

# Remove background color
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')

# Set equal aspect ratio
ax.set_box_aspect([1, 1, 1])

# Show plot
plt.show()

import plotly.graph_objects as go


# Create 3D plot with Plotly
fig = go.Figure()

# Plot Earth
u, v = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
r_earth = 6378.137 * 1000  # Earth's radius in meters
earth_x = r_earth * np.cos(u) * np.sin(v)
earth_y = r_earth * np.sin(u) * np.sin(v)
earth_z = r_earth * np.cos(v)

# fig.add_trace(go.Surface(x=earth_x, y=earth_y, z=earth_z, colorscale='Greys', opacity=0.75, showscale=False))
fig.add_trace(go.Surface(x=earth_x, y=earth_y, z=earth_z, surfacecolor=np.ones_like(earth_x), colorscale=[[0, 'cornflowerblue'], [1, 'cornflowerblue']], opacity=0.6, showscale=False))
# Plot satellite orbits
# x_list = np.zeros((len(line1), 1000))
# y_list = np.zeros((len(line1), 1000))
# z_list = np.zeros((len(line1), 1000))
for i in range(len(line1)):
    x, y, z = plot_3d_orbit(line1[i], line2[i])
    # x_list[i,:] = x
    # y_list[i,:] = y
    # z_list[i,:] = z
    # if i == 0, make line black and thick
    # fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', name=f'Satellite {i+1}', line=dict(width=8 if i == 0 else 4)))
    # repeat above line, but make the colors i = 1,black, i = 2, red, i = 3, blue
    fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', name=f'Satellite {i+1}', line=dict(width=10 if i == 0 else 1, color='black' if i == 0 else 'darkblue')))

# conj_point = []
# conj_count = 0
# # check for any instances when one satellite is within 100 m of another satellite (norm of difference of position vectors at any time is less than 100 m)
# for i in range(len(line1)):
#     for j in range(i+1, len(line1)):
#         for k in range(1000):
#             if np.linalg.norm(np.array([x_list[i,k], y_list[i,k], z_list[i,k]]) - np.array([x_list[j,k], y_list[j,k], z_list[j,k]])) < 1000:
#                 # record the position of the point
#                 conj_point.append([x_list[i,k], y_list[i,k], z_list[i,k]])
#                 conj_count = conj_count + 1
#                 print(conj_count)

# plot a red dot at the point of conjunction
# conj_point = np.array(conj_point)
# fig.add_trace(go.Scatter3d(x=conj_point[:,0], y=conj_point[:,1], z=conj_point[:,2], mode='markers', name='Conjunction Point', marker=dict(size=5, color='red')))


# Labels and title
fig.update_layout(
    scene=dict(
        # xaxis_title='X (m)',
        # yaxis_title='Y (m)',
        # zaxis_title='Z (m)',
        xaxis=dict(showbackground=False),
        yaxis=dict(showbackground=False),
        zaxis=dict(showbackground=False),

    ),
    # make axis labels white

    title='3D Satellite Orbit from TLE',
    showlegend=True
)
# save figure as pdf
# fig.write_image("3D_Satellite_Orbit.pdf")
fig.write_image("3d_sat_orbit_real.pdf", engine="kaleido", scale = 5)

fig.show()