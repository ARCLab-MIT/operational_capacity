# read csv collision_results
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
import ast
from skyfield.api import Distance, load, wgs84
from skyfield.positionlib import Geocentric
import datetime as dt
from tabulate import tabulate
# make figures arial font
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
import matplotlib.colors as mcolors
import matplotlib as mpl
import colorsys
import matplotlib.colors as colors 
# matplotlib.use('Agg') # Set a non-interactive backend like 'Agg' 

# don't display any plots 
plt.ioff()

def read_csv(file_name):
    with open(file_name, 'r') as f:
        reader = csv.reader(f)
        data = list(reader)
    return data

def eci2lla(x,y,z, t_dt):
    ts = load.timescale()
    year, month, day, hour, minute = t_dt.year, t_dt.month, t_dt.day, t_dt.hour, t_dt.minute
    t = ts.utc(year, month, day, hour, minute)
    d = Distance(m=[x, y, z])
    p = Geocentric(d.au, t=t)
    g = wgs84.subpoint(p)
    return g.latitude.degrees, g.longitude.degrees, g.elevation.m

# data = read_csv("data/collision_results_wp_0p1_1.csv")  
#data = read_csv('data/collision_results_monthly_prop/only_real/collision_results_012020.csv')
data = read_csv('col_result_6month/collision_results_012025.csv')

# for each line in the data file, check to see if the synodic period is less than 90 days. If so, keep it, otherwise discard it.
# data_shortened = []
# for i in range(1, len(data)):
#     # data[i][9] = float(data[i][9])
#     if float(data[i][9]) < 90:
#         data_shortened.append(data[i])

# # save data_shortened to a new csv file
# with open('collision_results_90min_100m_shortened.csv', 'w', newline='') as f:
#     writer = csv.writer(f)
#     writer.writerows(data_shortened)

r_e = 6378.137 # radius of the Earth in km
# data is a list of lists. Make it into a 2d np array
data_np = np.array(data)
name1 = data_np[1:, 0]
norad1 = (data_np[1:, 1]).astype(int)
alt_range1_str = data_np[1:, 2]
alt_range1 = np.zeros((len(alt_range1_str), 2))
incl1 = (data_np[1:, 13]).astype(float)
for i in range(len(alt_range1_str)):
    items = (alt_range1_str[i][1:-1]).split()
    alt_range1[i,:] = np.array([float(items[0][:-1])-r_e, float(items[1])-r_e])
# assign alt_range to a dictionary where the key is the corresponding norad ID. If the norad ID is already in the dictionary, skip it. 
name2 = data_np[1:, 3]
norad2 = (data_np[1:, 4]).astype(int)
alt_range2_str = data_np[1:, 5]
alt_range2 = np.zeros((len(alt_range2_str), 2))
incl2 = (data_np[1:, 14]).astype(float)
for i in range(len(alt_range2_str)):
    items = (alt_range2_str[i][1:-1]).split()
    alt_range2[i,:] = np.array([float(items[0][:-1])-r_e, float(items[1])-r_e])

alt_range_dict = {}
incl_dict = {}
for i in range(len(norad1)):
    if norad1[i] not in alt_range_dict:
        alt_range_dict[norad1[i]] = alt_range1[i]
    else:
        continue
for i in range(len(norad2)):
    if norad2[i] not in alt_range_dict:
        alt_range_dict[norad2[i]] = alt_range2[i]
    else:
        continue

for i in range(len(norad1)):
    if norad1[i] not in incl_dict:
        incl_dict[norad1[i]] = incl1[i]
    else:
        continue
for i in range(len(norad2)):
    if norad2[i] not in incl_dict:
        incl_dict[norad2[i]] = incl2[i]
    else:
        continue


cnode_pos_str = data_np[1:, 6]
cnode_pos = np.zeros((len(cnode_pos_str), 3))
for i in range(len(cnode_pos_str)):
    items = (cnode_pos_str[i][1:-1]).split()
    cnode_pos[i,:] = np.array([float(items[0]), float(items[1]), float(items[2])])
# Convert the list of lists to a 2D numpy array
cnode_alt = np.where(data_np[1:, 8] == '', '0.0', data_np[1:, 8]).astype(float)
collision_yn = data_np[1:, 9]
ttc = np.where(data_np[1:, 10] == '', '0', data_np[1:, 10]).astype(float)
t_s = (data_np[1:, 11]).astype(float)
print(np.argmin(t_s[t_s != 0]))
freq_col = (data_np[1:, 12]).astype(float)*(365.25/12) # collision frequency per month

# plot a histogram of t_s to see the distribution of synodic periods
plt.figure()
plt.hist(t_s[ttc != 0], bins = 100)
plt.xlabel('LCM Period (days)')
#plt.show()

t_s_no = t_s[t_s != 0]
t_s_no_s = t_s_no[t_s_no < 90]
print(len(t_s_no_s)/len(t_s)*100, '% of conjunctions have synodic period less than 90 days')
# Compute the distance from the center of the Earth to each point
r = np.sqrt(cnode_pos[:, 0]**2 + cnode_pos[:, 1]**2 + cnode_pos[:, 2]**2)

# Compute the latitude
latitude = np.arcsin(cnode_pos[:, 2] / r) * (180 / np.pi)  # Convert from radians to degrees

lat_dict = {}
for i in range(len(norad1)):
    if norad1[i] not in lat_dict:
        lat_dict[norad1[i]] = latitude[i]
    else:
        continue
for i in range(len(norad2)):
    if norad2[i] not in lat_dict:
        lat_dict[norad2[i]] = latitude[i]
    else:
        continue


# convert cnode_pos to cnode_pos_lla, and make a histogram of cnodes by latitude
# cnode_pos_lla = np.zeros((len(cnode_pos), 3))
# t_dt = dt.datetime(2021, 1, 1, 0, 0) # random date -- we don't care about longitude so it's unimportant
# for i in range(len(cnode_alt)):
#     cnode_pos_lla[i, :] = eci2lla(cnode_pos[i, 0], cnode_pos[i, 1], cnode_pos[i, 2], t_dt)
#     print(i/len(cnode_alt)*100, '% done')


# when including ETTC
# remove zero entries from ttc
idx_y_col = np.where(ttc != 0 )
# idx_y_col_sml = np.where((ttc != 0))
idx_y_col_sml = np.where((ttc != 0) & (ttc < t_s))
ttc_filt = ttc[idx_y_col_sml]
t_s_y_col = t_s[idx_y_col_sml]

# # when not including ETTC
# idx_y_col = np.arange(len(ttc))
# idx_y_col_sml = np.arange(len(ttc))
# ttc_filt = ttc[idx_y_col_sml]
# t_s_y_col = t_s[idx_y_col_sml]

print("latitude")
print(latitude[idx_y_col_sml])

plt_figs = True
if plt_figs == True:

    plt.figure()
    plt.hist(ttc_filt, bins = 100)
    # plt.yscale('log')
    plt.xlabel('Time to collision (s)')
    plt.ylabel('Number of conjunctions')
    plt.savefig('figures/1.pdf')
    plt.show()


    # t_s_y_col = t_s_y_col[t_s_y_col < 365] # we don't care about the conjunctions with large synodic periods -- they're in resonance.
    col_freq = 1/t_s_y_col*(365.25/12) # make this number of conjunctions per month
    # print(np.argmax(col_freq))

    print('Number of conjunctions with collisions and synodic period less than 90 days:', len(t_s_y_col))

    # make a 3d scatter plot of cnode_pos, where the color is determined by col_freq
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(cnode_pos[:, 0], cnode_pos[:, 1], cnode_pos[:, 2], s = 0.5, c='r')
    # plot a 3D sphere for the earth
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:500j]
    x = r_e*np.cos(u)*np.sin(v)
    y = r_e*np.sin(u)*np.sin(v)
    z = r_e*np.cos(v)
    ax.plot_surface(x, y, z, color = 'cornflowerblue', alpha = 0.5)
    # make the plot background white
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    plt.axis('equal')
    # plt.colorbar(sc, label = 'Conjunction frequency (per month)')
    ax.set_facecolor('white')
    # make the background of the plot white
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    # remove axis color
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    ax.grid(False)
    plt.savefig('figures/2.pdf')
    plt.tight_layout()
    plt.show()

    # make a histogram of conjunction frequency
    #cutoff = np.percentile(col_freq, 99)  # 99th percentile
    plt.figure(figsize=(4, 3))
    plt.hist(col_freq, bins=100, color='firebrick', alpha=0.5)
    # plt.xlim(0, 4)
    plt.yscale('log')
    plt.xlabel('Nodal conjunction frequency (monthly)')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.show()

    plt.figure(figsize = (4,3))
    plt.plot(latitude[idx_y_col_sml], col_freq, '.', markersize = 1, color = 'firebrick', alpha = 0.5)
    plt.xlabel('Latitude (deg)')
    plt.ylim(0,2)
    plt.ylabel('Nodal conjunction frequency (monthly)')
    plt.tight_layout()
    plt.savefig('figures/3.pdf')
    plt.show()

    plt.figure(figsize = (4, 1.5))
    plt.hist(latitude[idx_y_col_sml], bins=75, color = 'k', alpha = 0.3)
    plt.xlabel('Latitude [degrees]')
    # make x axis labels every 20 degrees
    plt.xticks(np.array([-80, -60, -40, -20, 0, 20, 40, 60, 80]))
    plt.xlim(-90, 90)
    plt.grid(axis = 'both', color = 'lightgray', linewidth = 0.25)
    # flip x axis direction
    plt.gca().invert_xaxis()
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig('figures/4.pdf')
    plt.show()

    # reapeat plot above but make it a 2d histogram
    # plt.hist2d(latitude[idx_y_col_sml], col_freq, bins = 100)
    # # plt.yscale('log')
    # plt.xlabel('Latitude (deg)')
    # plt.ylabel('Conjunctions per month')
    # plt.colorbar()
    # plt.show()


    # plot the histogram of cnode_pos_lla by latitude
    # plt.figure(figsize = (4,2))
    # plt.grid(axis = 'y', color = 'k')
    # # make sure grid is behind the bars
    # plt.gca().set_axisbelow(True)
    # plt.hist(latitude[idx_y_col_sml], bins=100, color = 'r')
    # plt.xlabel('Latitude (deg)')
    # plt.ylabel('Count')
    # # make all outside borders transparent
    # plt.gca().spines['top'].set_visible(False)
    # plt.gca().spines['right'].set_visible(False)
    # plt.gca().spines['left'].set_visible(False)
    # # make y axis on the right
    # plt.gca().yaxis.tick_right()
    # # same for label
    # plt.gca().yaxis.set_label_position("right")
    # plt.tight_layout()
    # plt.savefig('figures/5.png')
    # plt.show()

    plt.figure(figsize = (4,3))
    plt.plot(cnode_alt[idx_y_col_sml], col_freq, '.', markersize = 1, color = 'darkcyan', alpha = 0.5)
    plt.xlabel('Altitude (km)')
    plt.ylim(0,2)
    plt.ylabel('Conjunction frequency (monthly)')
    plt.tight_layout()
    plt.savefig('figures/6.pdf')
    plt.show()

    plt.figure(figsize = (4, 1.5))
    plt.hist(cnode_alt[idx_y_col_sml], bins=75, color = 'k', alpha = 0.3)
    plt.xlabel('Altitude [km]')
    # make x axis labels every 250 km
    plt.xticks(np.arange(250, 2000, 250))
    plt.xlim(250, 2000)
    plt.grid(axis = 'both', color = 'lightgray', linewidth = 0.25)
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig('figures/7.pdf')
    plt.show()

    # make imshow of cnode_alt[idx_y_col_sml] vs latitude[idx_y_col_sml] with count per bin as color
    # create a 2D histogram of latitude and altitude
    heatmap, xedges, yedges = np.histogram2d(latitude[idx_y_col_sml], cnode_alt[idx_y_col_sml], bins=75)
    # normalize the heatmap
    heatmap = np.nan_to_num(heatmap)  # Replace NaNs with 0
    # vmin, vmax = 0, 10
    # cmap = plt.cm.viridis
    # cmap_colors = cmap(np.linspace(0, 1, 256))
    # custom_cmap = mcolors.ListedColormap(cmap_colors)
    fig = plt.figure(figsize=(5, 3))
    heatmap_plot = plt.pcolormesh(yedges, xedges, heatmap, cmap=plt.cm.jet, shading='auto', norm = mcolors.LogNorm(vmin=1, vmax=heatmap.max()))
    cbar = plt.colorbar(heatmap_plot, label='Number of conjunction nodes', location = 'bottom')
    # make cbar sit on the bottom
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation=0)
    plt.xlim(250, 2000)  # Set x-axis limits for altitude
    plt.xticks(np.arange(250, 2000, 250))  # Set x-axis ticks for altitude
    # put colorbar on the bottom
    # add thin grid
    plt.grid(True, linestyle="-", alpha=0.25, color = 'gray')
    # Format ticks as integers
    plt.ylabel('Latitude [degrees]')
    plt.xlabel('Altitude [km]')
    # plt.title('Heatmap of Latitude vs Altitude - 01, 2025', fontsize=14)
    plt.tight_layout()
    # plt.grid(True, linestyle="--", alpha=0.5)
    # plt.savefig('figures/8.pdf')
    plt.show()


    # do the same for altitude
    # plot the histogram of cnode_pos_lla by latitude
    plt.figure(figsize = (4,2))
    plt.grid(axis = 'y', color = 'k')
    # make sure grid is behind the bars
    plt.gca().set_axisbelow(True)
    plt.hist(cnode_alt[idx_y_col_sml], bins=100, color = 'b')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Count')
    # make all outside borders transparent
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    # make y axis on the right
    plt.gca().yaxis.tick_right()
    # same for label
    plt.gca().yaxis.set_label_position("right")
    plt.tight_layout()
    plt.savefig('figures/8.pdf')
    plt.show()

    # # plot the CDF of ttc
    # plt.hist(col_freq, bins = 1000, color = 'k')
    # plt.yscale('log')
    # plt.xlabel('Conjunction frequency (per month)')
    # plt.ylabel('Count')
    # plt.show()

    # make ecdf of plot above
    # plt.figure(figsize = (4,3))
    # plt.ecdf(col_freq, color = 'k', linewidth = 2)
    # plt.grid()
    # plt.xlabel('Nodal Conjunction frequency (monthly)')
    # plt.ylabel('Cumulative Probability')
    # plt.ylim(0.8,1)
    # plt.grid(axis = 'both', color = 'gray', linewidth = 0.25)
    # plt.tight_layout()
    # plt.show()

    # plot the ecdf on col_freq, considering that the total number of satellites is 10854. 
    # add zeros to col_freq to make it 10854 long
    # col_freq_full = np.zeros(10854)
    # col_freq_full[:len(col_freq)] = col_freq
    # plt.figure(figsize = (3,4))
    # plt.ecdf(col_freq_full, color = 'k', linewidth = 2)
    # plt.grid(axis = 'both', color = 'gray', linewidth = 0.25)
    # plt.tight_layout()
    # plt.show()
else: 
    pass



# now, let's find out what percentage of conjunctions involve starlink
starlink_conj_name1 = np.char.find(name1[idx_y_col_sml], 'STARLINK') >= 0
starlink_conj_name2 = np.char.find(name2[idx_y_col_sml], 'STARLINK') >= 0
debris_conj_name1 = np.char.find(name1[idx_y_col_sml], 'DEB') >= 0
debris_conj_name2 = np.char.find(name2[idx_y_col_sml], 'DEB') >= 0

# what percentage of conjunctions involve at least one starlink?
starlink_conj = np.logical_or(starlink_conj_name1, starlink_conj_name2)
print(np.sum(starlink_conj)/len(starlink_conj)*100, '% of conjunctions involve at least one Starlink satellite')

# what percentage of conjunctions involve two starlinks?
starlink_conj_both = np.logical_and(starlink_conj_name1, starlink_conj_name2)
print(np.sum(starlink_conj_both)/len(starlink_conj)*100, '% of conjunctions involve two Starlink satellites')

# what percentage of conjunctions involve at least one debris object?
debris_conj = np.logical_or(debris_conj_name1, debris_conj_name2)
print(np.sum(debris_conj)/len(debris_conj)*100, '% of conjunctions involve at least one debris object')

# what percentage of conjunctions involve two debris objects?
debris_conj_both = np.logical_and(debris_conj_name1, debris_conj_name2)
print(np.sum(debris_conj_both)/len(debris_conj)*100, '% of conjunctions involve two debris objects')

# get a list of the unique norad ids across norad1 and norad2
norad = np.unique(np.concatenate((norad1[idx_y_col_sml], norad2[idx_y_col_sml])))
# get alt_range for each unique norad
alt_range = np.zeros((len(norad), 2))

# for each norad in norad, find the name associated with that norad from either norad1 or norad2
name = np.zeros(len(norad), dtype = object)
for i in range(len(norad)):
    idx = np.where(np.logical_or(norad1[idx_y_col_sml] == norad[i], norad2[idx_y_col_sml] == norad[i]))
    name[i] = name1[idx_y_col_sml][idx[0][0]] if norad1[idx_y_col_sml][idx[0][0]] == norad[i] else name2[idx_y_col_sml][idx[0][0]]

# get all the altitudes for conjunctions that involve starlink and plot the distribution
# alt_starlink = cnode_alt[np.logical_or(starlink_conj_name1, starlink_conj_name2)]

# for each norad, sum the number of conjunctions it is involved in
conj_count = np.zeros(len(norad))
sum_conj_freq = np.zeros(len(norad))
name_sec = []
alt_node = []
perc_sec_starlink = np.zeros(len(norad))
perc_sec_deb = np.zeros(len(norad))
for i in range(len(norad)):
    # get idxs of all conjunctions involving norad[i]
    idx_conj = np.logical_or(norad1[idx_y_col_sml] == norad[i], norad2[idx_y_col_sml] == norad[i])
    # count the number of conjunctions
    conj_count[i] = np.sum(np.logical_or(norad1[idx_y_col_sml] == norad[i], norad2[idx_y_col_sml] == norad[i]))  
    # get the name of the secondary object that is not norad[i]
    name_sec_item = [name2[idx_y_col_sml][j] if norad1[idx_y_col_sml][j] == norad[i] else name1[idx_y_col_sml][j] for j in np.where(idx_conj)[0]]
    # get the altitude of the conjunction nodes
    alt_node.append(cnode_alt[idx_y_col_sml][idx_conj])
    name_sec.append(name_sec_item)
    # perc_sec_starlink[i] = np.sum(np.char.find(name_sec_item, 'STARLINK') >= 0)/len(name_sec_item)*100
    # perc_sec_deb[i] = np.sum(np.char.find(name_sec_item, 'DEB') >= 0)/len(name_sec_item)*100
    # sum the col_freq for each conjunction
    sum_conj_freq[i] = np.sum(freq_col[idx_y_col_sml][idx_conj])
    c = 3

# rank norads/names based on sum_conj_freq
idx_rank = np.argsort(sum_conj_freq)[::-1]
norad_rank = norad[idx_rank]
name_rank = name[idx_rank]
sum_conj_freq_rank = sum_conj_freq[idx_rank]
conj_count_rank = conj_count[idx_rank]
# perc_sec_deb = perc_sec_deb[idx_rank]
# perc_sec_starlink = perc_sec_starlink[idx_rank]
name_sec_rank = [name_sec[i] for i in idx_rank]
alt_node_rank = [alt_node[i] for i in idx_rank]
alt_range_rank = np.concatenate((alt_range1, alt_range2), axis = 0)[idx_rank]


# print the top 10
# for i in range(10):
#     print('Norad:', norad_rank[-i], 'Conjunctions:', conj_count_rank[-i], 'Sum of conjunction frequency:', sum_conj_freq_rank[-i])

# print above in table format with column names and formatting
# Data for the table
table_data = []
for i in range(100):
    table_data.append([i + 1, name_rank[i], norad_rank[i], conj_count_rank[i],  sum_conj_freq_rank[i]])

# Column titles
headers = ["Rank", "Name", "Norad", "Conjunction Nodes", "Sat. Conjunction Frequency (monthly)"]

# Create the table
table = tabulate(table_data, headers=headers, tablefmt="pipe")

# Print the title and the table
print("**10 satellites with highest conjunction frequency**")
print(table)

# if norad id is 90000 or greater, extract the satellite's altitude and inclination and make a heat map of conjunction frequency across altitude and inclination

print("artificial objects")

# Data containers
altitudes = []
inclinations = []
latitudes = []
conjunction_frequencies = []

for i in range(len(norad)):
    if norad[i] in alt_range_dict and norad[i] in incl_dict and norad[i] in lat_dict:
        avg_altitude = np.mean(alt_range_dict[norad[i]])  # Compute average altitude
        
        # Store extracted data
        altitudes.append(avg_altitude)
        inclinations.append(incl_dict[norad[i]])
        latitudes.append(lat_dict[norad[i]])
        conjunction_frequencies.append(sum_conj_freq[i])

# Convert to NumPy arrays
altitudes = np.array(altitudes)
inclinations = np.array(inclinations)
latitudes = np.array(latitudes)
conjunction_frequencies = np.array(conjunction_frequencies)

# Define grid resolution for heatmap
alt_bins = np.linspace(300, 2000, 50)  # 50 altitude bins
inc_bins = np.linspace(min(inclinations), max(inclinations), 50)  # 50 inclination bins
lat_bins = np.linspace(-90, 90, 90)  # 50 latitude bins

# Compute 2D histogram (heatmap data)
heatmap, xedges, yedges = np.histogram2d(altitudes, latitudes, bins=[alt_bins, lat_bins], weights=conjunction_frequencies)

# Normalize heatmap values
heatmap = np.nan_to_num(heatmap)  # Replace NaNs with 0

vmin, vmax = 0, 10
cmap = plt.cm.viridis
cmap_colors = cmap(np.linspace(0, 1, 256))
custom_cmap = mcolors.ListedColormap(cmap_colors)

plt.figure(figsize=(5, 3))
heatmap_plot = plt.pcolormesh(alt_bins, lat_bins, heatmap.T, cmap=custom_cmap, shading='auto', vmin=vmin, vmax=vmax)
cbar = plt.colorbar(heatmap_plot, label='Conjunction Frequency')
cbar.ax.set_yticklabels([f"{int(tick)}" for tick in cbar.get_ticks()])  # Format ticks as integers

# Labels and title
plt.xlabel('Altitude [km]')
plt.ylabel('Latitude [degrees]')
# plt.title(f'Heatmap of Altitude vs Latitude - 01, 2025', fontsize=14)
plt.tight_layout()
# plt.grid(True, linestyle="--", alpha=0.5)
plt.show()

### **Plot Scatter Plot**
plt.figure(figsize=(9, 6))
sc = plt.scatter(altitudes, inclinations, c=conjunction_frequencies, cmap='plasma', edgecolors='k', s=40, alpha=0.8)
plt.colorbar(sc, label='Conjunction Frequency')
plt.xlabel('Altitude [km]')
plt.ylabel('Inclination [degrees]')
plt.title('Scatter Plot of Altitude vs Inclination', fontsize=14)
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()

print("end artificial objects")







# plot a histogram of the number of conjunctions each satellite is involved in
plt.hist(conj_count, bins = 100)
plt.yscale('log')
plt.xlabel('Number of conjunctions per satellite')
plt.ylabel('Number of satellites')
plt.savefig('figures/9.pdf')
plt.show()

sum_conj_freq_full = np.zeros(26802)
sum_conj_freq_full[:len(sum_conj_freq)] = sum_conj_freq
#np.savetxt('sum_conj_freq_full_2019.txt', sum_conj_freq_full)
#np.savetxt('sum_conj_freq_full_2022.txt', sum_conj_freq_full)

plt.figure(figsize = (3.5,4))
loaded_freq1 = np.loadtxt('sum_conj_freq_full_2019.txt')
loaded_freq2 = np.loadtxt('sum_conj_freq_full_2022.txt')
#plt.ecdf(loaded_freq1, color='purple', linewidth=2, label='January 2019')
#plt.ecdf(loaded_freq2, color='orange', linewidth=2, label='January 2022')
plt.ecdf(sum_conj_freq_full, color = 'darkcyan', linewidth = 2, label='January 2025')
plt.grid(axis = 'both', color = 'gray', linewidth = 0.25)
plt.axvline(x = 10, color = 'r', linestyle = '--', linewidth = 2)
plt.ylabel('Cumulative probability')
plt.xlabel('Sat. conj. freq. (monthly)')
plt.ylim(0,1)
plt.xlim(-0.25,5.5)
plt.legend()
plt.tight_layout()
plt.savefig('figures/10.pdf')
plt.show()

print("CDF Values")
print(len(sum_conj_freq_full))
print(sum_conj_freq_full)

plt.figure(figsize = (3.7,4))
#plt.ecdf(loaded_freq1, color='purple', linewidth=2, label='January 2019')
#plt.ecdf(loaded_freq2, color='orange', linewidth=2, label='January 2022')
plt.ecdf(sum_conj_freq_full, color = 'darkcyan', linewidth = 2, label='January 2025')
plt.grid(axis = 'both', color = 'gray', linewidth = 0.25)
plt.ylim(0.99,1)
plt.ylabel('Cumulative probability')
plt.xlabel('Sat. conj. freq. (monthly)')
plt.legend()
plt.tight_layout()
plt.savefig('figures/11.pdf')
plt.show()

# plot a histogram of the sum of the col_freq for each satellite
# plt.figure(figsize = (4,3))
# plt.ecdf(sum_conj_freq, linewidth = 3, color = 'firebrick')
# plt.xlabel('Satellite conjunction frequency (monthly)')
# plt.ylabel('Cumulative Probability')
# plt.ylim(0.95,1)
# plt.grid(axis = 'both', color = 'gray', linewidth = 0.25)
# plt.tight_layout()
# plt.show()

# print a table of the top 10 arrays in name_sec_rank. i.e. 1: [], 2:[]...
# Data for the table
num_obj_plot = 15
for i in range(num_obj_plot):
    print(i + 1, name_sec_rank[i])
# print a 

# plot top 10 by per and apogee
# r_p = np.array([304.6, 531.6, 388, 560, 383.8, 484.4, 573.6, 421.9, 528.9, 320.6])
# r_a = np.array([3119.1, 2031.8, 1721.7, 2913.4, 1719.4, 1716.1, 829.3, 1564.7, 561.7, 1602.9])
# r_p = [517.4, 543.7, 545.9, 391.4, 544.6, 403.8, 565.4, 581.7, 584.1, 553.8]
# r_a = [538.2, 559.4, 547.8, 504.4, 564.1, 568.8, 567, 660.1, 587.3, 561.3]

# r_p = alt_range_rank[:10, 0]
# r_a = alt_range_rank[:10, 1]

top_norads = norad_rank[:num_obj_plot]
r_p = []
r_a = []
for i in range(num_obj_plot):
    r_p.append(alt_range_dict[top_norads[i]][0])
    r_a.append(alt_range_dict[top_norads[i]][1])

# for 
# matplotlib.use('TkAgg')

plt.figure(figsize = (6,3))
for i in range(num_obj_plot):
    # plot a vertical line at x = i between y = per[i] and y = ap[i]
    plt.plot([i+1, i+1], [r_p[i], r_a[i]], 'k', linewidth = 3)
    # make a violin plot of the altitude of the conjunction nodes
    # plt.violinplot(alt_node_rank[i], positions = [i+1], showmedians = True, showextrema = False)
    plt.plot(np.ones(len(alt_node_rank[i]))*(i+1), alt_node_rank[i], 'x', markersize = 6, color = 'r', alpha = 0.6)
# plt.axhline(y = 550, color = 'firebrick', linestyle = '--', alpha = 0.5)
# plt.axhline(y = 1200, color = 'darkcyan', linestyle = '--', alpha = 0.5)
# make sure every 1 x is labeled    
plt.xticks(np.arange(num_obj_plot+1))
plt.xlabel('Rank')
# add a second x axis with the satellite names
# set text size for xtick labels to be smaller
mpl.rcParams.update({'font.size': 8})
plt.gca().set_xticklabels(name_rank[:num_obj_plot+1], rotation = 90)
plt.ylabel('Altitude (km)')
# plt.legend()
plt.tight_layout()
plt.savefig('figures/12.pdf')
plt.show()


    
