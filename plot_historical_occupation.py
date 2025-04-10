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
import pandas as pd
import matplotlib.dates as mdates

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

def count_population(file_path):
    count = 0

    with open(file_path, 'r') as file:
        lines = file.readlines()

    for i in range(0, len(lines), 3):
        if i + 2 >= len(lines):
            break  # Avoid index error if file isn't perfectly divisible by 3

        name_line = lines[i].strip()
        line1 = lines[i+1].strip()
        line2 = lines[i+2].strip()

        try:

            norad_number = int(line1[2:7])

            inclination = float(line2[8:16])

            if norad_number >= 90000:
                mean_motion = float(line2[46:54])
            else:
                mean_motion = float(line2[52:63])

            if 49.99 <= inclination <= 55.01 and 15.022 <= mean_motion <= 15.088:
                print(f"Match: {name_line}")  # Print satellite name
                count += 1

        except ValueError:
            continue  # Skip malformed lines

    return count

annual_occupation_550 = []
annual_occupation_800 = []
annual_occupation_900 = []
annual_population_550 = []
total_objects = []
start_year = 2000
end_year = 2025
for year in range(start_year, end_year):
    for month in range(1, 2):
        # data = read_csv("data/collision_results_wp_0p1_1.csv")  
        data = read_csv('data/collision_results_annual_prop/collision_results_' + str(year) + '.csv')
        print("year: ", year, " month: ", month)
        #data = read_csv(f"data/collision_results_monthly_prop/collision_results_{month:02d}{year}.csv")


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
        #print(np.argmin(t_s[t_s != 0]))
        freq_col = (data_np[1:, 12]).astype(float)*(365.25/12) # collision frequency per month

        t_s_no = t_s[t_s != 0]
        t_s_no_s = t_s_no[t_s_no < 90]

        # Compute the distance from the center of the Earth to each point
        r = np.sqrt(cnode_pos[:, 0]**2 + cnode_pos[:, 1]**2 + cnode_pos[:, 2]**2)

        # Compute the latitude
        latitude = np.arcsin(cnode_pos[:, 2] / r) * (180 / np.pi)  # Convert from radians to degrees

        # # when including ETTC
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

        # now, let's find out what percentage of conjunctions involve starlink
        starlink_conj_name1 = np.char.find(name1[idx_y_col_sml], 'STARLINK') >= 0
        starlink_conj_name2 = np.char.find(name2[idx_y_col_sml], 'STARLINK') >= 0
        debris_conj_name1 = np.char.find(name1[idx_y_col_sml], 'DEB') >= 0
        debris_conj_name2 = np.char.find(name2[idx_y_col_sml], 'DEB') >= 0

        # what percentage of conjunctions involve at least one starlink?
        starlink_conj = np.logical_or(starlink_conj_name1, starlink_conj_name2)
        #print(np.sum(starlink_conj)/len(starlink_conj)*100, '% of conjunctions involve at least one Starlink satellite')

        # what percentage of conjunctions involve two starlinks?
        starlink_conj_both = np.logical_and(starlink_conj_name1, starlink_conj_name2)
        #print(np.sum(starlink_conj_both)/len(starlink_conj)*100, '% of conjunctions involve two Starlink satellites')

        # what percentage of conjunctions involve at least one debris object?
        debris_conj = np.logical_or(debris_conj_name1, debris_conj_name2)
        #print(np.sum(debris_conj)/len(debris_conj)*100, '% of conjunctions involve at least one debris object')

        # what percentage of conjunctions involve two debris objects?
        debris_conj_both = np.logical_and(debris_conj_name1, debris_conj_name2)
        #print(np.sum(debris_conj_both)/len(debris_conj)*100, '% of conjunctions involve two debris objects')

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
            sum_conj_freq[i] = np.sum(freq_col[idx_y_col_sml][idx_conj])
            c = 3

        # rank norads/names based on sum_conj_freq
        idx_rank = np.argsort(sum_conj_freq)[::-1]
        norad_rank = norad[idx_rank]
        name_rank = name[idx_rank]
        sum_conj_freq_rank = sum_conj_freq[idx_rank]
        conj_count_rank = conj_count[idx_rank]
        name_sec_rank = [name_sec[i] for i in idx_rank]
        alt_node_rank = [alt_node[i] for i in idx_rank]
        alt_range_rank = np.concatenate((alt_range1, alt_range2), axis = 0)[idx_rank]


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
        #print("**10 satellites with highest conjunction frequency**")
        #print(table)

        collision_550 = False
        collision_800 = False
        collision_900 = False

        for i in range(len(norad)):
            if norad[i] == 91906:
                collision_550 = True
                a = i
            if norad[i] == 92615:
                collision_800 = True
                b = i
            if norad[i] == 93480:
                collision_900 = True
                c = i
            
        if collision_550 == True:
            #alt_avg = (sum_conj_freq[a-1] + sum_conj_freq[a] + sum_conj_freq[a+1])/3
            #annual_occupation_550.append(alt_avg)
            annual_occupation_550.append(sum_conj_freq[a])
        else:
            annual_occupation_550.append(0)
        if collision_800 == True:
            #incl_avg = (sum_conj_freq[b-171] + sum_conj_freq[b] + sum_conj_freq[b+171])/3
            #annual_occupation_800.append(incl_avg)
            annual_occupation_800.append(sum_conj_freq[b])
        else:
            annual_occupation_800.append(0)
        if collision_900 == True:
            annual_occupation_900.append(sum_conj_freq[c])
        else:
            annual_occupation_900.append(0)

        def count_lines(filename):
            with open(filename, "r") as file:
                return sum(1 for line in file)

        # Example usage
        filename = "data/batch_tles_annual/" + str(year) + ".txt"
        #filename = f"data/batch_tles_monthly/{month:02d}{year}.txt"
        num_lines = count_lines(filename)
        total_objects.append(num_lines/3)

        count = count_population(filename)
        annual_population_550.append(count)

        # if norad id is 90000 or greater, extract the satellite's altitude and inclination and make a heat map of conjunction frequency across altitude and inclination
        print("artificial objects")

        # Data containers
        altitudes = []
        inclinations = []
        conjunction_frequencies = []

        # Extract satellite data for NORAD IDs ≥ 90000
        for i in range(len(norad)):
            if norad[i] >= 90000 and norad[i] in alt_range_dict and norad[i] in incl_dict:
                avg_altitude = np.mean(alt_range_dict[norad[i]])  # Compute average altitude
                rounded_inclination = np.round(incl_dict[norad[i]])  # Get inclination

                # Store extracted data
                altitudes.append(avg_altitude)
                inclinations.append(rounded_inclination)
                #inclinations.append(incl_dict[norad[i]])
                conjunction_frequencies.append(sum_conj_freq[i])

        # Convert to NumPy arrays
        altitudes = np.array(altitudes)
        inclinations = np.array(inclinations)
        conjunction_frequencies = np.array(conjunction_frequencies)

        # Define grid resolution for heatmap
        alt_bins = np.arange(300, 2010, 10)  # altitude bins
        inc_bins = np.arange(0, 115.1, 5)  # inclination bins

        # Compute 2D histogram (heatmap data)
        heatmap, xedges, yedges = np.histogram2d(altitudes, inclinations, bins=[alt_bins, inc_bins], weights=conjunction_frequencies)
        hist_counts, _, _ = np.histogram2d(altitudes, inclinations, bins=[alt_bins, inc_bins])
        print("Histogram Counts per Bin:", hist_counts.sum(axis=0))  # Sum along altitude axis to check inclination bins
        print("Unique inclinations in dataset:", np.unique(inclinations))

        # Normalize heatmap values
        heatmap = np.nan_to_num(heatmap)  # Replace NaNs with 0

        # Create a custom colormap where 0 is a light blue and the rest follow 'plasma'
        vmin, vmax = 0, 10
        cmap = plt.cm.plasma
        cmap_colors = cmap(np.linspace(0, 1, 256))
        #cmap_colors[0] = np.array([173/255, 216/255, 230/255, 1])  # Light blue for 0 conjunctions
        custom_cmap = mcolors.ListedColormap(cmap_colors)

    # #Plot Heatmap
    #     plt.figure(figsize=(9, 6))
    #     heatmap_plot = plt.pcolormesh(alt_bins, inc_bins, heatmap.T, cmap=custom_cmap, shading='auto', vmin=vmin, vmax=vmax)
    #     cbar = plt.colorbar(heatmap_plot, label='Conjunction Frequency')
    #     cbar.ax.set_yticklabels([f"{int(tick)}" for tick in cbar.get_ticks()])  # Format ticks as integers

    #     # Labels and title
    #     plt.xlabel('Altitude (km)', fontsize=12)
    #     plt.ylabel('Inclination (degrees)', fontsize=12)
    #     #plt.title(f'Heatmap of Altitude vs Inclination - {month:02d},{year}', fontsize=14)
    #     plt.title(f'Heatmap of Altitude vs Inclination - {year}', fontsize=14)
    #     plt.grid(True, linestyle="--", alpha=0.5)
    #     plt.show()

        print("end artificial objects")

        # print a table of the top 10 arrays in name_sec_rank. i.e. 1: [], 2:[]...
        # Data for the table
        num_obj_plot = 15
        for i in range(num_obj_plot):
            (i + 1, name_sec_rank[i])

        top_norads = norad_rank[:num_obj_plot]
        r_p = []
        r_a = []
        for i in range(num_obj_plot):
            r_p.append(alt_range_dict[top_norads[i]][0])
            r_a.append(alt_range_dict[top_norads[i]][1])

plot_months = pd.date_range(start=f'{start_year}-01-01', end=f'{end_year-1}-12-01', freq='MS')
plot_years = pd.date_range(start=f'{start_year}-01-01', end=f'{end_year-1}-01-01', freq='YS')

# Your data must match the length of plot_months (72 months)
# annual_occupation_550 = [...]
# annual_population_550 = [...]

print("annual_occupation_550: ", annual_occupation_550)
print("annual_population_550: ", annual_population_550)

fig, ax1 = plt.subplots(figsize=(12, 6))

# Plot first data series (left y-axis)
color1 = 'tab:blue'
ax1.set_xlabel('Date')
ax1.set_ylabel('Conjunction Frequency at 550 km', color=color1)
#ax1.plot(plot_months, annual_occupation_550, color=color1, label='Conjunction Frequency')
ax1.plot(plot_years, annual_occupation_550, color=color1, label='Conjunction Frequency')
ax1.tick_params(axis='y', labelcolor=color1)

# Plot second data series (right y-axis)
ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel('Number of Objects at 550 km', color=color2)
#ax2.plot(plot_months, annual_population_550, color=color2, label='Number of Objects')
ax2.plot(plot_years, annual_population_550, color=color2, label='Number of Objects')
ax2.tick_params(axis='y', labelcolor=color2)

# Format x-axis
#ax1.set_xlim([pd.Timestamp('2019-01-01'), pd.Timestamp('2024-12-31')])
ax1.set_xlim([pd.Timestamp(f'{start_year}-01-01'), pd.Timestamp(f'{end_year-1}-1-01')])
#ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=12))
ax1.xaxis.set_major_locator(mdates.YearLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

#ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))


# Final plot formatting
plt.title('Monthly Conjunction Frequency vs Number of Objects at 550 km and 53° Inclination')
fig.autofmt_xdate()
fig.tight_layout()
ax1.grid(True)
plt.show()


# plot_years = list(range(start_year, end_year))
# plt.figure()
# plt.plot(plot_years, annual_occupation_800, label='800 km')
# plt.xlabel('Year')
# plt.ylabel('Conjunction Frequency')
# plt.title('Conjunction Frequency at 800 km (Location of Iridium-Cosmos Collision in 2009)')
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# plot_years = list(range(start_year, end_year))
# plt.figure()
# plt.plot(plot_years, annual_occupation_900, label='900 km')
# plt.xlabel('Year')
# plt.ylabel('Conjunction frequency')
# plt.title('Conjunction Frequency at 900 km (Location of Chinese AntiSat Test in 2007)')
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# plt.figure()
# plt.plot(plot_years, total_objects)
# plt.xlabel('Year')
# plt.ylabel('Number of Tracked Objects')
# plt.title('Total Number of Tracked Objects per Year')
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()
