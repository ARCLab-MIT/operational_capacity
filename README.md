# operational_capacity
Framework for investigating the operationally feasible orbital carrying capacity. Assesses operational feasibility of the existing TLE catalog assuming Keplerian dynamics. 

The main file is tle_conjunctions_clean.py. This file takes in a set of TLEs and outputs a csv file with information about object pairs with conjunction nodes with likely collisions. 

The data folder contains tle files (both raw ones pulled from space-track.org and ones that have been cleaned to remove duplicates) for each year between 2000 and 2025 and each month between 01/2019 and 12/2024. If you want to use different tle files, you can use the import spacetrack.py file in the helper functions folder.

The tle_conjunctions_clean.py file also creates satellite parameter and conjunction parameter files to use. If these files exist already, it will just use those, otherwise it will create more. Some existing files that were generated for this project are included in the data folder. However, it is recommended to create your own conj_params and sat_params files for your own work.

plot_collision_results.py visualizes results from a single month or year and plot_historical_occupation plots results over time. These files can both be found in the visualization_code folder.

# Acknowledgments
Research was sponsored by the Department of the Air Force Artificial Intelligence Accelerator and was accomplished under Cooperative Agreement Number FA8750-19-2-1000. The views and conclusions contained in this document are those of the authors and should not be interpreted as representing the official policies, either expressed or implied, of the Department of the Air Force or the U.S. Government. The U.S. Government is authorized to reproduce and distribute reprints for Government purposes notwithstanding any copyright notation herein.
