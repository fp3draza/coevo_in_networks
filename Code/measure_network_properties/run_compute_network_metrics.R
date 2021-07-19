# This script runs the functions to compute network metrics
# Source required functions
source("Code/measure_network_properties/compute_network_metrics_functions.R")

# Specify parameters
# Define folder where networks are stored
# It will currently read data from networks of 60 species with a 30-30 resource-consumer ratio. This is specified in the directory path (networks_30_30).
# If you want to work with networks of different size, change the numbers (networks_X_X)
network_folder = "~/Documents/PhD/Work/Projects/coevo_in_networks/Data/networks_30_30/"
# Define folder where output from MODULAR software is stored
# It will currently read data from networks of 60 species with a 30-30 resource-consumer ratio. This is specified in the directory path (for_modular_30_30).
# If you want to work with networks of different size, change the numbers (for_modular_X_X)
modular_file = "~/Documents/PhD/Work/Projects/coevo_in_networks/Data/for_modular_30_30/resultsSA/OUT_MOD.txt"
# Define where output should be stored
# It will currently store data from networks of 60 species with a 30-30 resource-consumer ratio. This is specified in the directory path (processed_output_30_30).
# If you want to work with networks of different size, change the numbers (processed_output_X_X)
save_folder =  "~/Documents/PhD/Work/Projects/coevo_in_networks/Results/processed_output_30_30/network_metrics/"
# Define network size
# Change these numbers if you want to work with networks of other size
num_sp_in_A = 30
num_sp_in_B = 30

# Run function with specified parameters
build_dataframe_on_nested_modular_metrics(network_folder,modular_file,save_folder,num_sp_in_A,num_sp_in_B)

