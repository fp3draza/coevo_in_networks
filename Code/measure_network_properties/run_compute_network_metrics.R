# This script runs the functions to compute network metrics
# Source required functions
source("Code/measure_network_properties/compute_network_metrics_functions.R")

# Specify parameters
# Define folder where networks are stored
network_folder = "Data/networks/"
# Define folder where output from MODULAR software is stored
modular_file = "Data/for_modular/resultsSA/OUT_MOD.txt"
# Define where output should be stored
save_folder =  "Results/processed_output/network_metrics/"
# Define network size
num_sp_in_A = 30
num_sp_in_B = 30

# Run function with specified parameters
build_dataframe_on_nested_modular_metrics(network_folder,modular_file,save_folder,num_sp_in_A,num_sp_in_B)

