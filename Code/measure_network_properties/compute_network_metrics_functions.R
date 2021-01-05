############################################################################################
# Compute network metrics
# This script contains the code necessary to measure network descriptors

############################################################################################
# First load required functions and libraries
# Source the function to measure nestedness
source("Code/measure_network_properties/nestedness_function.R")
# Load library
require(igraph)

############################################################################################
# To measure modularity and determine modules and membership, first run the MODULAR software on each networks
# with the following configuration:
	# Networks: Bipartite
	# Input File: Matrix
	# Metric to optimize: Newman & Girvan 2004
	# Method: Simulated Annealing

extract_modularity_from_network = function(modularity_file,network_file,metric){
  # This function extracts the modularity results for each network after running the MODULAR software.
  # It requires the following parameters:
	# - modularity_file: file containing the output form the MODULAR software
	# - network_file: file containing the network in question
	# - metric: string, metric to extract from modularity results


  # Read modularity details
  modularity_table = read.table(modularity_file, header=TRUE)
  modularity_table = modularity_table[unique(modularity_table$File),]

  # Network to read
  network_name = tools::file_path_sans_ext(network_file)
  network_name_txt = paste(network_name,".txt",sep="")

  # Extract modularity values
  return_variable = modularity_table[modularity_table$File == network_name_txt,metric]

  # Return
  return(return_variable)
}

############################################################################################
# Finally, I build a dataframe where all network metrics will be stored. The function
# build_dataframe_on_nested_modular_metrics takes as parameters:
# network_folder: a folder to look for the built networks
# modular_file: a file where the results from the MODULAR software are stored
# save_folder: a place to store the dataframe
# num_sp_in_A : number of species in group A
# num_sp_in_B : number of species in group B

build_dataframe_on_nested_modular_metrics = function(network_folder,modular_file,save_folder,num_sp_in_A,num_sp_in_B){
   # This function fills a dataframe where all network metrics are stored. 
   # The function requires the following parameters:
	# - network_folder: folder where networks are stored
	# - modular_file: output file from the MODULAR software
	# - save_folder: folder where data frame will be stored
	# - num_sp_in_A: number of species in guild A
	# - num_sp_in_B: number of species in guild B

  # Intialize dataframe
  data_frame_for_nested_modular_metrics = data.frame(network_full_name = character(),
                                                     connectance = double(), type = character(), number_species = double(), number_species_A = double(),
                                                     number_species_B = double(), row_nestedness = double(), column_nestedness = double(),
                                                     nestedness = double(), modularity_SA = double(), number_of_modular_groups = double(), stringsAsFactors = FALSE)

  # List all network files
  list_of_networks = list.files(network_folder)

  # Loop through networks
  for (network in list_of_networks) {

    # Extract basic information
    network_full_name = tools::file_path_sans_ext(network)
    connectance = connectance(paste(network_folder,network,sep=""))
    replica = as.numeric(strsplit(network_full_name,"_")[[1]][4])
    number_species = num_sp_in_A + num_sp_in_B
    number_species_A = num_sp_in_A
    number_species_B = num_sp_in_B

    # Calculate nestedness
    print(paste(network_folder,network,sep=""))
    nested = nestedness_s(paste(network_folder,network,sep=""))
    row_nestedness = nested[[1]]
    column_nestedness = nested[[2]]
    nestedness = nested[[3]]

    # Calculate modularity and number of modular groups
    modularity_results = modularity_and_roles(paste(network_folder,network,sep=""))
    modularity_SA = modularity_results[[2]]
    number_of_modular_groups = max(modularity_results[[1]]$module)

    # Type of network
    if (grepl("random", network, fixed = TRUE)) {
      type = "random"
    }

    else{
      type = 'niche'
    }


    # Fill current row of the dataframe
    data_frame_for_nested_modular_metrics[nrow(data_frame_for_nested_modular_metrics)+1,] = c(network_full_name,connectance,
                                                                                            type, number_species, number_species_A,
                                                                                              number_species_B, row_nestedness, column_nestedness,
                                                                                              nestedness, modularity_SA,number_of_modular_groups)
  }
  # Save data frame as csv file
  write.csv(data_frame_for_nested_modular_metrics,paste(save_folder,"/modular_nested_metrics.csv",sep=""))

}
