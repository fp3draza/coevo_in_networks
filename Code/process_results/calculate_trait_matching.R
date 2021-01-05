# Load required packages
require(doSNOW)
require(parallel)
require(vegan)
require(dplyr)
require(tidyr)

compute_trait_matching_high_level = function(connectance,replica,mutualistic_selection,sim_replica){
 
  # This function computes the degree of trait matching in a network, module and species. It requires the following parameters:
	# - connectance : the degree of connectance of the network
	# - replica: the replica value of the network
	# - mutualistic_selection: the strength of mutualistic selection
   	# - sim_replica: the replica of the simulation

  # Load required libraries
  require(vegan)
  require(dplyr)

  # Output the current network being processed
  print(paste(connectance,"_",replica,"_",sim_replica,sep=''))

  # Read the simulation result for the current network with specific connectance, replica, strength of mutualistic selection and simulation replica
  result_data = read.csv(paste('Results/model_output/all_networks_m_',mutualistic_selection,'/C_',connectance,'_R_',replica,'/C_',connectance,'_R_',replica,'.csv_m_',mutualistic_selection,'_alpha02_phi05_theta0-10_sim',sim_replica,'.csv',sep = ''),row.names = 1)

  # Read the membership output of the MODULAR software for the given network
  membership = read.table(paste('Data/for_modular/resultsSA/MEMBERS_C_',connectance,'_R_',replica,'.txt',sep =''))
  membership = membership[-1,]
  membership$index = as.numeric(substring(membership$V1,2))
  membership = membership[order(membership$index),]
  membership_index = as.numeric(substring(membership$V1,2))
  membership_rows = subset(membership, substring(membership$V1,1,1) == 'R')
  membership_cols = subset(membership, substring(membership$V1,1,1) == 'C')
  mem_vec = as.numeric(c(membership_rows$V2,membership_cols$V2))

  # Define the network size 
  n_sp = as.numeric(dim(result_data)[2])
  n_row = n_sp / 2
  n_col = n_sp / 2

  # Define the final trait value of the species after coevolution
  z = as.numeric(result_data[3,])
  
  # Define the sensitivity parameter
  alpha = 0.2

  # Calculate the degree of trait matching between species of different guilds
  matching_diff = dist(z,upper=TRUE,diag=TRUE)
  matching = exp(-alpha*(matching_diff)^2)

  # Calculate the mean trait matching value for the network
  network_matching = mean(as.vector(matching))

  # Calculate an alternative metric of trait matching (THIS WAS NOT USED IN ANALYSES)
  network_high_matching = sum(as.vector(matching) > 0.8)/length(as.vector(matching))

  # Convert into dataframe
  matching = data.frame(as.matrix(matching))

  # Compute the degree of trait matching for each species
  species_matching = rowMeans(matching)
  # This is an alternative metric which was not used in analyses
  species_high_matching = rowSums(matching > 0.8)/60

  # Name columns and rows of the dataframe
  colnames(matching) = c(paste('R',1:30,sep=''),paste('C',1:30,sep=''))
  rownames(matching) = c(paste('R',1:30,sep=''),paste('C',1:30,sep=''))

  # Initialize empty vector for analysing trait matching in modules
  module_vector = NULL
  module_matching_vector = NULL
  module_high_matching_vector = NULL

  # Cycle through the modules in the network
  for (module in unique(membership$V2)) {
    # Extract the module membership
    index = membership %>% filter(V2 == module) %>% .$V1
    # Subset species by module membership
    sub_matching = matching[as.character(index),as.character(index)]
    module_vector = c(module_vector,module)
    # Calculate degree of trait matching in the current module
    module_matching_vector = c(module_matching_vector,mean(as.vector(as.matrix(sub_matching))))
    # Alternative metric of trait matching at the module level, this was not used in analyses
    module_high_matching_vector = c(module_high_matching_vector,sum(as.vector(as.matrix(sub_matching) > 0.8))/length(as.vector(as.matrix(sub_matching))))
  }

  # Create dataframe containing the module membership and trait matching at the module
  module_df = data.frame(as.numeric(module_vector),module_matching_vector,module_high_matching_vector,stringsAsFactors = FALSE)
  # Name dataframe
  module_df = module_df %>% rename(module = as.numeric.module_vector.)
  col_names = c(paste('R',1:30,sep = ''),paste('C',1:30,sep = ''))
  network_name = paste('C_',connectance,'_R_',replica, sep = '')
 
  # Store additional information 
  m_val = mutualistic_selection
  type = 'niche'
  module = as.numeric(mem_vec)

  # Create empty dataframe to store all trait matching information for the current network
  df = data.frame(network_name = network_name, mutualistic_selection = m_val, type = type, names = col_names,module = module, network_matching = network_matching,
                  network_high_matching = network_high_matching, species_matching = species_matching,
                  species_high_matching = species_high_matching,stringsAsFactors = FALSE)

  # Merge dataframes into a single one
  df = left_join(df,module_df %>% mutate(module = module + 1))

  # Return datafarme containing trait matching information for current network
  return(df)
}


run_tm_high_analysis = function(m,connectance_from,connectance_to,connectance_by,replica_from,replica_to,replica_by,replica_result_from,replica_result_to,replica_result_by,results_directory){
 
 # This function computes the trait matching metrics for all netowrks. It requires the following parameters:
	# - m: strength of mutualistic selection
 	# - connectance_from : min connectance
	# - connectance_to: max connectance
	# - connectance_by: increment of connectance by
	# - replica_from: min replica
	# - replica_to: max replica
	# - replica_by: increment of replica by
	# - replica_result_from: min simulation replica
	# - replica_result_to: max simulation replica
	# - replica_result_by: increment simulation replica by
	# - results_directory: directory where results are stored

  # Create empty dataframe
  data_frame_analysis = NULL

  # Cycle through all connectance values
  for (connectance in seq(connectance_from,connectance_to,connectance_by)) {
    # Cycle through all replica values
    for (replica in seq(replica_from,replica_to,replica_by)) {
      # Cycle through all simulation replica values
      for(replica_result in seq(replica_result_from,replica_result_to,replica_result_by)){
        # Compute trait matching results of current network
	current_result = compute_trait_matching_high_level(connectance,replica,m,replica_result)
	# Add trait matching results to dataframe
        data_frame_analysis = rbind(data_frame_analysis, current_result)
      }
    }
  }
  # Store dataframe with all trait matching results
  write.csv(data_frame_analysis,paste("Results/processed_output/trait_matching/",m,"_r_sq_trait_matching_analysis_niche_no_noise.csv",sep=""))
}

tm_RUN = function(m){

 # This function computes trait matching in all networks across the gradient of mutualistic selection strength.
 # The function requires the following parameters:
  	# - m : mutualistic selection strength

  # Define parameters to pass to function
  connectance_from = 0.05
  connectance_to = 0.49
  connectance_by = 0.02
  replica_from = 1
  replica_to = 20
  replica_by = 1
  replica_result_from = 1
  replica_result_to = 20
  replica_result_by = 1
  results_directory = 'Results/model_output/'

  # Run trait matching analysis for the specified parameters
  run_tm_high_analysis(m,connectance_from,connectance_to,connectance_by,replica_from,replica_to,replica_by,replica_result_from,replica_result_to,replica_result_by,results_directory)

}

###################################################################################################
# RUN CODE

# Define the directory where Results should be stored
results_directory = 'Results/model_output/'
# List the directories contained in the Results directory
directories = list.dirs(results_directory,recursive = FALSE)
# Extract the different mutualistic selection levels
m_vals = gsub("[^0-9.]", "",  directories)

# This code runs the trait matching analyis in parallel
# Detect number of cores available
cores = detectCores()
# Create cluster with the number of cores
cl =  makeCluster(cores)
registerDoSNOW(cl)
# For each mutualistic selection value calculate trait matching analysis for all networks
foreach(m=m_vals) %dopar% tm_RUN(m)
# Stop cluster
stopCluster(cl)


