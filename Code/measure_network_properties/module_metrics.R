# This script measures and records the module descriptors of a set of networks
# Load required libraries
require(vegan)
require(dplyr)
require(tidyr)


measure_module_connectance = function(connectance_og,replica){
  # This function computes the connectance of the modules in the network. 
  # It requires the following parameters:
  # - connectance_og: target connectance value of the network
  # - replica: number of replica of the network
  
  # Load required libraries
  require(vegan)
  require(dplyr)
  require(igraph)
  
  # Read Membership data from MODULAR software output
  # It will currently read data from networks of 60 species with a 30-30 resource-consumer ratio. This is specified in the directory path (for_modular_40_60).
  # If you want to work with networks of different size, change the numbers (for_modular_X_X)
  membership = read.table(paste('Data/for_modular_40_60/resultsSA/MEMBERS_C_',connectance_og,'_R_',replica,'.txt',sep =''))
  membership = membership[-1,]
  membership$index = as.numeric(substring(membership$V1,2))
  membership = membership[order(membership$index),]
  membership_index = as.numeric(substring(membership$V1,2))
  membership_rows = subset(membership, substring(membership$V1,1,1) == 'R')
  membership_cols = subset(membership, substring(membership$V1,1,1) == 'C')
  
  # Load network file
  # It will currently read data from networks of 60 species with a 30-30 resource-consumer ratio. This is specified in the directory path (networks_30_30).
  # If you want to work with networks of different size, change the numbers (networks_X_X)
  network = read.csv(paste('Data/networks_30_30/C_',connectance_og,'_R_',replica,'.csv',sep = ''),header=FALSE)
  n_col = ncol(network)
  n_row = nrow(network)
  # Name rows and columns of dataframe
  colnames(network) = c(paste('C',1:n_col,sep=''))
  rownames(network) = c(paste('R',1:n_row,sep=''))
  
  # Initialize empty vectors for filling module descriptors
  module_vector = NULL
  module_size = NULL
  module_connectance = NULL
  
  # Cycle through modules
  for (module in unique(membership$V2)) {
    # Subset network by modules
    index_row = membership_rows %>% filter(V2 == module) %>% .$V1
    index_col = membership_cols %>% filter(V2 == module) %>% .$V1
    sub_net = network[as.character(index_row),as.character(index_col)]
    module_vector = c(module_vector,module)
    m = sum(dim(sub_net))
    
    # Compute connectance
    if (m == 0) {
      m = length(sub_net)
      module_size = c(module_size,m)
      connectance = sum(sub_net)/length(sub_net)
      module_connectance = c(module_connectance,connectance)
    }
    
    else if (m == 1) {
      m = 2
      module_size = c(module_size,m)
      connectance = 1
      module_connectance = c(module_connectance,connectance)
    }
    
    else if(m > 2){
      module_size = c(module_size,m)
      connectance = sum(sub_net)/(nrow(sub_net)*ncol(sub_net))
      module_connectance = c(module_connectance,connectance)
    }
  }
  
  # Extract number of modules
  num_modules = max(as.numeric(unique(membership$V2))) + 1
  network_name = paste('C_',connectance_og,'_R_',replica, sep = '')
  
  # Save results
  df = data.frame(network_name = network_name, module = as.numeric(module_vector), module_size = module_size,
                  module_connectance = module_connectance, number_of_modules = num_modules,
                  stringsAsFactors = FALSE)
  
  
  return(df)
}

measure_module_metrics = function(connectance_from,connectance_to,connectance_by,replica_from,replica_to,replica_by){
  # This function computes the module metrics for a specified set of networks.
  # It requires the following parameters:
  # - connectance_from: float, lowest connectance value
  # - connectance_to: float, highest connectance value
  # - connectance_by: float, increment in connectance
  # - replica_from: int, lowest replica value
  # - replica_to: int, highest replica value
  # - replica_by: int, increment in replica
  
  # Initialize empty dataframe
  data_frame_analysis = NULL
  # loop through connectance values
  for (connectance in seq(connectance_from,connectance_to,connectance_by)) {
    # loop through replcias
    for (replica in seq(replica_from,replica_to,replica_by)) {
      # Measure module properties
      current_result = measure_module_connectance(connectance,replica)
      # Append to  dataframe
      data_frame_analysis = rbind(data_frame_analysis, current_result)
    }
  }
  # Write data to csv 
  # It will currently save data for networks of 60 species with a 30-30 resource-consumer ratio. This is specified in the directory path (processed_output_30_30).
  # If you want to work with networks of different size, change the numbers (processed_output_X_X)
  write.csv(data_frame_analysis,paste("Results/processed_output_30_30/module_metrics/module_metrics.csv",sep=""))
}


###################################################################################################
# Run functions to measure module metrics for networks ranging in connectance from 
# 0.05 to 0.49 with increments of 0.02. Measure properties for all 20 replicas per network.
measure_module_metrics(0.05,0.49,0.02,1,20,1)


