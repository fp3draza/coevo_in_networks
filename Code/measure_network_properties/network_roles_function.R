# Load required libraries
require(igraph)
require(brainGraph)

compute_network_roles = function(connectance, replica){
  # This function computes the network roles of species
  # The funciton requires the following parameters:
	# - Connectance: float, the connectance value of the network
	# - Replica: int, the replica value of the network

  # Load network file
  # It will currently read data from networks of 60 species with a 30-30 resource-consumer ratio. This is specified in the directory path (networks_30_30).
  # If you want to work with networks of different size, change the numbers (networks_X_X)
  network_data = read.csv(paste('Data/networks_30_30/C_',connectance,'_R_',replica,'.csv',sep=''),header = FALSE)
 
  # Name columns and rows of network
  colnames(network_data) = paste('C',1:ncol(network_data),sep = '')
  rownames(network_data) = paste('R',1:nrow(network_data),sep = '')
  
  # Load membership file form MODULAR output
  # It will currently read data from networks of 60 species with a 30-30 resource-consumer ratio. This is specified in the directory path (for_modular_30_30).
  # If you want to work with networks of different size, change the numbers (for_modular_X_X)
  membership = read.table(paste('Data/for_modular_30_30/resultsSA/MEMBERS_C_',connectance,'_R_',replica,'.txt',sep =''))
  membership = membership[-1,]
  membership$index = as.numeric(substring(membership$V1,2))
  membership = membership[order(membership$index),]
  membership_results = membership_vector(membership,'Greys')
  membership_rows = membership_results[[2]]
  membership_cols = membership_results[[1]]

  # Reorder network file 
  network_data = network_data[as.character(membership_rows$V1),as.character(membership_cols$V1)]

  # Convert incidence matrix to graph 
  network = graph_from_incidence_matrix(network_data)

  # Define module membership
  mem = c(membership_rows$V2,membership_cols$V2)

  # Calculate the within module degree of species
  within_module_degree_z_score = within_module_deg_z_score(network,as.numeric(mem))

  # Calculate among-module connectivity
  among_module_connectivity = part_coeff(network,as.numeric(mem))

  # Define network metadata
  network_name = paste('C_',connectance,'_R_',replica, sep = '')
  type = 'niche'
  role = 'none'
  col = 'none'

  # Organize results in dataframe
  df_result = as.data.frame(cbind(network_name,type,names(among_module_connectivity),mem,within_module_degree_z_score,among_module_connectivity,role,col),stringsAsFactors = FALSE)
  df_result$mem = as.numeric(df_result$mem) - 1
  # Name dataframe columns
  colnames(df_result) = c('network_name','type','names','module','within_module_degree','among_module_connectivity','role','color')
  rownames(df_result) = NULL

  # Classify species depending on among-module connectivity and within module degree
  for (row in 1:nrow(df_result)) {
    if (is.nan(as.numeric(df_result$within_module_degree[row])) | is.na(as.numeric(df_result$within_module_degree[row]))) {
      df_result$role[row] = 'none'
      df_result$color[row] = '#######'
    }
    
    else{
    # Classify species as peripherals
      if(as.numeric(df_result$within_module_degree[row]) <= 2.5 & as.numeric(df_result$among_module_connectivity[row]) <= 0.62) {
        df_result$role[row] = 'peripheral'
        df_result$color[row] = '#CB2314'
      }
    
      # Classify species as connector
      else if (as.numeric(df_result$within_module_degree[row]) <= 2.5 & as.numeric(df_result$among_module_connectivity[row]) > 0.62) {
        df_result$role[row] = 'connector'
        df_result$color[row] = '#273046'
      }
  
     # Classify species as network hubs
      else if (as.numeric(df_result$within_module_degree[row]) > 2.5 & as.numeric(df_result$among_module_connectivity[row]) > 0.62) {
        df_result$role[row] = 'network hub'
        df_result$color[row] = '#1E1E1E'
      }
  
     # Classify species as module hubs
      else if (as.numeric(df_result$within_module_degree[row]) > 2.5 & as.numeric(df_result$among_module_connectivity[row]) <= 0.62) {
        df_result$role[row] = 'module hub'
        df_result$color[row] = '#354823'
      }
    }
  }

  # Return dataframe
  return(df_result)

}
