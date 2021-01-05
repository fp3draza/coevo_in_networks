# Load required libraries
require(dplyr)
require(ggfortify)
require(FactoMineR)
require(factoextra)

###############################################################
# Load data
data = read.csv('Results/processed_output/wrangled_data/complete_data.csv',row.names = 1)
# Remove NAs
data = data %>% na.omit()

##############################################################
# Run PCA at network level
network_PCA = PCA(data[,c("connectance","modularity_SA","number_of_modular_groups")],graph = FALSE)
# Add the network pca coordinates to the dataframe
data$network_pca_score = as.vector(network_PCA$ind$coord[,1])

##############################################################
# Run PCA at module level
module_PCA = PCA(data[,c("module_size","module_connectance")],graph = FALSE)
# Add the module pca coordiantes to the dataframe
data$module_pca_score = as.vector(module_PCA$ind$coord[,1])

##############################################################
# Run PCA at species level
species_PCA = PCA(data[,c("within_module_degree","among_module_connectivity")],graph = FALSE)
# Add the species pca coordinates to the dataframe
data$species_pca_score = as.vector(species_PCA$ind$coord[,1])

# Run PCA with all netowrk descriptors
overall_PCA =  PCA(data[,c("within_module_degree","among_module_connectivity",
                           "module_size","module_connectance",
                           "connectance","modularity_SA","number_of_modular_groups")],graph = FALSE)

# Add PCA coordinates to dataframe
data$overall_pca_score_dim1 = as.vector(overall_PCA$ind$coord[,1])
data$overall_pca_score_dim2 = as.vector(overall_PCA$ind$coord[,2])

#############################################################
# Save output
write.csv(data,'Results/processed_output/wrangled_data/complete_data_with_pca.csv')
