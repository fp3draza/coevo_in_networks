# This scripts joins all processed output results into a single dataframe
require(dplyr)
# Read data
# This piece of code will read in the results of networks consisting of 60 species with 30-30 resource and consumer species. 
# This is specified by setting processed_output_30_30. If you want to work with the results of different network sizes/ratios,
# change the numbers (processed_output_X_X).
# Read trait matching results
trait_matching_data = read.csv('~/Documents/PhD/Work/Projects/coevo_in_networks/Results/processed_output_30_30/trait_matching/merged_results.csv',row.names = 1)
# Read species role data
species_role_data = read.csv('~/Documents/PhD/Work/Projects/coevo_in_networks/Results/processed_output_30_30/network_roles/network_roles_niche.csv',row.names = 1)
# Read network properties data
network_data = read.csv('~/Documents/PhD/Work/Projects/coevo_in_networks/Results/processed_output_30_30/network_metrics/modular_nested_metrics.csv',row.names = 1)
# Read module properties data
module_data = read.csv('~/Documents/PhD/Work/Projects/coevo_in_networks/Results/processed_output_30_30/module_metrics/module_metrics.csv',row.names = 1)

# Calculate standardized trait matching values

# Zscore for networks
trait_matching_data = trait_matching_data %>%
  group_by(mutualistic_selection, network_name) %>%
  mutate(mean_network_matching = mean(network_matching,na.rm = TRUE),
         mean_network_matching_initial = mean(network_matching_initial, na.rm = TRUE),
         mean_network_matching_theta = mean(network_matching_theta, na.rm = TRUE))

trait_matching_data = trait_matching_data %>%
  group_by(mutualistic_selection) %>%
  mutate(z_score_network = (mean_network_matching - mean(mean_network_matching))/sd(mean_network_matching),
         z_score_network_initial = (mean_network_matching_initial - mean(mean_network_matching_initial))/sd(mean_network_matching_initial),
         z_score_network_theta = (mean_network_matching_theta - mean(mean_network_matching_theta))/sd(mean_network_matching_theta)) %>%
  mutate(z_score_network_level = if_else(z_score_network > 1, 'high',if_else(z_score_network < -1,'low','average')))

# Zscore for modules
trait_matching_data = trait_matching_data %>%
  group_by(mutualistic_selection,network_name,module) %>%
  mutate(mean_module_matching = mean(module_matching_vector,na.rm = TRUE),
         mean_module_matching_initial = mean(module_matching_vector_initial, na.rm = TRUE),
         mean_module_matching_theta = mean(module_matching_vector_theta, na.rm = TRUE))

trait_matching_data = trait_matching_data %>%
  group_by(mutualistic_selection) %>%
  mutate(z_score_modules = (mean_module_matching - mean(mean_module_matching))/sd(mean_module_matching),
         z_score_modules_initial = (mean_module_matching_initial - mean(mean_module_matching_initial))/sd(mean_module_matching_initial),
         z_score_modules_theta = (mean_module_matching_theta - mean(mean_module_matching_theta))/sd(mean_module_matching_theta)) %>%
  mutate(z_score_modules_level = if_else(z_score_modules > 1, 'high',if_else(z_score_modules < -1,'low','average')))

# Zscore for species
trait_matching_data = trait_matching_data %>%
  group_by(mutualistic_selection, network_name, module, names) %>%
  mutate(mean_species_matching = mean(species_matching, na.rm =  TRUE),
         mean_species_matching_initial = mean(species_matching_initial, na.rm = TRUE),
         mean_species_matching_theta = mean(species_matching_theta, na.rm = TRUE))

trait_matching_data = trait_matching_data %>%
  group_by(mutualistic_selection) %>%
  mutate(z_score_species = (mean_species_matching - mean(mean_species_matching))/sd(mean_species_matching),
         z_score_species_initial = (mean_species_matching_initial - mean(mean_species_matching_initial))/sd(mean_species_matching_initial),
         z_score_species_theta = (mean_species_matching_theta - mean(mean_species_matching_theta))/sd(mean_species_matching_theta)) %>%
  mutate(z_score_species_level = if_else(z_score_species > 1, 'high',if_else(z_score_species < -1,'low','average')))

# Merge dataframes

# Rename column from species_role_data
species_role_data = species_role_data %>% mutate(module = module + 1)

# Merge species_role_data with trait_matching_data
merged_data = left_join(trait_matching_data,species_role_data)

# Rename column from network_data
network_data = network_data %>% rename(network_name = network_full_name)

# Merge network_data with merged_data
merged_data = left_join(merged_data,network_data)

# Merge module_data with merged_data
#module_data = module_data %>% mutate(module = module + 1)
complete_data = left_join(merged_data,module_data)

# Remove dataframes no longer needed
rm(merged_data, module_data, network_data, species_role_data, trait_matching_data)

# Create separate files for data
data_for_network_pca <- complete_data %>% select(network_name, mutualistic_selection, names, module, network_matching, network_matching_initial, network_matching_theta,
                                                 connectance, number_species, number_species_A, number_species_B,
                                                 nestedness, modularity_SA, number_of_modular_groups)

data_for_network_pca <- data_for_network_pca %>% distinct(network_name,mutualistic_selection, .keep_all = TRUE)


data_for_module_pca <- complete_data %>% select(network_name, mutualistic_selection, names, module, module_matching_vector, module_matching_vector_initial, module_matching_vector_theta,
                                                module_size, number_species, number_species_A, number_species_B,
                                                module_connectance)

data_for_module_pca <- data_for_module_pca %>% distinct(network_name,mutualistic_selection,module, .keep_all = TRUE)

data_for_species_pca <- complete_data %>% select(network_name, mutualistic_selection, names, module, mean_species_matching, 
                                                 within_module_degree, number_species, number_species_A, number_species_B,
                                                 among_module_connectivity,z_score_species)

data_for_overall_pca <- complete_data %>% select(network_name, mutualistic_selection, names, module, mean_species_matching, mean_species_matching_initial, mean_species_matching_theta,
                                                 within_module_degree, number_species, number_species_A, number_species_B,
                                                 among_module_connectivity,z_score_species,module_size,module_connectance,
                                                 connectance,nestedness,modularity_SA,number_of_modular_groups)

# Save merged dataframe
# This piece of code will store the results from the networks consisting of 60 species with 30-30 resource and consumer species. 
# The results are stored in a directory that contains the information of number and ratio of species in the networks (processed_output_30_30).
# If you work with networks of different size/ratios, change the numbers (processed_output_X_X)
write.csv(data_for_network_pca, '~/Documents/PhD/Work/Projects/coevo_in_networks/Results/processed_output_30_30/wrangled_data/complete_data_for_network_analysis.csv')
write.csv(data_for_module_pca, '~/Documents/PhD/Work/Projects/coevo_in_networks/Results/processed_output_30_30/wrangled_data/complete_data_for_module_analysis.csv')
write.csv(data_for_species_pca, '~/Documents/PhD/Work/Projects/coevo_in_networks/Results/processed_output_30_30/wrangled_data/complete_data_for_species_analysis.csv')
write.csv(data_for_overall_pca, '~/Documents/PhD/Work/Projects/coevo_in_networks/Results/processed_output_30_30/wrangled_data/complete_data_for_overall_analysis.csv')
