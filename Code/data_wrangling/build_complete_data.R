# This scripts joins all processed output results into a single dataframe

# Read data

# Read trait matching results
trait_matching_data = read.csv('Results/processed_output/trait_matching/merged_results.csv',row.names = 1)
# Read species role data
species_role_data = read.csv('Results/processed_output/network_roles/network_roles_niche.csv',row.names = 1)
# Read network properties data
network_data = read.csv('Results/processed_output/network_metrics/modular_nested_metrics.csv',row.names = 1)
# Read module properties data
module_data = read.csv('Results/processed_output/module_metrics/module_metrics.csv',row.names = 1)

# Calculate standardized trait matching values

# Zscore for networks
trait_matching_data = trait_matching_data %>%
                      group_by(mutualistic_selection, network_name) %>%
                      mutate(mean_network_matching = mean(network_matching,na.rm = TRUE))

trait_matching_data = trait_matching_data %>%
                      group_by(mutualistic_selection) %>%
                      mutate(z_score_network = (mean_network_matching - mean(mean_network_matching))/sd(mean_network_matching)) %>%
                      mutate(z_score_network_level = if_else(z_score_network > 1, 'high',if_else(z_score_network < -1,'low','average')))

# Zscore for modules
trait_matching_data = trait_matching_data %>%
                      group_by(mutualistic_selection,network_name,module) %>%
                      mutate(mean_module_matching = mean(module_matching_vector,na.rm = TRUE))

trait_matching_data = trait_matching_data %>%
                      group_by(mutualistic_selection) %>%
                      mutate(z_score_modules = (mean_module_matching - mean(mean_module_matching))/sd(mean_module_matching)) %>%
                      mutate(z_score_modules_level = if_else(z_score_modules > 1, 'high',if_else(z_score_modules < -1,'low','average')))

# Zscore for species
trait_matching_data = trait_matching_data %>%
                      group_by(mutualistic_selection, network_name, module, names) %>%
                      mutate(mean_species_matching = mean(species_matching, na.rm =  TRUE))

trait_matching_data = trait_matching_data %>%
                      group_by(mutualistic_selection) %>%
                      mutate(z_score_species = (mean_species_matching - mean(mean_species_matching))/sd(mean_species_matching)) %>%
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
module_data = module_data %>% mutate(module = module + 1)
complete_data = left_join(merged_data,module_data)


# Save merged dataframe
write.csv(complete_data, 'Results/processed_output/wrangled_data/complete_data.csv')
