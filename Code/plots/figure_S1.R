# Load required libraries
require(ggplot2)
require(dplyr)
require(cowplot)
require(ggpubr)

##############################################################
# Read data
data = read.csv('Results/processed_output/wrangled_data/complete_data_with_pca.csv',row.names = 1)

# Wrangle data 
data_small = data %>% select(network_name,module,names,modularity_SA,connectance,
                             number_of_modular_groups,
                             module_size,module_connectance,within_module_degree,
                             among_module_connectivity) 
# Build plots
plot_modularity = ggplot(data_small %>% distinct(network_name, .keep_all = TRUE), aes(x = modularity_SA)) + geom_histogram(position='identity', bins = 80) + theme_minimal() + theme(text=element_text(size = 13),aspect.ratio = 1) + xlab('Modularity') + ylab('Count')
plot_connectance = ggplot(data_small %>% distinct(network_name, .keep_all = TRUE), aes(x = connectance)) + geom_histogram(position='identity', bins = 80) + theme_minimal() + theme(text=element_text(size = 13),aspect.ratio = 1) + xlab('Connectance') + ylab('Count')
plot_number_modules = ggplot(data_small %>% distinct(network_name, .keep_all = TRUE), aes(x = number_of_modular_groups)) + geom_histogram(position='identity', bins = 100) + theme_minimal() + theme(text=element_text(size = 13),aspect.ratio = 1) + xlab('Number of modules') + ylab('Count')
plot_module_size = ggplot(data_small %>% distinct(network_name, module, .keep_all = TRUE), aes(x = module_size)) + geom_histogram(position='identity', bins = 80) + theme_minimal() + theme(text=element_text(size = 13),aspect.ratio = 1) + xlab('Module size') + ylab('Count')
plot_module_connectance =ggplot(data_small %>% distinct(network_name, module, .keep_all = TRUE), aes(x = module_connectance)) + geom_histogram(position='identity', bins = 80) + theme_minimal() + theme(text=element_text(size = 13),aspect.ratio = 1) + xlab('Module connectance') + ylab('Count')
plot_within_module_degree = ggplot(data_small %>% distinct(network_name, module, names, .keep_all = TRUE), aes(x = within_module_degree)) + geom_histogram(position='identity', bins = 80) + theme_minimal() + theme(text=element_text(size = 13),aspect.ratio = 1) + xlab('Within-module connectance') + ylab('Count')
plot_among_module_connectivity = ggplot(data_small %>% distinct(network_name, module, names, .keep_all = TRUE), aes(x = among_module_connectivity)) + geom_histogram(position='identity',bins = 100) + theme_minimal() + theme(text=element_text(size = 13),aspect.ratio = 1) + xlab('Among-module connectivity') + ylab('Count')

# Merge plots
graphics.off()
final_plot = plot_grid(plot_modularity,plot_connectance,plot_number_modules,plot_module_size,
                       plot_module_connectance,plot_within_module_degree,NULL,
                      plot_among_module_connectivity,NULL)

# In case you wish to save plot
#ggsave('Results/figures/figure_S1.pdf')

