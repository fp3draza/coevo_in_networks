# Load required libraries
require(ggplot2)
require(dplyr)
require(ggpubr)

# Load data
data = read.csv('Results/processed_output/wrangled_data/complete_data_with_pca.csv',row.names = 1)

###########################################################
# Network Scale

# Graph correlation of networks PCA vs network trait matching
data_networks = data %>% select(network_pca_score,z_score_network,
                                mutualistic_selection,network_name) %>% filter(mutualistic_selection %in% c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) %>%
  distinct %>% group_by(network_name, mutualistic_selection) %>%
  mutate(mean_match = mean(z_score_network)) %>%
  distinct()


# Build and save plot
graphics.off()
p = ggscatter(data_networks, x = 'network_pca_score', y = 'mean_match', 
              add = 'reg.line', add.params = list(color = '#FC4E07', fill = 'lightgray'), conf.int = TRUE,
              cor.coef = TRUE, cor.method = 'spearman', size = 0.85, alpha = 0.7,
              font.label = c(size = 8,'plain')) + theme_minimal() + 
  facet_wrap(~mutualistic_selection) + theme(text=element_text(size = 11),aspect.ratio = 1) + xlab('Network Structure') + 
  ylab('Trait Matching') 
# In case you want to save plot
#ggsave('Results/figures_v2/figure_S3.pdf')

###########################################################
# Module Scale

# Graph correlation of networks PCA vs network trait matching
data_module = data %>% select(module_pca_score,z_score_modules,
                              mutualistic_selection,network_name,module) %>% filter(mutualistic_selection %in% c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) %>%
  distinct %>% group_by(network_name, module, mutualistic_selection) %>%
  mutate(mean_match = mean(z_score_modules)) %>%
  distinct()


# Build and save plot
graphics.off()
p = ggscatter(data_module, x = 'module_pca_score', y = 'mean_match', 
              add = 'reg.line', add.params = list(color = '#FC4E07', fill = 'lightgray'), conf.int = TRUE,
              cor.coef = TRUE, cor.method = 'spearman', size = 0.85, alpha = 0.7,
              font.label = c(size = 8,'plain'), cor.coeff.args = list(label.x.npc = 'left', label.y.npc = 'bottom')) + theme_minimal() + 
  facet_wrap(~mutualistic_selection) + theme(text=element_text(size = 11),aspect.ratio = 1) + xlab('Module Structure') + 
  ylab('Trait Matching') 
# In case you want to save plot
#ggsave('Results/figures_v2/figure_S4.pdf')

###########################################################
# Species Scale

# Graph correlation of networks PCA vs network trait matching
data_species = data %>% select(species_pca_score,z_score_species,
                               mutualistic_selection,network_name,module,names) %>% filter(mutualistic_selection %in% c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) %>%
  distinct %>% group_by(network_name, module, names, mutualistic_selection) %>%
  mutate(mean_match = mean(z_score_species)) %>%
  distinct()


# Build and save plot
graphics.off()
p = ggscatter(data_species, x = 'species_pca_score', y = 'mean_match', 
              add = 'reg.line', add.params = list(color = '#FC4E07', fill = 'lightgray'), conf.int = TRUE,
              cor.coef = TRUE, cor.method = 'spearman', size = 0.85, alpha = 0.7,
              font.label = c(size = 8,'plain'), cor.coeff.args = list(label.x.npc = 'left', label.y.npc = 'bottom')) + theme_minimal() + 
  facet_wrap(~mutualistic_selection) + theme(text=element_text(size = 11),aspect.ratio = 1) + xlab('Species Role') + 
  ylab('Trait Matching') 
# In case you want to save plot
#ggsave('Results/figures_v2/figure_S5.pdf')


