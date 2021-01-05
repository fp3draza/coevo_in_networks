# Load required libraries
require(ggplot2)
require(dplyr)
require(cowplot)


correlation_across_m_network = function(data, method_cor){
  # This function performs the correlations between the degree of trait matching
  # at the network level with the PCA score of networks. The function requires
  # the following parameters:
	# - data: dataframe where results are stored
	# - method_cor: method of the correlation to fit

  # Obtain the vector of mutualistic selection values  
  mut_selection = unique(data$mutualistic_selection)
  
  # Determine the type of correlation to fit, can be pearson or spearman
  if (method_cor == 'pearson') {
   # Create dataframe to store the output from the correlation
   df = matrix(nrow = length(mut_selection), ncol = 4) 
   colnames(df) = c('mutualistic_selection','r','c_low','c_upp')
  }
  
  else if(method_cor == 'spearman'){
  # Create dataframe to store the output from the correlation
  df = matrix(nrow = length(mut_selection),ncol = 2)
  colnames(df) = c('mutualistic_selection','rho')
  }
  
  # Set the counter
  count = 0
  
  # Loop throught the different levels of mutualistic selection
  for (mut in mut_selection) {
    
    # Update counter
    count = count + 1
    
    # Wrangle data
    data_transformed = data %>% select(network_pca_score,z_score_network,
                                       mutualistic_selection,network_name) %>% 
      filter(mutualistic_selection == mut) %>% distinct() %>% group_by(network_name) %>%
      mutate(mean_match = mean(z_score_network)) %>%
      distinct()
    
    # Fit Correlation
    cor = cor.test(data_transformed$network_pca_score,data_transformed$mean_match, method = method_cor)
    
    # Extract results
    if (method_cor == 'pearson') {
      rho = cor$estimate
      m = mut
      c_low = cor$conf.int[2]
      c_upp = cor$conf.int[1]
      df[count,1] = m
      df[count,2] = rho
      df[count,3] = c_low
      df[count,4] = c_upp
    }
    
    
    else if(method_cor == 'spearman'){
    rho = cor$estimate
    m = mut
    
    df[count,1] = m
    df[count,2] = rho
    }
    
   # Remove subset data for current mutualistic selection
    rm(data_transformed)
    
  }
  # Make a dataframe
  df = data.frame(df)

  # Return the output of the correlation as a dataframe
  return(df)
}

correlation_across_m_module = function(data, method_cor){

  # This function calculates the correaltion between module structure and the 
  # degree of trait matching at the module level. The function requires:
 	# - data : data where trait matching and module properites are stored
	# - method_cor : the type of correlation to fit

  # Obtain the levels of mutualistic selection  
  mut_selection = unique(data$mutualistic_selection)
  
  # Determine the type of correaltion to fit, can be pearson or spearman
  if (method_cor == 'pearson') {
    # Create dataframe to store output from correlation
    df = matrix(nrow = length(mut_selection), ncol = 4) 
    colnames(df) = c('mutualistic_selection','r','c_low','c_upp')
  }
  
  else if(method_cor == 'spearman'){
    # Create dataframe to store output from correlation
    df = matrix(nrow = length(mut_selection),ncol = 2)
    colnames(df) = c('mutualistic_selection','rho')
  }
  
  # Start counter
  count = 0
  
  # Loop throught the levels of mutualistic selection
  for (mut in mut_selection) {
    
    # Update counter
    count = count + 1
    
    # Wrangle data
    data_transformed = data %>% select(module_pca_score,z_score_modules,
                                       mutualistic_selection,network_name,module) %>% 
      filter(mutualistic_selection == mut) %>%
      distinct() %>% group_by(network_name,module) %>%
      mutate(mean_match = mean(z_score_modules)) %>%
      distinct()
    
    # Fit Correlation
    cor = cor.test(data_transformed$module_pca_score,data_transformed$mean_match, method = method_cor)
    
    # Extract results
    if (method_cor == 'pearson') {
      rho = cor$estimate
      m = mut
      c_low = cor$conf.int[2]
      c_upp = cor$conf.int[1]
      df[count,1] = m
      df[count,2] = rho
      df[count,3] = c_low
      df[count,4] = c_upp
    }
    
    
    else if(method_cor == 'spearman'){
      rho = cor$estimate
      m = mut
      
      df[count,1] = m
      df[count,2] = rho
    }
    
    # Remove subset data of current mutualistic selection
    rm(data_transformed)
    
  }

  # Create dataframe
  df = data.frame(df)
  # Return results of correaltion as a dataframe
  return(df)
}

correlation_across_m_species = function(data, method_cor){

  # This function computes the correlation between species role and 
  # species trait matching. The function requires:
	# - data: data where trait matching and species role information are stored
	# - method_cor: type of correlation to fit
 
  # Obtain the levels of mutualistic selection 
  mut_selection = unique(data$mutualistic_selection)
  
  # Determine the type of correlation to fit, can be pearson or spearman
  if (method_cor == 'pearson') {
    df = matrix(nrow = length(mut_selection), ncol = 4) 
    colnames(df) = c('mutualistic_selection','r','c_low','c_upp')
  }
  
  else if(method_cor == 'spearman'){
    df = matrix(nrow = length(mut_selection),ncol = 2)
    colnames(df) = c('mutualistic_selection','rho')
  }
  
  # Start counter
  count = 0
  
  # Loop throught the levels of mutualsitic selection
  for (mut in mut_selection) {
    
    # Update counter
    count = count + 1
    
    # Wrangle data
    data_transformed = data %>% select(species_pca_score,species_matching,
                                       mutualistic_selection,network_name,module,names) %>% 
      filter(mutualistic_selection == mut) %>%
      distinct() %>% group_by(network_name,module,names) %>%
      mutate(mean_match = mean(species_matching)) %>%
      select(-species_matching) %>% distinct()
    
    # Fit Correlation
    cor = cor.test(data_transformed$species_pca_score,data_transformed$mean_match, method = method_cor)
   
    # Extract results
    if (method_cor == 'pearson') {
      rho = cor$estimate
      m = mut
      c_low = cor$conf.int[2]
      c_upp = cor$conf.int[1]
      df[count,1] = m
      df[count,2] = rho
      df[count,3] = c_low
      df[count,4] = c_upp
    }
    
    
    else if(method_cor == 'spearman'){
      rho = cor$estimate
      m = mut
      
      df[count,1] = m
      df[count,2] = rho
    }

   # Remove subset dataframe with current mutualsitic selection level 
    rm(data_transformed)    
  }
  
  # Create dataframe
  # colnames(df) = c('mutualistic_seelction','rho')
  df = data.frame(df)
  # Return correlation output as dataframe
  return(df)
}


# Load data
data = read.csv('Results/processed_output/wrangled_data/complete_data_with_pca.csv',row.names = 1)

# Fit correlation at network scale
network_correlation = correlation_across_m_network(data, 'spearman')

# Plot correlation result at network scale
network_correlation_plot = ggplot(data = network_correlation, aes(x = mutualistic_selection, y = rho)) + 
  geom_point(size = 4) + theme_classic() + labs(x = 'Mutualistic selection',y = '\u03c1') +
  theme(aspect.ratio = 1) + theme(text = element_text(size = 11)) + ggtitle('Network structure and trait matching') +
  geom_hline(yintercept = 0, linetype = 'dashed') + ylim(-1,1)

# Fit correlation at module scale 
module_correlation = correlation_across_m_module(data, 'spearman')

# Plot correlation result at module scale
module_correlation_plot = ggplot(data = module_correlation, aes(x = mutualistic_selection, y = rho)) + 
  geom_point(size = 4) + theme_classic() + labs(x = 'Mutualistic selection',y = '\u03c1') +
  theme(aspect.ratio = 1) + theme(text = element_text(size = 11)) + ggtitle('Module structure and trait matching') +
  geom_hline(yintercept = 0, linetype = 'dashed') + ylim(-1,1)

# Fit correlation at species scale
species_correaltion = correlation_across_m_species(data, 'spearman')

# Plot correlation result at species scale
species_correlation_plot = ggplot(data = species_correaltion, aes(x = mutualistic_selection, y = rho)) + 
  geom_point(size = 4) + theme_classic() + labs(x = 'Mutualistic selection', y = '\u03c1') +
  theme(aspect.ratio = 1) + theme(text = element_text(size = 11)) + ggtitle('Species role and trait matching') +
  geom_hline(yintercept = 0, linetype = 'dashed') + ylim(-1,1)

# Join plots
graphics.off()
plot_grid(network_correlation_plot,module_correlation_plot,species_correlation_plot, ncol = 3)

# In case you wish to save the plot1
#ggsave("Results/figures/figure_2.pdf")
