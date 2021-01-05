# Load required libraries
require(ggplot2)
require(vegan)
require(corrplot)
require(ggpubr)

compute_rda = function(data, m){
  
  # This function computes the Redundancy Analysis.
  # It requires the following parameters:
    # - data: data where trait matching and network descriptor data is stored
    # - m : levels of mutualistic selection
  
  # Fit RDA with the degree of trait matching achieved at the species level
  # and the network descriptors
  pca_result = rda(data_sub %>% 
                     select(modularity_SA,connectance,
                            number_of_modular_groups,
                            module_size,module_connectance,
                            within_module_degree,
                            among_module_connectivity) ~ data_sub$z_score_species,scale = TRUE)
  
  # Extract the eigenvalue of the RDA
  eig = pca_result$CCA$eig
  
  # Return the value
  return(eig)
}

# Load data
data = read.csv('Results/processed_output/wrangled_data/complete_data_with_pca.csv',row.names = 1)
# Extract the levels of mutualistic selection
m_seq = sort(unique(data$mutualistic_selection))
# Initialize an empty vector where the eigenvalue of RDAs will be stored
eig_vector = rep(0,length(m_seq))
# Start counter
count = 1
# Loop through the mutualistic selection values
for (m_val in m_seq) {
  # Subset data by the current mutualistic selection value
  data_sub = data %>% filter(mutualistic_selection == as.numeric(m_val))
  # Extract the eigenvalue of RDA and store in vector
  eig_vector[count] = compute_rda(data_sub,m)
  # Update counter
  count = count + 1
}

# Create a tibble from the mutualistic selection values and eigenvalue
data_plot = tibble(m = m_seq, eig = eig_vector)

# Create Plot
ggplot(data=data_plot, aes(x = m, y = eig)) + geom_point(size = 4) + theme_classic() + labs(x = 'Mutualistic selection',y = 'Role of network descriptor') +
  theme(aspect.ratio = 1) + theme(text = element_text(size = 18)) 


compute_rda_relation = function(data, m){
  # This function extracts the two largest coordinates of the RDA
  # at each level of mutualistic selection. It requires:
    # - data: data where trait matching and network descriptor information is stored
 
  # Load required libraries
  require(corrplot)

  # Fit RDA 
  pca_result = rda(data_sub %>% 
                     select(modularity_SA,connectance,
                            number_of_modular_groups,
                            module_size,module_connectance,
                            within_module_degree,
                            among_module_connectivity) ~ data_sub$z_score_species,scale = TRUE)
  
  # Extract the scores of the network descriptors
  scores = scores(pca_result)$species
  
  # Extract list of names of network descriptors
  names = rownames(scores)
  
  # Extract the index, names and values of the largest coordinate
  max_index = which.max(scores[,1])
  max_name = names[max_index]
  max_val = scores[max_index,1]
  
  # Extract the index, names and values of the smallest coordinate
  min_index = which.min(scores[,1])
  min_name = names[min_index]
  min_val = scores[min_index,1]
  
  # Extract the coordiantes of the Trait Matching Variables
  tm_name = 'trait matching'
  tm_val = pca_result$CCA$biplot
  
  # If the coordinate of Trait Matching is negative, convert it to positive
  # This is done so in the plot, positive values always indicate an increase
  # of trait matching
  if (tm_val < 0){
    tm_val = tm_val * -1
    max_val = max_val * -1
    min_val = min_val * -1
  }
  
  # Create a vector with the max value
  max = c(max_name, max_val,m)
  # Create a vector with the min value
  min = c(min_name, min_val,m)
  # Bind min and max
  full = rbind(max,min)
  # Convert to dataframe
  full = as.data.frame(full)
  
  # Return the coordinates as a dataframe
  return(full)
}


# Load data
data = read.csv('Results/processed_output/wrangled_data/complete_data_with_pca.csv',row.names = 1)
# Obtain the vector of mutualistic selection values
m_seq = sort(unique(data$mutualistic_selection))
# Create empty dataframe
df = NULL
# Start count
count = 0 
# Loop through the mutualistic selection values
for (m_val in m_seq) {
  # Subset data by the current mutualistic selection values
  data_sub = data %>% filter(mutualistic_selection == as.numeric(m_val))
  # Fit RDA and extract the min and max coordinates
  res = compute_rda_relation(data_sub,m_val)
  # Update the dataframe
  df = rbind(df,res)
  # Update the counter
  count = count + 1
}

# Name the columns of the dataframe
colnames(df) = c('var','val','m')
# Convert to numeric type
df$val = as.numeric(levels(df$val))
# Fill level of mutualistic selection
df$m = rep(as.numeric(levels(df$m)),each = 2)

# Define color of coordinate depending on if it increases or decreases trait matching
df = df %>% mutate(co = if_else(val < 0, 'darkblue', 'darkred'))
# Change names of coordinates for plot
df = df %>% mutate(lab = if_else(var == 'connectance','C',if_else(var == 'modularity_SA', 'M', if_else(var == 'number_of_modular_groups', 'NM', if_else(var == 'within_module_degree', 'WMD','NA')))))
# Generate plot
ggdotchart(df, x = "m", y = "val",   
           color = "co",
           palette = c("#00AFBB","#FC4E07"),
           add = "segments", 
           sorting = 'none',
           dot.size = 4,  
           ggtheme = theme_pubr()
) +  geom_hline(yintercept=0, linetype="dashed", color = "black") +  theme_classic() +
  theme(aspect.ratio = 1) + theme(text = element_text(size = 14), legend.position = 'none')  + labs(x = 'Mutualistic selection',y = 'Coordinates of network descriptors')
ggsave('Results/figures/figure_3.pdf')