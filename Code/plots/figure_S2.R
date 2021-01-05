#########################################
# Require
require(corrplot)
require(vegan)
require(factoextra)
require(RColorBrewer)

#########################################
# Load data
data = read.csv('Results/processed_output/wrangled_data/complete_data_with_pca.csv',row.names = 1)

#########################################
# Fit PCA
data = data %>% rename(`C` = connectance, `M` = modularity_SA, `NM` = number_of_modular_groups,
                       `MS` = module_size, `MC` = module_connectance,
                       `WMD` = within_module_degree, `AMC` = among_module_connectivity)

network_PCA = PCA(data[,c("C","M","NM")], graph = FALSE)
module_PCA = PCA(data[,c("MS","MC")], graph = FALSE)
species_PCA = PCA(data[,c("WMD","AMC")], graph = FALSE)

# Extract PCA data
net_vars = get_pca_var(network_PCA)
mod_vars = get_pca_var(module_PCA)
sp_vars = get_pca_var(species_PCA)

# Define color palette
col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Plot coordinates of PCA for networks
corrplot(net_vars$coor, is.corr=FALSE, tl.col = "black", tl.srt = 45, addCoef.col = 'black', col = col(200), cl.align="l")

# Plot coordinates of PCA for modules
corrplot(mod_vars$coor, is.corr=FALSE, tl.col = "black", tl.srt = 45, addCoef.col = 'black', col = col(200), cl.align="l")

# Plot coordinates of PCA for species
corrplot(sp_vars$coor, is.corr=FALSE, tl.col = "black", tl.srt = 45, addCoef.col = 'black', col = col(200), cl.align="l")

# These plots cant be joined using cowplot, I joined them by 'hand'
# And saved them by hand