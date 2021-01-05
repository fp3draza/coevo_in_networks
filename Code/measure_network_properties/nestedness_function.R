# Load required packages
require(bipartite)
require(rnetcarto)

nestedness_s <- function(x){ 

  # This function computes the nestedness of a network. It follows Bastolla et al. 2009 implementation.
  # The function requires the following parameters:
	# - x : incidence matrix of the network

  # Read incidence matrix
  x = read.csv(x,header = FALSE)

  # Compute row nestedness
  nested_rows <- 0
  for(i in 1:nrow(x)){
    for(j in i:nrow(x)){
      if(j>i){
	# Sum of common interactions
        shared <- sum(x[i,]*x[j,]) 
	# Min of the degrees
        min_shared <- min(sum(x[i,]),sum(x[j,])) 
        nested_rows <- nested_rows+(shared/min_shared)
      }
    }
  }
  nestedness_rows <- nested_rows/(nrow(x)*(nrow(x)-1)/2)

  # Compute column nestedness
  nested_columns <- 0
  for(i in 1:ncol(x)){
    for(j in i:ncol(x)){
      if(j>i){
	# Sum of common interactions
        shared <- sum(x[,i]*x[,j]) 
	# Min of the degrees
        min_shared <- min(sum(x[,i]),sum(x[,j])) 
        nested_columns <- nested_columns+(shared/min_shared)
      }
    }
  }
  nestedness_columns <- nested_columns/(ncol(x)*(ncol(x)-1)/2)

  # Compute network nestedness
  nestedness_network <- (nested_rows+nested_columns)/((nrow(x)*(nrow(x)-1)/2)+(ncol(x)*(ncol(x)-1)/2))
  
  # Return row, column and network nestedness
  return(list(nestedness_rows, nestedness_columns, nestedness_network))
}

connectance = function(x){
  # This function computes the connectance of a given network
  # The function requires the following parameters:
	# - x : incidence matrix of the network

  # Read incidence matrix
  x = read.csv(x,header = FALSE)
  # Convert incidence matrix to graph
  x = graph_from_incidence_matrix(x)
  # Compute connectance
  connectance = edge_density(x)
  # Return connectance value
  return(connectance)
}

modularity_and_roles = function(x){
  # This function computes the modularity and species role of a network
  # The function requires the following parameters:
  	# - x : incidence matrix of the network

  # Read incidence matrix
  x = read.csv(x,header = FALSE)
  # Name rows and columns of matrix
  rownames(x) = paste("R",1:nrow(x),sep = '')
  colnames(x) = paste("C",1:ncol(x),sep = '')
  # Convert incidence matrix to graph
  x = graph_from_incidence_matrix(x)
  # Convert graph to adjacency matrix
  x = as.matrix(as_adjacency_matrix(x))
  # Obtain roles and modularity of adjacency matrix
  res = netcarto(x,bipartite = FALSE)
  # Return value
  return(res)
}
