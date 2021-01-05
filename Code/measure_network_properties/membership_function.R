membership_vector = function(membership_file,color_pal){
  # This function builds a vector of species organized by module membership   
  # The function requires:
	# membership_file: output file from MODULAR software
	# color_pal : color palette

  # Load required packages
  require(RColorBrewer)
  members = membership_file

  # Build vector of membership for column species
  members_cols = members[which(grepl('C', members$V1, fixed = TRUE)),]
  members_cols = members_cols[order(members_cols$V2),]
  colors_by_module_cols = brewer.pal(nlevels(as.factor(members_cols$V2)),name=color_pal)
  members_cols$color = rep(NA,dim(members_cols)[1])

  # Build vector of membership for row species
  members_rows = members[which(grepl('R', members$V1, fixed = TRUE)),]
  members_rows = members_rows[order(members_rows$V2),]
  colors_by_module_rows = brewer.pal(nlevels(as.factor(members_rows$V2)),name=color_pal)
  members_rows$color = rep(NA,dim(members_rows)[1])
  
  for (i in 1:dim(members_rows)[1]) {
    members_cols$color[i] = colors_by_module_cols[members_cols$V2[i]]
    members_rows$color[i] = colors_by_module_rows[members_rows$V2[i]]
  }

  return(list(members_cols,members_rows))
}
