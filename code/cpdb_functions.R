#some functions
make_cpdb_df <- function(cpdb_dir,senders, receivers, cluster_interactions = TRUE, 
                         filter_secreted  = FALSE, remove_auto = TRUE, filter_integrins = TRUE,
                         scaling = c("min_max", "scale"),
                         annotations_drop = c("guidetopharmacology.org")){
  message("reading in files")
  cpdb_files <- pblapply(list.files(cpdb_dir), function(f){
    read.table(file.path(cpdb_dir, f), sep = "\t", quote = NULL, header = TRUE)
  })
  names(cpdb_files) <- gsub(".txt", "", list.files(cpdb_dir))
  cpdb_index <- data.frame(cpdb_files$means[, 1:11],row.names = cpdb_files$means[, 1])
  means <- reshape::melt(data.frame(cpdb_files$means[, 12:ncol(cpdb_files$means)], row.names = cpdb_index$id_cp_interaction, "interaction" = cpdb_index$id_cp_interaction))
  p_vals <- reshape::melt(data.frame(cpdb_files$pvalues[, 12:ncol(cpdb_files$pvalues)], row.names = cpdb_index$id_cp_interaction, "interaction" = cpdb_index$id_cp_interaction))
  
  cell_a <- unlist(lapply(strsplit(as.character(means$variable), "[.]"), function(x){x[1]}))
  cell_b <- unlist(lapply(strsplit(as.character(means$variable), "[.]"), function(x){x[2]}))
  interactions <- cpdb_index[means$interaction, ]
  molecule_a <- unlist(lapply(strsplit(as.character(interactions$interacting_pair), "[_]"), function(x){x[1]}))
  molecule_b <- unlist(lapply(strsplit(as.character(interactions$interacting_pair), "[_]"), function(x){x[2]}))
  cpdb_df <- data.frame(interaction = means$interaction, cells = means$variable, 
                        mean = means$value, pval = p_vals$value, cell_a = cell_a, cell_b = cell_b, cpdb_index[means$interaction, ],
                        molecule_a = molecule_a, molecule_b = molecule_b)
  #now we have a bound dataframe
  
  #now we want to get the ligands and receptors the right way round because inexplicably cpdb get's this all askew
  cpdb_df$receptor_cell <- NA
  cpdb_df$receptor_cell[cpdb_df$receptor_a %in% "True"] <- cpdb_df$cell_a[cpdb_df$receptor_a %in% "True"]
  cpdb_df$receptor_cell[cpdb_df$receptor_b %in% "True"] <- cpdb_df$cell_b[cpdb_df$receptor_b %in% "True"]
  cpdb_df$ligand_cell <- NA
  cpdb_df$ligand_cell[cpdb_df$receptor_a %in% "False"] <- cpdb_df$cell_a[cpdb_df$receptor_a %in% "False"]
  cpdb_df$ligand_cell[cpdb_df$receptor_b %in% "False"] <- cpdb_df$cell_b[cpdb_df$receptor_b %in% "False"]
  cpdb_df$receptor <- NA
  cpdb_df$receptor[cpdb_df$receptor_a %in% "True"] <- cpdb_df$molecule_a[cpdb_df$receptor_a %in% "True"]
  cpdb_df$receptor[cpdb_df$receptor_b %in% "True"] <- cpdb_df$molecule_b[cpdb_df$receptor_b %in% "True"]
  cpdb_df$ligand <- NA
  cpdb_df$ligand[cpdb_df$receptor_a %in% "False"] <- cpdb_df$molecule_a[cpdb_df$receptor_a %in% "False"]
  cpdb_df$ligand[cpdb_df$receptor_b %in% "False"] <- cpdb_df$molecule_b[cpdb_df$receptor_b %in% "False"]
  
  cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "ligand_cell"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "cell_a"]
  cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "receptor_cell"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "cell_b"]
  cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "ligand"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "molecule_a"]
  cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "receptor"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "molecule_b"]
  
  cpdb_df$LR <- paste0(cpdb_df$ligand, "->", cpdb_df$receptor)
  cpdb_df$SR <- paste0(cpdb_df$ligand_cell, "->", cpdb_df$receptor_cell)
  
  #get these ordered by how we provided them
  egrid <- expand.grid(receivers, senders)
  cpdb_df$SR <- factor(cpdb_df$SR, levels =  paste0(egrid$Var2, "->", egrid$Var1))
  
  #subset to just secreted interactions
  if(filter_secreted){
    cpdb_df <- cpdb_df[cpdb_df$secreted %in% "True", ]
  }
  
  #remoove annotation strategy
  for(i in annotations_drop){
    cpdb_df <- cpdb_df[!cpdb_df$annotation_strategy %in% i, ]
  }
  
  #remove auto interactions
  if(remove_auto){
    cpdb_df <- cpdb_df[!cpdb_df$cell_a == cpdb_df$cell_b, ]
    cpdb_df <- cpdb_df[!cpdb_df$molecule_a == cpdb_df$molecule_b, ]
  }
  
  if(filter_integrins){
    cpdb_df <- cpdb_df[cpdb_df$is_integrin %in% "False", ]
  }
  
  #subset to the cells we are interested in 
  cell_types <- union(senders, receivers)
  cpdb_df <- cpdb_df[cpdb_df$cell_a %in% cell_types + cpdb_df$cell_b %in% cell_types == 2, ]
  
  #subset to senders and receivers
  cpdb_df <- cpdb_df[cpdb_df$ligand_cell %in% senders, ]
  cpdb_df <- cpdb_df[cpdb_df$receptor_cell %in% receivers, ]
  
  #remove totally nonsignificant interactions
  #now find the interacting pairs which are all non significant 
  message("filtering non significant interactions")
  nonsig_interactions <- unique(cpdb_df$interaction)[unlist(lapply(unique(cpdb_df$interaction), function(x){
    return(sum( cpdb_df[cpdb_df$interaction %in% x, "pval"]) == sum(cpdb_df$interaction %in% x))
  }))]
  
  cpdb_df <- cpdb_df[!cpdb_df$interaction %in% nonsig_interactions, ]
  #and render nonsignificant p values null
  cpdb_df$pval[cpdb_df$pval == 1] =NA
  
  #make sender receiver a factor
  cpdb_df$SR <- factor(cpdb_df$SR)
  
  #scale the data
  mat <- matrix(0, ncol = length(unique(cpdb_df$LR)), nrow = length(unique(cpdb_df$SR)))
  colnames(mat) <- unique(cpdb_df$LR)
  rownames(mat) <- unique(cpdb_df$SR)
  message("clustering and/or scaling step...")
  pb = txtProgressBar(min = 1, max = length(unique(cpdb_df$LR)), initial = 1) 
  for(i in unique(cpdb_df$LR)){
    setTxtProgressBar(pb, grep(i, unique(cpdb_df$LR)))
    
    for(j in unique(cpdb_df$SR)){
      mat[j, i] <- log1p(cpdb_df[cpdb_df$SR %in% j & cpdb_df$LR %in% i, "mean"])
    }
  }
  
  if(scaling == "min_max"){
    mat <- apply(mat, 2, min_max_normalisation)
  }
  if(scaling == "scale"){
    mat <- scale(mat)
  }
  melted_mat <- reshape::melt(mat)
  rownames(melted_mat) <- paste0(melted_mat$X1, melted_mat$X2)
  cpdb_df$scaled <- melted_mat[paste0(cpdb_df$SR, cpdb_df$LR), "value"]
  
  
  
  #now cluster the interactions
  if(cluster_interactions){
    message("clustering interactions")
    cpdb_df$LR <- factor(cpdb_df$LR, levels = colnames(mat)[hclust(dist(t(mat)), method = "ward.D")$order])
  }else{
    cpdb_df$LR <- factor(cpdb_df$LR)
  }
  return(cpdb_df)
  
}
cellphone_db_plot <- function(cpdb_dir,senders, receivers, cluster_interactions = TRUE, 
                              filter_secreted  = FALSE, scaling = c("min_max", "scale"), remove_auto = TRUE,
                              filter_integrins = TRUE, 
                              annotations_drop = c("guidetopharmacology.org")){
  cpdb_df <- make_cpdb_df(cpdb_dir = cpdb_dir, senders = senders, receivers = receivers, cluster_interactions = cluster_interactions,
                          filter_secreted = filter_secreted, scaling = scaling, remove_auto=remove_auto,
                          filter_integrins = filter_integrins, 
                          annotations_drop = annotations_drop)
  #plot
  message("plotting result")
  if(scaling == "scale"){
    scaled_colors <- c("palegreen4", "grey80", "plum3")
    max_scaled <- max(abs(cpdb_df$scaled))
    pl <- ggplot(cpdb_df, aes(x = SR, y = LR, size = -log10(pval +1e-3 ), fill = scaled )) + geom_point(pch = 21) + theme_bw() +
      scale_size_continuous(limits = c(0, 4)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank()) + scale_fill_gradient2(low = scaled_colors[1], mid = scaled_colors[2], high = scaled_colors[3],
                                                                 limits = c(-max_scaled, max_scaled ))
    
  }
  if(scaling == "min_max"){
    pl <- ggplot(cpdb_df, aes(x = SR, y = LR, size = -log10(pval +1e-3 ), fill = scaled)) + geom_point(pch = 21) + theme_bw() +
      scale_size_continuous(limits = c(0, 4)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank()) + scale_fill_gradientn(colors = viridis::magma(100), limits = c(0, 1)) 
  }
  if(scaling == FALSE){
    pl <- ggplot(cpdb_df, aes(x = SR, y = LR, size = -log10(pval +1e-3 ), fill = log1p(mean))) + geom_point(pch = 21) + theme_bw() +
      scale_size_continuous(limits = c(0, 4)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank()) + scale_fill_gradientn(colors = viridis::magma(100), limits = c(0, max(log1p(cpdb_df$mean)))) 
    
    
  }
  message("done")
  return(pl)
  
}
