#import functions #####
source("code/functions.R")
ad <- import("anndata")
sc <- import("scanpy")
#read in data ####
adata <- sc$read_h5ad("data/scvi_trained_adata.h5ad")
#calculate stress scores ####
vdb_stress_genes <- c("JUNB", "FOSB", "ZFP36", "EGR1", "FOS",
                      "HSPA1A", "HSPA1B", "HSPA8", "HSPE1", "JUND",
                      "CEBPB", "ATF3", "SOCS3", "HSP90AB1", "HSPB1") #these are genes from van den brink
vdb_stress_genes <- adata$var_names$values[adata$var$Symbol %in% vdb_stress_genes]
sc$tl$score_genes(adata, vdb_stress_genes, score_name = "stress_genes")

#doublet and QC filtering ####
#do neighbors and cluster
sc$pp$neighbors(adata, use_rep = "X_scVI", n_neighbors = as.integer(20))
sc$tl$leiden(adata, resolution = 1.5, key_added = "qc_filtering_clusters")
#qc markers
gs <- import("geosketch")
qc_gs_idx <- gs$gs(X = adata$obsm$get("X_scVI"), N = as.integer(50000))
qc_gs_idx <- 1:nrow(adata$obs) %in% unlist(qc_gs_idx)
qc_markers <- tfidf_all_markers(adata=adata[qc_gs_idx], groups=  "qc_filtering_clusters",
                             gene_subset = adata$var_names$values[adata$var$highly_variable]) #do markers 
#find clusters enriched for doublets
tab = table(adata$obs$qc_filtering_clusters, adata$obs$scrublet_cuts)
tab = tab/rowSums(tab)

#find clusters marked out by poor QC
qc_metrics <- c("log1p_n_genes_by_counts", "pct_counts_MT", "log1p_total_counts","log1p_total_counts_MT", "pct_counts_in_top_500_genes", "pct_counts_in_top_50_genes", "stress_genes")
qc_df <- scale(data.frame(aggregate(adata$obs[, qc_metrics], by =list(adata$obs$qc_filtering_clusters), mean), "fraction_doublets" = tab[, 2],  row.names = 1))
qc_df_long <- reshape::melt(qc_df)
qc_df_long$X1 <- factor(qc_df_long$X1)
ggplot(qc_df_long, aes(x= X1, y=X2, fill = value)) + geom_tile() + coord_fixed() + 
  scale_fill_gradient2(low=  "darkblue", mid = "grey95", high = "red", midpoint = 0,
                       limit = c(-6,6)) + theme_bw() + theme(axis.text = element_text(color = "black"), axis.title = element_blank()) 

cluster_qc_remove <- c("16", "23") #these are the clusters we need to remove because they look low quality.
adata <- adata[!adata$obs$qc_filtering_clusters %in% cluster_qc_remove] #these are T-B doublets

#neighbors cluster umap  #####
sc$pp$neighbors(adata, use_rep = "X_scVI", n_neighbors = as.integer(20))
sc$tl$umap(adata, min_dist = 0.3)
sc$tl$leiden(adata, resolution = 1)
annotated_raster_plot(adata = adata, group= "leiden") #plot by cluster
leiden_cl <- adata$obs$leiden
#marker finding per cluster ####
markers <- tfidf_all_markers(adata=adata,
                             gene_subset = adata$var_names$values[adata$var$highly_variable]) #do markers 

#subsplit clusters if they need it ####
#split apart some clusters which look like they have more unappreciated heterogeneity.. 
adata <- recluster_leiden(adata, clusters = "14", resolution = 0.2) 
adata <- recluster_leiden(adata, clusters = "13", resolution = 0.2) 
adata <- recluster_leiden(adata, clusters = "16", resolution = 0.2) 
adata <- recluster_leiden(adata, clusters = "15", resolution = 0.2) 
annotated_raster_plot(adata = adata, group= "leiden") #plot by cluster

#get some more markers on the subsetted data #####
markers <- tfidf_all_markers(adata=adata, gene_subset = adata$var_names$values[adata$var$highly_variable]) #do markers 

#save this clustered object  - this is an intermediate stopping point :)  #####
adata$write_h5ad("data/clustered_HLN_adata.h5ad")
#read back the clustered object #####
adata <- sc$read_h5ad("data/clustered_HLN_adata.h5ad")

#do celltypist predictions ####
celltypist <- import("celltypist")
celltypist_in <- adata$copy()
celltypist_in$X <- celltypist_in$layers$get("counts")
sc$pp$normalize_total(celltypist_in, target_sum = 1e4) #this is how celltypist like it 
sc$pp$log1p(celltypist_in)
celltypist_in$var_names <- celltypist_in$var$Symbol
celltypist_in$var_names_make_unique()
#compute an overclustering
sc$tl$leiden(adata, resolution = 2.5, key_added = "celltypist_overclusters")
#write the celltypist object out.
celltypist_in$write_h5ad("data/celltypist_in.h5ad")
write.table(celltypist_in$obs$celltypist_overclusters, "data/celltypist_overclustering.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
ct_predictions = celltypist$annotate("data/celltypist_in.h5ad", model = 'Immune_All_Low.pkl')
saveRDS(ct_predictions$predicted_labels, "data/celltypist_predictions.RDS")
saveRDS(ct_predictions$probability_table, "data/celltypist_probabilities.RDS")

#take a quick look at celltypist outputs.. these are sensible.. 
celltypist_colors <- colorRampPalette(c("grey90", "grey50", "black"))(10)
ggplot(data.frame(x=celltypist_in$obsm$get("X_umap")[, 1],
                  y=  celltypist_in$obsm$get("X_umap")[, 2], prediction = ct_predictions$predicted_labels$'predicted labels'),
       aes(x=x,y=y, color = prediction) ) + ggrastr::geom_point_rast(pch = 19, cex=0.01)  +  theme_classic() + 
  guides(colour = guide_legend(override.aes = list(size=5)))

#plot a table of celltypist outputs
ct_tab <- table(celltypist_in$obs$leiden, ct_predictions$predicted_labels$'predicted labels')
ct_df <- reshape::melt(ct_tab/rowSums(ct_tab))
sc$tl$dendrogram(celltypist_in, groupby = "leiden", use_rep = "X_scVI")
ct_df$Var.2 <- factor(ct_df$Var.2, levels = colnames(ct_tab)[hclust(dist(t(ct_tab)), method = "complete")$order])
ct_df$Var.1 <- factor(ct_df$Var.1,  levels =celltypist_in$uns$get("dendrogram_leiden")$categories_ordered)
ggplot(ct_df, aes(x= Var.1, y=  Var.2, fill = value)) + ggrastr::geom_tile_rast() + 
  scale_fill_gradientn(colours =celltypist_colors) + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))

#plot on single UMAP maps
ag_pred <- data.frame(aggregate(ct_predictions$probability_table, by = list(celltypist_in$obs$leiden), mean), row.names = 1)
cell_types_plot <- c("aDC", "pDC", "DC1", "DC2", "DC3", "Macrophages", "ILC", "Naive B cells", "Memory B cells", "Regulatory T cells")
ag_plot_df <- do.call(rbind, lapply(cell_types_plot, function(x){
  message(x)
  p_table <- ct_predictions$probability_table
  df  <- data.frame(x=celltypist_in$obsm$get("X_umap")[, 1],
                    y=  celltypist_in$obsm$get("X_umap")[, 2], prediction = p_table[, colnames(p_table) %in% x], 
                    celltype = x)
  return(df)
}))
ggplot(ag_plot_df, aes(x=x,y=y, color = prediction) ) + ggrastr::geom_point_rast(pch = 19, cex=0.01)  +  theme_classic() + 
  scale_color_gradientn(colours = celltypist_colors) + facet_wrap(~celltype)

#annotate ####
#make some annotatons
annotations <- c(#b cells
                      "7" = "naive_b", #naive B D+M
                      "2" = "naive_b", #naive B D+M
                      "1" = "memory_b", #mem markers
                      "17" = "germinal_centre_b", #GC AICDA SH/clswitch
                      '14,3' = "proliferating_b", #prl
                      "19" = "plasmablast", #pb
                    #t cells
                 #innate lymphocytes
                 "13,2" = "ILC", #characteristic KLRB sig' innate
                 "13,3" = "NK_1",  #essentially  brights and dims
                 "13,1" = "NK_2", #essentially  brights and dims
                 "13,0" = "NK_2", #essentially  brights and dims
                "5" = "CD8_Tem", #Tem
                "14,2" = "CD8_Tem_proliferating", #prl + Tem
                 "14,1" = "CD8_Texh", #exhaustion markers
                 "12" = "CD8_Tem_terminal",
                 "14,0" = "CD4_Texh_proliferating", #prol and exh markers
                 "8" = "CD4_Texh", #numerous checkpoints
                 "3" = "CD4_Tfh", #CXCL13 & CXCR5 expression etc.
                 "10" = "CD4_Treg", #FOXP3 and CTLA4 amongst others
                 "4" = "CD4_Th17", #prominent KLRB1 & RORA expression
                 "9" = "CD4_Th1", #looks like inflammatory NFkB actiavtion signature, NR4A2, TNFAIP3, SOCS3
                 "11" = "CD4_Th0", #naive
                 "0" = "CD4_Th0",
                  "6" = "CD8_Tcm", #central mem
                 #myeloid
                 "20" = "mast", #characteristic
                 "15,0" = "classical_monocyte",
                 "15,1" = "activated_DC", #FSCN, CCR7, LAMP3 etc.
                 "15,2" = "cDC2", #CLEC10A, CD1c
                 "15,3" = "cDC1", #CLEC9A, IDO etc.
                 "15,4" = "macrophage", #mphage M2 and tissue mac C1Q, RNase, FOLR2
                 "15,5" = "macrophage",
                 "18" = "plasmacytoid_DC", #pDC signature
                 #stroma
                  "16,2" = "fibroblast", #Very characteristic, LUM DCN etc..
                 "16,1" = "lymphatic_endothelium", #NTS, PROX1 ACKR4,
                 "16,0" = "fenestrated_endothelium", #PLVAP
                 "16,3" = "nonfenestrated_endothelium" #non PLVAP
)
#check we have everything
unique(adata$obs$leiden)[!unique(adata$obs$leiden) %in% names(annotations)]
#now annotate the cells
adata$obs$cell_type <-factor( annotations[match(adata$obs$leiden, names(annotations))])
#plot this
annotated_raster_plot(adata = adata[gs_idx], group= "cell_type") #plot by cluster


#plot the celltypist labels and also aoki labels according to celltype #####
celltypist <- import("celltypist")
celltypist_in <- sc$read_h5ad("data/celltypist_in.h5ad")
celltypist_in$obs$cell_type <-factor( annotations[match(celltypist_in$obs$leiden, names(annotations))])
#replot this as jaccard distance
ct_predictions <- readRDS("data/celltypist_predictions.RDS")
celltypist_in$obs$predicted_label <- ct_predictions$'predicted labels'
celltypist_jmat <- do.call(cbind, pblapply(unique(celltypist_in$obs$predicted_label), function(i){
  message(i)
  jvec <- list()
  for(j in unique(celltypist_in$obs$cell_type)){
    i_vec <- celltypist_in$obs_names$values[celltypist_in$obs$predicted_label %in% i]
    j_vec <- celltypist_in$obs_names$values[celltypist_in$obs$cell_type %in% j]
    jvec[j] <- length(intersect(i_vec, j_vec))/length(union(i_vec, j_vec))
  }
  jvec <- unlist(jvec)
  return(jvec)
}))
colnames(celltypist_jmat) <- unique(celltypist_in$obs$predicted_label)
celltypist_jmat_df <- reshape::melt(celltypist_jmat)
celltypist_jmat_df$prediction_set <- "celltypist"
colnames(celltypist_jmat_df) <- c("manual_annotations", "external_annotations", "value", "prediction_set")
ggplot(celltypist_jmat_df, aes(x= celltypist_annotations, y=  manual_annotations, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colours = celltypist_colors) + coord_fixed() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) + ylab("celltypist_labels") + xlab("manual_annotation")

#read in the stiedl data
load("data/merged_sc_obj_CancerDiscovery.Rdata")
sce <- sc.obj
#match these labels and the stiedl labels.
sce <- sce[, adata$obs_names$values[adata$obs$centre %in% "bc_cancer"]]
adata_stiedl <- adata[adata$obs$centre %in% "bc_cancer"]
adata_stiedl$obs$cell_type <- as.character(adata_stiedl$obs$cell_type)
adata_stiedl$obs$stield_labels[match(colnames(sce), adata_stiedl$obs_names$values)] <- as.character(sce$cluster_name)
stiedl_jmat <- do.call(cbind, pblapply(unique(adata_stiedl$obs$stield_labels), function(i){
  message(i)
  jvec <- list()
  for(j in unique(adata_stiedl$obs$cell_type)){
    i_vec <- adata_stiedl$obs_names$values[adata_stiedl$obs$stield_labels %in% i]
    j_vec <- adata_stiedl$obs_names$values[adata_stiedl$obs$cell_type %in% j]
    jvec[j] <- length(intersect(i_vec, j_vec))/length(union(i_vec, j_vec))
  }
  jvec <- unlist(jvec)
  return(jvec)
}))
colnames(stiedl_jmat) <- unique(adata_stiedl$obs$stield_labels)
stiedl_jmat_df <- reshape::melt(stiedl_jmat)
stiedl_jmat_df$prediction_set <- "stiedl"
colnames(stiedl_jmat_df) <- c("manual_annotations", "external_annotations", "value", "prediction_set")

#bind these together
jmat_df <- rbind(celltypist_jmat_df, stiedl_jmat_df)
jmat_df$manual_annotations <- factor(jmat_df$manual_annotations, levels = unique(annotations)) #order these by our annotation order.
ggplot(jmat_df, aes(x=external_annotations , y=  manual_annotations, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colours = colorRampPalette(colors = c("grey90", "black"))(100), limits = c(0,1)) + coord_fixed() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) + facet_wrap(~prediction_set, nrow = 1) + 
  xlab("external annotation") + ylab("manual annotation")

#make a colored embedding ####
annotated_raster_plot(adata = adata, group= "cell_type") #plot by celltype
#make celltype a factor
adata$obs$cell_type <- factor(adata$obs$cell_type, 
                              levels = unique(annotations))
#make some celltype colors
cell_type_colors <- colorRampPalette(c("seagreen4", "sandybrown", 
                                       "royalblue4", "paleturquoise3","purple4",
                                       "springgreen4","red4","peachpuff3","tan4","orchid4",
                                       "yellow4","orange","turquoise4","wheat3"))(length(unique(adata$obs$cell_type)))
#get a disimilarity map #attempt to get colors far apart which is not very successful really 
pos <- do.call(rbind, pblapply(levels(adata$obs$cell_type), function(x){
  Matrix::colMeans(adata$obsm$get("X_umap")[adata$obs$cell_type %in% x, ])
})) 
dmat <- 1-dist(pos) 
hcl <- hclust(dmat)
names(cell_type_colors) <- levels(adata$obs$cell_type)[hcl$order]
#refine some colors
cell_type_colors$CD8_Tem_proliferating <- "forestgreen"
cell_type_colors$CD4_Tfh <- "salmon2"
cell_type_colors$cDC2 <- "darkblue"
cell_type_colors$cDC1 <- "plum3"
cell_type_colors$activated_DC <- "orangered3"
#save celltype colors
saveRDS(cell_type_colors, file = "data/cell_type_colors.RDS")
#plot this - we can use this for figure 1..
annotated_raster_plot(layout = "X_umap", group = "cell_type",adata = adata) + 
  scale_color_manual(values = cell_type_colors) + theme_void() + theme(legend.position = "none")

#save adata ####
adata$write_h5ad("data/annotated_adata.h5ad")

#get cell type markers ####
celltype_markers <- tfidf_all_markers(adata, groups = "cell_type")

#save markers ####
saveRDS(celltype_markers, "data/celltype_markers.RDS")

#look at cell type proportions by disease ####
adata <- sc$read_h5ad("data/annotated_adata.h5ad")
adata$obs$disease_state <- as.character(adata$obs$disease)
adata$obs$disease_state[adata$obs$disease_state %in% c("HDLR", "HDMC","HDNLP","HDNS")] <- "cHL"
adata$obs$disease_state <- factor(adata$obs$disease_state, levels = c("healthy", "reactive_lymph_node", "cHL"))
freq_df <- data.frame(table( adata$obs$disease_state, adata$obs$cell_type)/rowSums(table( adata$obs$disease_state, adata$obs$cell_type)))
freq_df[freq_df$Freq == 0, "Freq"] <- NA
ggplot(freq_df, aes(x = Var1, y = Var2, size = Freq, fill = Var1)) + geom_point(pch=21) +  theme_bw() +
  scale_fill_manual(values = c("healthy" = "darkblue", "reactive_lymph_node" = "purple", "cHL" = "darkred")) + 
  scale_size_continuous(breaks = c(0.01, 0.05, 0.1, 0.2 ))

#get the MNP proportions in disease and health #####
adata <- sc$read_h5ad("data/annotated_adata.h5ad")
adata <- adata[adata$obs$cell_type %in% c("cDC2", "cDC1", "activated_DC", "macrophage", "classical_monocyte", "mast", "plasmacytoid_DC")]
adata$obs$disease_state <- as.character(adata$obs$disease)
adata$obs$disease_state[adata$obs$disease_state %in% c("reactive_lymph_node", "healthy")] <- "healthy"
adata$obs$disease_state[adata$obs$disease_state %in% c("HDLR", "HDMC","HDNLP","HDNS")] <- "cHL"
adata$obs$disease_state <- factor(adata$obs$disease_state, levels = c("healthy", "cHL"))
sub <- min(table(adata$obs$disease_state))
gs <- import("geosketch")
df <- do.call(rbind, lapply(unique(adata$obs$disease_state), function(x){
  adata_sub <- adata[adata$obs$disease_state %in% x]
  idx = gs$gs(adata_sub$obsm$get("X_scVI"), sub)
  df <- data.frame(adata_sub$obs[unlist(idx), c("disease_state", 'cell_type')])
  return(df)
}))
tab <- data.frame(table(df$cell_type, df$disease_state)/rowSums(table(df$cell_type, df$disease_state)))
tab$Var1 <- factor(tab$Var1, levels = c("plasmacytoid_DC", "mast", "activated_DC", "cDC1", "cDC2", "macrophage", "classical_monocyte"))
ggplot(tab, aes(x=Var2, y= Var1, fill = Freq)) + geom_tile() + scale_fill_gradientn(colours = c("white", "black"), limits = c(0, 1)) + theme_bw() + coord_fixed()

#get the disease proportions in disease and health #####
adata <- sc$read_h5ad("data/annotated_adata.h5ad")
sub <- min(table(adata$obs$disease))
gs <- import("geosketch")
df <- do.call(rbind, lapply(unique(adata$obs$disease), function(x){
  adata_sub <- adata[adata$obs$disease %in% x]
  idx = gs$gs(adata_sub$obsm$get("X_scVI"), sub)
  df <- data.frame(adata_sub$obs[unlist(idx), c("disease", 'cell_type')])
  return(df)
}))
tab <- data.frame(table(df$cell_type, df$disease)/rowSums(table(df$cell_type, df$disease)))
tab$Var1 <- factor(tab$Var1, levels = levels(adata$obs$cell_type))
tab$Var2 <- factor(tab$Var2, levels = c("healthy", 'reactive_lymph_node', 'HDLR', "HDMC", "HDNS", "HDNLP"))
tab <- tab[!tab$Freq == "NaN", ]
ggplot(tab, aes(x=Var2, y= Var1, fill = Freq)) + geom_tile() + scale_fill_gradientn(colours = c("white", "black"), limits = c(0, 1)) + 
  theme_bw() + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1))




