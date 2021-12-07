#import functions #####
source("code/functions.R")
ad <- import("anndata")
sc <- import("scanpy")

#import adata object #####
adata <- sc$read_h5ad("data/annotated_adata.h5ad")
#subset to myeloid #####
adata <- adata[adata$obs$cell_type %in% c("cDC2", "activated_DC", "cDC1", "plasmacytoid_DC", "macrophage", "classical_monocyte", "mast")]
#myeloid coexpression ####
genes_plot <- c('ITGAX', "CD1C", "IL3RA", "CADM1", "LAMP3", "CD68")
long_genes <- c('HLA-DRA', 'CD52', 'ITGAX', 'ITGAM', 'CPA3', 'TPSAB1', 'KIT',
'CLEC4C', 'IL3RA', 'IRF4', 'CXCR3','S100A8', 'S100A9', 'VCAN', 'FCN1',
'CD14','C1QA', 'C1QC', 'FOLR2', 'MRC1', 'CD68', 'FCGR3A',
'FSCN1', 'CCL17', 'LAMP3', 'CCR7',
'FCER1A', 'CLEC10A', 'CD1C',
'CLEC9A', 'IDO1', 'XCR1', 'CADM1')
coex_df <- pblapply(as.character( unique(adata$obs$cell_type)), function(i){
  message(paste0("celltype", i))
  cnames <- adata$obs_names$values[adata$obs$cell_type %in% i]
  c <- 1
  coex_list <- list()
  for(n in genes_plot){
    message(n)
    for(p in genes_plot){
      n_vec <- cnames[adata[adata$obs$cell_type %in% i]$X[, which(adata$var$Symbol == n)] > 0]
      p_vec <- cnames[adata[adata$obs$cell_type %in% i]$X[, which(adata$var$Symbol == p)] > 0]
      fraction_coexp <- length(intersect(n_vec, p_vec))/sum(adata$obs$cell_type %in% i)
      coex_list[[c]] <- c("gene_1" = n, "gene_2" = p, "fraction_coexp" = as.numeric(fraction_coexp), "celltype" = i)
      c = c+1
    }
  }
  coex_df <- do.call(rbind, coex_list)
  return(coex_df)
})
coex_df <- do.call(rbind, coex_df)
coex_dummy  <- do.call(rbind, lapply(as.character( unique(adata$obs$cell_type)), function(i){
  grd <- expand.grid(long_genes[!long_genes %in% genes_plot],long_genes[!long_genes %in% genes_plot])
  return(data.frame("gene_1" = as.character(grd[, 1]), "gene_2" =as.character( grd[, 2]), "fraction_coexp" = as.character(NA), "celltype" = i))
}))

full_coex_df <- data.frame(do.call(rbind, list(coex_df, coex_dummy)))
full_coex_df$fraction_coexp <- as.numeric(full_coex_df$fraction_coexp)
coex_df_plot <- full_coex_df
celltype_colors <- readRDS("data/cell_type_colors.RDS")
coex_df_plot$gene_2 <- factor(coex_df_plot$gene_2, levels = genes_plot)
coex_df_plot$gene_1 <- factor(coex_df_plot$gene_1, levels = long_genes)
ggplot(coex_df_plot, aes(x = gene_1, y = gene_2, fill = celltype,  size = fraction_coexp)) + geom_point(shape = 22) + theme_bw()  + 
  theme(axis.text.x =  element_text(angle = 45, hjust =1)) + scale_fill_manual(values = celltype_colors)  + 
  scale_size(limits =  c(0.2,1)) +   facet_wrap(~celltype, ncol = 1 )










