#differential abundance analysis
#import functions #####
source("code/functions.R")
ad <- import("anndata")
sc <- import("scanpy")

#sort out some colors #####
library(rcartocolor)
display_carto_all(type = 'diverging', colorblind_friendly = TRUE)
#carto <- rcartocolor::carto_pal(n = 7, "ArmyRose")
#da_colors <-  scale_color_gradient2(low = carto[1], mid ="grey90", high = carto[7] , limits = c(-7, 7))
da_colors <-  scale_color_gradient2(low = "darkblue", mid ="grey90", high = "darkred" , limits = c(-7, 7))
da_fills <-  scale_fill_gradient2(low = "darkblue", mid ="grey70", high = "darkred" , limits = c(-7, 7))
#make an mobj ####
library(miloR)
#read in the data
adata <- sc$read_h5ad("data/annotated_adata.h5ad")
#make an sce and convert to milo
sce <- SingleCellExperiment(assays = list('logcounts' = t(adata$X)), 
                            reducedDims = list("PCA" = adata$obsm$get("X_pca"),
                                               "SCVI" = adata$obsm$get("X_scVI"),
                                               "UMAP" = adata$obsm$get("X_umap")),
                            colData = data.frame(adata$obs),
                            rowData = data.frame(adata$var))
rownames(sce) <- rowData(sce)$Symbol #make rownames into symbols 
mobj <- Milo(sce)
saveRDS(mobj, "data/da_analysis_mobj.RDS")
#geosketch to downsample ####
mobj <- readRDS("data/da_analysis_mobj.RDS")
#gs <- import("geosketch")
#gs_idx <- gs$gs(X =reducedDim(mobj, "SCVI"), N = as.integer(75000))
#gs_idx <- 1:ncol(mobj) %in% unlist(gs_idx)
#mobj <- mobj[, gs_idx] #this downsamples it to a manageable size for memory.

#build KNN graph ####
#build the knn  graph using the milo method
mobj <- buildGraph(mobj, k = 30, d = 15, reduced.dim = "SCVI")
#define neighbourhoods
mobj <- makeNhoods(mobj, reduced_dims = "SCVI", d = 15)
#n size should be 50-100
plotNhoodSizeHist(mobj) + theme(axis.text = element_text(color = "black", size = 5),
                                axis.title = element_text(color = "black",size = 5))
mobj <- countCells(mobj, meta.data = colData(mobj), samples="channel")
head(nhoodCounts(mobj))
#now calculate neighbourhood distances.
mobj <- calcNhoodDistance(mobj, d=15, reduced.dim = "scvi")

#do differential abundance analysis ####
#set "disease state"
colData(mobj)$disease_state <- "healthy"
colData(mobj)$disease_state[colData(mobj)$disease %in% c("HDLR", "HDMC", 
                                                         "HDNLP",  
                                                         "HDNS")] <- "cHL"
colData(mobj)$disease_state <- factor(colData(mobj)$disease_state)
design <- data.frame(colData(mobj))[,c("channel", "centre", "disease_state")]
design <- dplyr::distinct(design)
design$disease_state <- factor(design$disease_state, levels = c("healthy", "cHL"))
rownames(design) <- design$channel
da_results <- testNhoods(mobj, design = ~centre + disease_state,
                         design.df = design, norm.method="TMM")
da_results <- annotateNhoods(mobj, da_results, coldata_col = "cell_type")
ggplot(da_results, aes(cell_type_fraction)) + geom_histogram(bins=50)
#annotate mixed cells neighbourhoods
da_results$cell_type <- ifelse(da_results$cell_type_fraction < 0.6, "Mixed", da_results$cell_type)

#plot this
ggplot(da_results, aes(PValue)) + 
  geom_histogram(bins=50, fill = "forestgreen") + theme_minimal() + 
  theme(axis.text = element_text(color = "black"))
ggplot(da_results, aes(x = logFC, y = -log10(PValue))) + theme_minimal() + 
  geom_point() + theme(axis.text = element_text(color = "black")) + geom_point(data = da_results[da_results$cell_type %in% "Plasmacytoid_DC", ],
                                                                               aes(x=logFC, y = -log10(PValue)), size = 2, color  = "red")


#plot the neighborhoods ####
plotNhood <- function(x, redDim = "UMAP", color_by = LFC, alpha = 1, 
                      size_scale = c(0.1, 5)){
  layout <- reducedDim(x, redDim)[unlist(x@nhoodIndex), ]
  size <- colSums(x@nhoods)
  df <- data.frame(x= layout[, 1], y = layout[, 2], size = size)
  pl <- ggplot(df, aes(x = x, y=y, 
                       size = size, fill = size)) + geom_point(pch = 21, alpha = alpha, color ="black") +
    theme_void() + scale_size(range = size_scale, name = "Nhood size")
  return(pl)
}
plotNhoodDA <- function(x, redDim = "UMAP",da_results, alpha=1, size_scale = c(0.1, 5),
                        color_scheme = scale_color_gradient2(low = "darkblue", mid ="grey90", high = "darkred" , limits = c(-6, 6))){
  layout <- reducedDim(x, redDim)[unlist(x@nhoodIndex), ]
  size <- colSums(x@nhoods)
  df <- data.frame(x= layout[, 1], y = layout[, 2], size = size)
  LFC <- da_results$logFC
  df <- data.frame(x= layout[, 1], y = layout[, 2], size = size, 
                   LFC = LFC)
  pl <- ggplot(df, aes(x = x, y=y, 
                       size = size, color = LFC)) + ggrastr::geom_point_rast(pch = 19, alpha = alpha) +
    theme_void() + scale_size(range = size_scale, name = "Nhood size") + color_scheme
  return(pl)
}
plotNhood(mobj, alpha = 0.7, size_scale = c(0.1, 3))
plotNhoodDA(mobj,redDim = "UMAP",  da_results, alpha=0.7, size_scale = c(0.1, 3), color_scheme = da_colors)

#plot the celltypes ####
#order these by average LFC
avLFC <- unlist(lapply(unique(da_results$cell_type), function(x){
  mean(da_results[da_results$cell_type %in% x, "logFC"])
}))
names(avLFC) <- unique(da_results$cell_type)
avLFC <- avLFC[order(avLFC)]
#then set this as a factor
da_results$cell_type <- factor(da_results$cell_type, levels = names(avLFC))
#then plot according to this factor
plotDAbeeswarm(da_results, group.by = "cell_type",
) + theme(axis.text = element_text(size =5, color = "black"),
          axis.title = element_text(size = 5, color = "black")) + da_colors + coord_polar()
#make a per-cell plot #####
#do this handmade version
da_plot <- da_results 
da_plot$size <- colSums(mobj@nhoods)
da_plot$is_sign <- da_plot$SpatialFDR < 0.05
da_plot$cell_type <- factor(da_plot$cell_type, levels =  rev(levels(da_plot$cell_type)))
da_plot <- da_plot[!da_plot$cell_type %in% "Mixed", ] #remove the "mixed" neighborhoods for plot clarity
ggplot(da_plot, aes(x= cell_type, y = logFC, 
                    color = logFC)) + 
  ggrastr::geom_quasirandom_rast(pch = 19, size = 0.1) +
  da_colors + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#get the celltype colors and plot with these colors
cell_type_colors <- readRDS("data/cell_type_colors.RDS")
cell_type_colors$mixed <- "grey50"
ggplot(da_plot, aes(x = logFC, y = cell_type, fill = cell_type) ) +geom_violin(scale= 'width') + theme_bw() + 
  scale_fill_manual(values = unlist(cell_type_colors)) + theme(legend.position = "none") 

#save objects #####
saveRDS(da_results, "data/differential_abundance_results.RDS")
saveRDS(mobj, "data/da_analysis_mobj.RDS")
