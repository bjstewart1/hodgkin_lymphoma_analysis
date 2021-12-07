# nanostring data analysis 
#functions #####
source("code/functions.R")
ad <- import("anndata")
sc <- import("scanpy")

#read in the normalised nanostring #####
nstrng <- read.csv("data/nanostring_data/HLN_nanostring_normalised.csv") #this is the normalised data from nanostring
nstrng <-nstrng[, 4:ncol(nstrng)]
colnames(nstrng) <- gsub("[.]", "-", colnames(nstrng))
nstrng <- nstrng[, !colnames(nstrng) %in% "Negative-probe"] #remove the negative probe
#get the annotations 
nstring_annotations <- read.csv("data/nanostring_data/nanostring_annotation.csv")
#plot positions on the slide
ggplot(nstring_annotations, aes(x = ROICoordinateX, y= ROICoordinateY, fill = PDL1.REGION)) + 
  geom_point(pch = 21, size = 3) + facet_wrap(~SlideName, nrow = 1) + theme_bw() + coord_fixed() + 
  scale_fill_manual(values= c("HIGH" = "red", "LOW" = "darkblue", "LOW_CTRL" = "grey50")) + theme(axis.text = element_blank(),
                                                                              axis.ticks = element_blank()) + 
  xlab("X coordinate") + ylab("Y coordinate")

ggplot(nstring_annotations, aes(x = ROICoordinateX, y= ROICoordinateY, fill = NOTE)) + 
  geom_point(pch = 21, size = 3) + facet_wrap(~SlideName, nrow = 1) + theme_bw() + coord_fixed() + theme(axis.text = element_blank(),
                                                                                                  axis.ticks = element_blank()) + 
  xlab("X coordinate") + ylab("Y coordinate")

#do some preprocessing - PCA #####
#make this an adata object
adata <- ad$AnnData(X = nstrng, obs = data.frame(nstring_annotations))
adata$obs_names <- nstring_annotations$SegmentDisplayName
adata$var_names <- colnames(nstrng)

#do pca
sc$pp$pca(adata)
#plot PCA
ggplot(data.frame(adata$obsm$get("X_pca"),adata$obs), aes(x=X1, y=X2, fill = PDL1.low)) + 
  geom_point(pch=21, size = 3) + scale_fill_manual(values = c("FALSE" = "grey50","TRUE" = "orange")) + 
  coord_fixed()
#construct a graph and cluster ####
#do shared nn graph 
knn_idx <- FNN::get.knn(adata$obsm$get("X_pca"))$nn.index
snn_el <- do.call(rbind, pblapply(1:nrow(knn_idx), function(x){
  nbrs <- knn_idx[x, ]
  nbr_list <- list()
  p = 1
  for(i in nbrs){
    nbr_1 <- nbrs 
    nbr_2 <- knn_idx[i, ]
    jdist <- length(intersect(nbr_1, nbr_2))/length(union(nbr_1, nbr_2))
    nbr_list[[p]] <- data.frame("node1" = x, "node2" = i, "jdist" = jdist)
    p = p+1
      }
  nbr_out <- do.call(rbind, nbr_list)
  return(nbr_out)
  
}))
gr <- graph_from_edgelist(as.matrix(snn_el[, 1:2]), directed = FALSE)
E(gr)$weight <- snn_el$jdist
plot(gr)
clusters <- leiden::leiden(gr)
cluster_colors <- paletteer::paletteer_d("nationalparkcolors::Arches", length(unique(clusters))) 
names(cluster_colors) <- unique(clusters)
#add clusters to the objects
adata$obs$cluster <- factor(clusters, levels = unique(clusters))

library(ggraph)
set.seed(100)
lay <- layout_with_fr(gr)
pl1 <- ggraph(gr, layout = lay) + geom_edge_link(edge_colour = "grey60", edge_width = 0.5) + 
  geom_node_point(pch = 21, aes(fill = factor(clusters)), size = 5) +
  theme_void() + 
  scale_fill_manual(values = cluster_colors) + coord_fixed()

pl2 <- ggraph(gr, layout = lay) + geom_edge_link(edge_colour = "grey60", edge_width = 0.5) + 
  geom_node_point(pch = 21, aes(fill = adata$obs$PDL1.REGION), size = 5) +
  theme_void()  + scale_fill_manual(values= c("HIGH" = "red", "LOW" = "darkblue", "LOW_CTRL" = "grey50")) + 
  coord_fixed() 

pl3 <- ggraph(gr, layout = lay) + geom_edge_link(edge_colour = "grey60", edge_width = 0.5) + 
  geom_node_point(pch = 21, aes(fill = adata$obs$NOTE), size = 5) +
  theme_void() + coord_fixed() 

pl4<- ggraph(gr, layout = lay) + geom_edge_link(edge_colour = "grey60", edge_width = 0.5) + 
  geom_node_point(pch = 21, aes(fill = factor(adata$obs$SlideName)), size = 5) +
  theme_void() + coord_fixed() 

pl5<- ggraph(gr, layout = lay) + geom_edge_link(edge_colour = "grey60", edge_width = 0.5) + 
  geom_node_point(pch = 21, aes(fill = factor(adata$obs$ScanLabel)), size = 5) +
  theme_void() + coord_fixed()

pl_list <- list(pl1, pl2, pl3, pl4, pl5)
pl_list <- lapply(pl_list, function(pl){
  pl <- pl +  theme(legend.position = "none")
  return(pl)
}) #remove the legends for sizing purposes
plot_grid(plotlist = pl_list, nrow = 2) 

#plot nuclei per cluster
ggplot(data.frame("nuclei_count" = nstring_annotations$AOINucleiCount,
                  "cluster" = adata$obs$cluster), aes(x=cluster, y = nuclei_count)) + geom_boxplot() + theme_bw()

#marker genes #####
#now get DE genes for each of these clusters..
#do limma
get_de_markers <- function(adata, cluster = "cluster"){
  de_list <- pblapply(as.character(levels(adata$obs[,cluster])), function(cl){
    library(limma)
    expr <- as.matrix(t(adata$X))
    gene_names <- rownames(expr) <- adata$var$Symbol
    ds <- adata$obs[, cluster] %in% cl
    groups <- as.integer(ds)
    groups <- data.frame(gps = groups)
    groups$gps[groups$gps %in% 0] <- "y"
    groups$gps[groups$gps %in% 1] <- "x"
    mm <- model.matrix(~0 + gps, data =groups) 
    fit <- lmFit(expr, mm)  
    contr <- makeContrasts(gpsx - gpsy, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrasts = contr)
    tmp <- eBayes(tmp)
    de <- topTable(tmp, sort.by = "P", n = nrow(expr)) # top 20 DE genes
    de$gene <- adata$var_names$values[as.numeric(rownames(de))]
    de <- de[order(-log10(de$P.Value)*sign(de$logFC), decreasing = TRUE), ]
    return(de)
  })
  names(de_list) <- levels(adata$obs[, cluster])
return(de_list)
}
cluster_markers <- get_de_markers(adata, cluster = "cluster")
cluster_markers <- cluster_markers[order(names(cluster_markers))]

#plot out the gene expression differences
cl <- as.character(adata$obs$cluster)
markers_use <- unique(unlist(lapply(cluster_markers, function(x){x[1:10, "gene"]}))) #we plot the top ten for each cluster
expr <- data.frame(t(adata$X), row.names = adata$var_names$values)[markers_use, order(cl)]
expr <- t(scale(t(expr)))
expr <- reshape2::melt(expr)
expr$Var1 <- factor(expr$Var1, levels = rev(markers_use))
ggplot(expr, aes(x=Var2, y = Var1, fill = value)) + ggrastr::geom_tile_rast() + coord_fixed() + 
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "grey90",midpoint = 0, 
                       limits = c(-5, 5)) + theme(axis.text.x = element_blank(), 
                                                  axis.text.y = element_text(face = "italic", color = "black", size = 5),
                                                  axis.ticks = element_blank(),
                                                  axis.title = element_blank())

#save the nanostring adata
adata$write("data/nanostring_adata.h5ad")

#now do GSVA on DE markers from microarray study #####
#now get DE markers from the microarray study
microarray_markers <- readRDS("data/HRS_DE.RDS")
microarray_markers <- list("HRSC" = microarray_markers[microarray_markers$significant %in% "HRS", "gene"],
                           "Germinal_centre" =  microarray_markers[microarray_markers$significant %in% "GC", "gene"])

#deconvolution of cell types #####
#nanostring deconvolution
#read in ns data
ns_adata <- sc$read_h5ad('data/nanostring_adata.h5ad')
#read in sc data
sc_adata <- sc$read_h5ad("data/annotated_adata.h5ad")
sc_adata$var_names <- sc_adata$var$Symbol
sc_adata$var_names_make_unique()
sc_adata$obs$broad_celltype <- as.character(sc_adata$obs$cell_type)
sc_adata$obs$broad_celltype[sc_adata$obs$broad_celltype %in% c("NK_1", "NK_2")] <- "NK_cell"
sc_adata$obs$broad_celltype[sc_adata$obs$broad_celltype %in% c("naive_b", "memory_b", "proliferating_b")] <- "b_cell"
sc_adata$obs$broad_celltype[sc_adata$obs$broad_celltype %in% c("CD4_Th0", "CD4_Th17", "CD4_Th1")] <- "CD4_T_cell"
sc_adata$obs$broad_celltype[sc_adata$obs$broad_celltype %in% c("CD8_Tem_terminal", "CD8_Tem", "CD8_Tem_proliferating")] <- "CD8_Tem"
sc_adata$obs$broad_celltype[sc_adata$obs$broad_celltype %in% c("nonfenestrated_endothelium", "lymphatic_endothelium", "fenestrated_endothelium")] <- "endothelium"
sc_adata$obs$broad_celltype[sc_adata$obs$broad_celltype %in% c("CD4_Texh_proliferating", "CD4_Texh")] <- "CD4_Texh"
#read in microarray data
microarray_data <- readRDS("data/GSE39133.RDS")

#get single cell markers
sc_markers <- tfidf_all_markers(sc_adata, "broad_celltype", gene_subset = sc_adata$var_names$values[sc_adata$var$highly_variable] )
HRS_signature <- readRDS( "data/HRSC_signature.RDS")
marker_list <-lapply(sc_markers, function(x){as.character(x[1:50, "gene"])}) #we take 50 markers here
marker_list[["HRSC"]] <- HRS_signature$gene

#create a binary matrix of marker genes membership
ns_genes <- ns_adata$var_names$values
binary_matrix <- do.call(cbind, lapply(marker_list, function(x){
  ns_genes %in% x
}))*1
rownames(binary_matrix) <- ns_genes
binary_matrix <- binary_matrix[!rowSums(binary_matrix) == 0, ] 
pheatmap(binary_matrix)
binary_mat_df <- reshape2::melt(binary_matrix)
binary_mat_df$Var1 <- factor(binary_mat_df$Var1, levels = rownames(binary_matrix)[hclust(dist(binary_matrix))$order])
ggplot(binary_mat_df, aes(x=Var2, y=Var1, fill = value %in% 1)) + geom_tile() + scale_fill_manual(values = c("TRUE" = 'black', 'FALSE' = "white")) +  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(color = "black",angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 4)) + xlab("") + ylab("")

#get the nanostring matrix
ns_matrix <- t(ns_adata$X)
rownames(ns_matrix) <- ns_adata$var_names$values
colnames(ns_matrix) <- ns_adata$obs_names$values
ns_matrix <- ns_matrix[rownames(binary_matrix), ] #subset to marker genes
#ns_matrix <- t(apply(ns_matrix, 1, min_max_normalisation))
#do the deconvolution using rlm (https://doi.org/10.1016/j.celrep.2019.01.041 https://github.com/giannimonaco/ABIS)
#rlm is meant to be less sensitive to outliers than lm 
est.res <- apply(ns_matrix, 2, function(bulk) {
  coefs <- coef(MASS::rlm(bulk~ binary_matrix, maxit = 10000))
  coefs[coefs < 0] <- 0 #make anything <0 = 0
  coefs <- coefs[2:length(coefs)] / sum(coefs[2: length(coefs)]) #fractions
  return(coefs)})
rownames(est.res) <- gsub("binary_matrix", "", rownames(est.res))
est.res <- t(est.res) #these sum to 1 and we will assume they are fractions 
rowSums(est.res)

#get nuclei counts
nstring_annotations <- read.csv("data/nanostring_data/nanostring_annotation.csv")
counts <- apply(est.res,2, function(x){round(x*nstring_annotations$AOINucleiCount)}) #we get the nuclei counts and * by the estimated fractions and round to get estimated cell numbers
counts <- counts[, !colSums(counts) == 0] #remove the two cell types that don't seem represented (CD8Tem and mast)
ag_expr <- data.frame(aggregate(counts, by = list(ns_adata$obs$cluster), mean), row.names=  1) #get the average per cluster

#then we can do differential expression of cell counts using ROI as replicates 
library("DESeq2")
da_results <- do.call(rbind, lapply(unique(ns_adata$obs$cluster), function(x){
  dds <- DESeqDataSetFromMatrix(countData = t(counts),
                                colData = data.frame("cluster" = ns_adata$obs$cluster %in% x), design = ~cluster )
  dds <- DESeq(dds)
  res <- results(dds)
  res$cluster <- x
  return(res)
}))
da_results$celltype <- rownames(da_results)
da_results <- da_results[!is.na(da_results$log2FoldChange), ]
da_results$cluster <- factor(da_results$cluster, levels = 1:5)
da_results <- as.data.frame(da_results)
#now get a matrix of log fold changes one against rest
lfc_matrix <- do.call(rbind, lapply(unique(da_results$celltype), function(x){
  inp <- da_results[da_results$celltype %in% x, ]
  inp$log2FoldChange
}))
rownames(lfc_matrix) <- unique(da_results$celltype)
colnames(lfc_matrix) <- unique(da_results$cluster)
da_results$celltype <- factor(da_results$celltype, levels = rownames(lfc_matrix)[hclust(dist(lfc_matrix))$order])
da_fills <-  scale_fill_gradient2(low = "darkblue", mid ="grey70", high = "darkred" , limits = c(-7, 7))
#and plot this
ggplot(da_results, aes(x = cluster, y = celltype, fill = log2FoldChange)) + geom_tile() + 
  scale_fill_gradient2(low = "darkblue", "mid" = "grey90", "high" = "red", midpoint = 0, 
                       limits = c(-max(abs(da_results$log2FoldChange)), max(abs(da_results$log2FoldChange))) ) + theme_bw() + coord_fixed() + scale_size_continuous(range = c(0.5,3))






