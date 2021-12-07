#Functions
library(scran)
library(scater)
library(cowplot)
library(reticulate)
use_python("/opt/conda/bin/python")
library(pheatmap)
library(DropletUtils)
library(igraph)
library(SoupX)
library(pbapply)
library(paletteer)

#helpful themes
figure_theme <- theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
                      axis.line = element_blank(),
                      axis.text.x  = element_text(color = "black", size = 7.5),
                      axis.text.y  = element_text(color = "black", size = 7.5),
                      axis.title.x  = element_text(color = "black", size = 10),
                      axis.title.y  = element_text(color = "black", size = 10)
)
figure_colors <- paletteer::scale_color_paletteer_d("ggsci", "default_igv")
figure_fills <- paletteer::scale_fill_paletteer_d("ggsci", "default_igv")

# Get density of points in 2 dimensions.
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @return The density within each square.
get_density <- function(x, y, n = 100) {
  require(MASS)
  dens <- kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Plot the x-y density with ggplot2 raster
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param title Title for the plot.
#' @return A plot of x and y colored by density.
density_plot_raster <- function(x, y, xlab = NULL, ylab = NULL, title = NULL, n = 100, colorscale = viridis::viridis(100)){
  require(ggrastr)
  dat <- data.frame(xlab = x, ylab = y)
  dat$density <- get_density(dat$x, dat$y, n=n)
  ggplot(dat, aes(x = x, y = y, col = density)) + ggrastr::geom_point_rast(pch=19, cex=0.05) + 
    scale_color_gradientn(colours = colorscale)+ xlab(xlab) + ylab(ylab) + ggtitle(title) + theme_bw() + theme(axis.text = element_text(color = "black"))
}

#function to get average coordinates based on groups 
#' @param layout  - a layout to use, be it PCA, tSNE, UMAP, graph layout
#' @param groups - a factor to group by
average.coords <- function(layout, groups){
  groups <- as.factor(groups)
  dat <- data.frame("X" = layout[, 1], "Y" = layout[, 2])
  centroids <- do.call(rbind, lapply(levels(groups), function(x){
    cm <- colMeans(dat[groups == x, ])
    return(cm)
  }))
  rownames(centroids) <- levels(groups)
  return(centroids)
}

#annotated raster plot
#' @param layout the layout to use from the adata object
#' @param adata - the adata object to use
#' @param group - the grouping to use in adata$obs
annotated_raster_plot <- function(layout = "X_umap", group = "leiden", adata, xlab = "Dim1", ylab = "Dim2", size = 0.05, dpi = 600){
  require(ggrastr)
  group <- adata$obs[, group]
  layout <- adata$obsm$get(layout)
  av.coords  <- average.coords(layout, groups = group)
  set.seed(100)
  idx <- sample(1:nrow(layout))
  ggplot(data.frame("x" = layout[idx, 1], "y" = layout[idx, 2], "Factor" = as.factor(group[idx])), 
         aes(x = x, y=y, col = Factor)) + ggrastr::geom_point_rast(pch=19, size=size, raster.dpi = dpi) + xlab(xlab) +ylab(ylab)+ theme_classic() + 
    theme(legend.position="none") + annotate("text", x = av.coords[, 1], y=av.coords[, 2],label= rownames(av.coords))
}

# raster plot
#' @param layout the layout to use from the adata object
#' @param adata - the adata object to use
#' @param group - the grouping to use in adata$obs
raster_plot <- function(layout = "X_umap", group = "leiden", adata, xlab = "Dim1", ylab = "Dim2", size = 0.05, dpi = 600){
  require(ggrastr)
  require(ggrastr)
  group <- adata$obs[, group]
  layout <- adata$obsm$get(layout)
  set.seed(100)
  idx <- sample(1:nrow(layout))
  ggplot(data.frame("x" = layout[idx, 1], "y" = layout[idx, 2], "Factor" = as.factor(group[idx])), 
         aes(x = x, y=y, col = Factor)) + ggrastr::geom_point_rast(pch=19, size=size, raster.dpi = dpi) + xlab(xlab) +ylab(ylab)+ theme_classic() + 
    guides(colour = guide_legend(override.aes = list(size=5)))
}

#min max normalisation
#' @param x - a vector of values
min_max_normalisation <- function(x){
  return((x - min(x))/(max(x) -min(x)))
}

#tfidf markers

#' Calculates the tf-idf for a set of target cells against a background of the cells given in "universe".
#' From Matt Youngs script for kidney cells. 10.1126/science.aat1699
#'
#' @param data The data matrix to use.
#' @param target Columns that are the target.
#' @param universe Columns that we should consider (target must be a subset).
tfidf = function(data,target,universe){
  require(Matrix)
  if(!all(target %in% universe))
    stop('Target must be a subset of universe')
  nObs = Matrix::rowSums(data[,colnames(data) %in% target]>0)
  nTot = Matrix::rowSums(data[,colnames(data) %in% universe]>0)
  tf = nObs/length(target)
  idf = log(length(universe)/nTot)
  score = tf*idf
  #Calculate p-value for significance based on using a hypergeometric distribution to simulate the results of infinite random sampling
  pvals = phyper(nObs-1,nTot,length(universe)-nTot,length(target),lower.tail=FALSE)
  qvals = p.adjust(pvals,method='BH')
  ntf = (exp(-idf)*length(universe)-tf*length(target))/(length(universe)-length(target))
  return(data.frame(geneFrequency=tf,
                    geneFrequencyOutsideCluster=ntf,
                    geneFrequencyGlobal=exp(-idf),
                    geneExpression=Matrix::rowMeans(data[,colnames(data) %in% target]),
                    geneExpressionOutsideCluster = Matrix::rowMeans(data[, !colnames(data) %in% target]),
                    geneExpressionGlobal = Matrix::rowMeans(data),
                    idf=idf,
                    tfidf=score,
                    qval=qvals)[order(score,decreasing=TRUE),])
}



#find all markers with tfidf
#' @param adata - an adata object
#' @param groups - a set of cluster identities in the adata
#' @param gene_subset - genes to subset - HVG
#' @param use_raw - boolean use raw?
tfidf_all_markers <- function(adata, groups = "leiden", gene_subset = NULL, use_raw = FALSE){
  groups <- adata$obs[, groups]
  groups <- factor(groups)
  
  if(use_raw){
    adata <- adata$raw
  }else{adata <-adata}
  if(is.null(gene_subset)){
    gene_subset <- adata$var_names$values
    message("using all genes")
  }else{
    message("subsetting some genes")
  }
  #matrix organisation
  np <- import("numpy")
  mat_in <- np$transpose(adata$X)
  cellnames <- adata$obs_names$values
  colnames(mat_in) <-  cellnames
  rownames(mat_in) <- adata$var_names$values
  #then subset
  mat_in <- mat_in[rownames(mat_in) %in% gene_subset, ]
  #library(pbmcapply)
  tfid.list <-  pblapply(levels(groups), function(i){
    tfid <- tfidf(data = mat_in, target = colnames(mat_in)[groups == i], universe = colnames(mat_in))
    tfid$gene <- adata$var[rownames(tfid), "Symbol"]
    return(tfid)
  })
  names(tfid.list) <- as.character(levels(groups))
  return(tfid.list)    
}


#scrublet wrapper
#' @adata the anndata object
#' @channel the obs column that corresponds to the 10X channel. You could potentially use "donor" here...
scrublet_wrapper <- function(adata, channel = "channel"){
  scrb <- import("scrublet")
  scrb_results <- pblapply(unique(adata$obs[, channel]), function(x){
    sub_adata <- adata[adata$obs[, channel] %in% x]
    scrub = scrb$Scrublet(sub_adata$X)
    scrub_out = scrub$scrub_doublets()
    return(scrub_out)
  })
  adata$obs$scrublet_scores <-  unlist(lapply(scrb_results, function(x){x[1]}))
  adata$obs$scrublet_cuts <-  unlist(lapply(scrb_results, function(x){x[2]}))
  return(adata)
}


#function to recluster existing leiden clusters on a graph.
#' @param adata an anndata object
#' @param clusters the cluster to recluster
#' @param resolution the resolution to recluster with
recluster_leiden <- function(adata, clusters, resolution = 0.2){
  clusters <- as.list(clusters)
  sc$tl$leiden(adata, 
               restrict_to = 
                 tuple("leiden", clusters), 
               resolution = as.numeric(resolution))
  adata$obs$leiden <- as.character(adata$obs$leiden)
  adata$obs$leiden_R <- as.character(adata$obs$leiden_R)
  adata$obs$leiden[adata$obs$leiden %in% clusters ] <- adata$obs$leiden_R[adata$obs$leiden %in% clusters]
  adata$obs$leiden <- as.factor(adata$obs$leiden)
  return(adata)
}

