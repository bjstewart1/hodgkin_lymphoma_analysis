#import functions #####
source("code/functions.R")
ad <- import("anndata")
sc <- import("scanpy")
# read in the stiedl lymph node data ####
load("data/merged_sc_obj_CancerDiscovery.Rdata")
sce <- sc.obj
#now make an adata object... 
adata <- ad$AnnData(t(counts(sce)))
adata$var_names <- rowData(sce)$id
adata$var$Symbol <- rownames(rowData(sce))
adata$var$ID <- rowData(sce)$id
adata$obs_names <- rownames(colData(sce))
adata$obs$empty_drops = TRUE
adata$obs$dd = TRUE
adata$obs$channel = colData(sce)$dataset
adata$obs$donor = colData(sce)$dataset
adata$obs$donor[adata$obs$donor %in% c("RLN-1_R1", "RLN-1_R2")] <- "RLN-1"
adata$obs$anatomy = "lymph_node"
colData(sce)$Pathology[is.na(colData(sce)$Pathology)] <- "reactive_lymph_node"
adata$obs$disease = colData(sce)$Pathology
adata$obs$centre = "bc_cancer"
adata$obs$cell_selection <- "unselected"
adata$obs$tenx_run <- sce$Chip
bc_hln <- adata
rm(adata)


# read in Sanger lymph node ####
dir <- "/nfs/team297/bs16/human_10X_data/lymph_node"
files <- list.files(dir)
#this will read in each sample, do default drops
adata_list <- pblapply(files, function(f){
  #this will read in each sample, do default drops
  message("sample ", f)
  adata <- sc$read_h5ad(file.path(dir, f, "busout", "adata.h5ad"))
  set.seed(100)
  e.out <- emptyDrops(t(adata$X))
  is.cell <- e.out$FDR <= 0.01
  is.cell <- is.cell %in% "TRUE"
  sum(is.cell)
  adata$obs$empty_drops <- is.cell
  #do default drops
  adata$obs$dd <- DropletUtils::defaultDrops(t(adata$X))
  adata <- adata[adata$obs$dd]
  adata$obs$channel <- f #add the channel name
  adata$obs_names <- paste0(adata$obs_names$values, "_", f) #and make the obs names unique
  return(adata)
})
adata <- sc$concat(adata_list)
#now get symbols for the genes
t2g <- read.table( "/nfs/team297/bs16/tools/cellranger_bus_human_ref/non_velocity/t2g.txt")
t2g <- t2g[!duplicated(t2g$V2), ]
table(t2g$V2 == adata$var_names$values)
adata$var$ID <- adata$var_names$values
adata$var$Symbol <- as.character(t2g$V3)
sanger_ln <- adata
rm(adata)

#read in Sanger HLN ####
dir <- "/nfs/team297/bs16/human_10X_data/Hodgkins_LN"
files <- list.files(dir)
files <- files[2:9]
#this will read in each sample, do default drops
adata_list <- pblapply(files, function(f){
  #this will read in each sample, do default drops
  message("sample ", f)
  adata <- sc$read_h5ad(file.path(dir, f, "busout", "adata.h5ad"))
  set.seed(100)
  e.out <- emptyDrops(t(adata$X))
  is.cell <- e.out$FDR <= 0.01
  is.cell <- is.cell %in% "TRUE"
  sum(is.cell)
  adata$obs$empty_drops <- is.cell
  #do default drops
  adata$obs$dd <- DropletUtils::defaultDrops(t(adata$X))
  adata <- adata[adata$obs$dd]
  adata$obs$channel <- f #add the channel name
  adata$obs_names <- paste0(adata$obs_names$values, "_", f) #and make the obs names unique
  return(adata)
})
adata <- sc$concat(adata_list)
#now get symbols for the genes
t2g <- read.table( "/nfs/team297/bs16/tools/cellranger_bus_human_ref/non_velocity/t2g.txt")
t2g <- t2g[!duplicated(t2g$V2), ]
table(t2g$V2 == adata$var_names$values)
adata$var$ID <- adata$var_names$values
adata$var$Symbol <- as.character(t2g$V3)
sanger_hln <- adata
rm(adata)

# sort out metadata #####
#sanger lymph node
sanger_ln_meta <- read.csv("data/sanger_LN_meta.csv")
sanger_ln_meta$Sanger.Sample.ID <- gsub("-", "_", sanger_ln_meta$Sanger.Sample.ID)
idx <- match(sanger_ln$obs$channel, sanger_ln_meta$Sanger.Sample.ID)
sanger_ln_meta <- sanger_ln_meta[idx, ]
sanger_ln$obs$disease = "healthy_lymph_node"
sanger_ln$obs$anatomy = sanger_ln_meta$organ.code
sanger_ln$obs$disease = "healthy"
sanger_ln$obs$centre = "WSI"
sanger_ln$obs$donor <- sanger_ln_meta$Source.ID
sanger_ln$obs$tenx_run <- sanger_ln$obs$donor #donors are run as one batch
#sort out levels of tissue processing
ln_processing <- sanger_ln_meta$explanation.of.processing
ln_processing[ln_processing %in% grep("total", ln_processing,value = TRUE, ignore.case =TRUE)] <- "unselected"
ln_processing[ln_processing %in% "Single cells for 10X"] <- "unselected"
ln_processing[ln_processing %in% grep("CD4 T cells",ln_processing, value = TRUE )] <- "CD4+ T cells"
ln_processing[ln_processing %in% grep("all other immune cells",ln_processing, value = TRUE )] <- "CD45+_non_T_cells"
sanger_ln$obs$cell_selection <- ln_processing

#sanger hodgkins lymph node
sanger_hln$obs$disease = "HDNS"
sanger_hln$obs$disease[sanger_hln$obs$channel %in% c("4710STDY7018928",
                                           "4710STDY7018929",
                                           "4710STDY7018930",
                                           "4710STDY7018931")] <- "HDNLP"
sanger_hln$obs$donor <- "WSI_HLN_donor_1"
sanger_hln$obs$donor[sanger_hln$obs$channel %in% c("4710STDY7018928",
                                         "4710STDY7018929",
                                         "4710STDY7018930",
                                         "4710STDY7018931")] <- "WSI_HLN_donor_2"
sanger_hln$obs$anatomy = "Hodgkin_LN"
sanger_hln$obs$centre = "WSI"
sanger_hln$obs$tenx_run <- sanger_hln$obs$donor #donors are run as one batch
sanger_hln$obs$cell_selection <- "unselected"

#sort out var for now #####
#now get symbols for the genes
t2g <- read.table( "/nfs/team297/bs16/tools/cellranger_bus_human_ref/non_velocity/t2g.txt")
t2g <- t2g[!duplicated(t2g$V2), ]
sanger_hln$var$ID <- t2g$V2
sanger_hln$var$Symbol <- t2g$V3
sanger_ln$var$ID <- t2g$V2
sanger_ln$var$Symbol <- t2g$V3

sanger_hln$var_names <- sanger_hln$var$Symbol
sanger_hln$var_names_make_unique()
sanger_ln$var_names <- sanger_ln$var$Symbol
sanger_ln$var_names_make_unique()
bc_hln$var_names <- bc_hln$var$Symbol
bc_hln$var_names_make_unique()

#save the raw objects.. #####
sanger_ln$write_h5ad("data/sanger_ln_raw.h5ad")
sanger_hln$write_h5ad("data/sanger_hln_raw.h5ad")
bc_hln$write_h5ad("data/bc_hln_raw.h5ad")

#read back the raw objects  ####
sanger_hln <- sc$read_h5ad("data/sanger_hln_raw.h5ad")
#sanger_hln$X <- sanger_hln$layers$get("corrected_counts")
sanger_ln <- sc$read_h5ad("data/sanger_ln_raw.h5ad")
#sanger_ln$X <- sanger_ln$layers$get("corrected_counts")
bc_hln <- sc$read_h5ad("data/bc_hln_raw.h5ad")

#harmonise var names #####
sanger_hln$var_names <- sanger_hln$var$ID <- unlist(pblapply(strsplit(as.character(sanger_hln$var$ID), "[.]"), function(x){x[1]}))
sanger_ln$var_names <- sanger_ln$var$ID <- unlist(pblapply(strsplit(as.character(sanger_ln$var$ID), "[.]"), function(x){x[1]}))
bc_hln$var_names <- bc_hln$var$ID
#get the intersection of these var names
var_names <- list(bc_hln$var_names$values, sanger_ln$var_names$values, sanger_hln$var_names$values)
var_intersect <- Reduce(intersect, var_names)
#now subset each object to these genes
sanger_hln <- py_eval("r.sanger_hln[:, r.var_intersect]")
sanger_ln <- py_eval("r.sanger_ln[:, r.var_intersect]")
bc_hln <- py_eval("r.bc_hln[:, r.var_intersect]")

#concatenate the data ####
concatenated_adata <- sc$concat(list(sanger_ln, sanger_hln, bc_hln))
concatenated_adata$layers <- list("counts" = concatenated_adata$X) #we use corrected counts for everything
concatenated_adata$var <- sanger_ln$var #borroww var
adata <- concatenated_adata
#do some global QC for now #####
adata$var$MT = 1:nrow(adata$var) %in% grep("^MT-", adata$var$Symbol)
sc$pp$calculate_qc_metrics(adata, qc_vars= list("MT"), inplace = TRUE)

n_gene_cut = 500
pct_mito_cut = 10
top_50_cut = 65
total_counts_cut = 600

ggplot(adata$obs, aes(x=  n_genes_by_counts)) + geom_histogram(bins = 200) + theme_bw() + 
  geom_vline(xintercept = n_gene_cut, color = 'red') + 
  theme(axis.text = element_text(color = "black")) + facet_wrap(~centre)

ggplot(adata$obs, aes(x=  pct_counts_MT)) + geom_histogram(bins = 200) + theme_bw() + 
  geom_vline(xintercept = pct_mito_cut, color = 'red') +
  theme(axis.text = element_text(color = "black")) + facet_wrap(~centre)

#density plot of genes by MT counts
dp1 <- density_plot_raster(x = adata$obs$log1p_n_genes_by_counts, y= adata$obs$total_counts_MT) + 
  theme_bw() + xlab("log1p_n_genes_by_counts") + 
  ylab("log1p_total_counts_MT") + geom_vline(xintercept = log1p(n_gene_cut), color = "red")

#density plot of genes by pct counts in top 50 genes
dp2 <- density_plot_raster(x = adata$obs$log1p_n_genes_by_counts, y= adata$obs$pct_counts_in_top_50_genes) + 
  theme_bw() + xlab("log1p_n_genes_by_counts") + 
  ylab("pct_counts_in_top_50_genes") + 
  geom_vline(xintercept = log1p(n_gene_cut), color = "red") + geom_hline(yintercept = top_50_cut, color = "pink") + 
  ylim(c(0, 100))

#density plot of genes by total counts
dp3 <- density_plot_raster(x = adata$obs$log1p_n_genes_by_counts, y= adata$obs$log1p_total_counts) + 
  theme_bw() + ylab("log1p_total_counts") + 
  xlab("log1p_n_genes_by_counts") + 
  geom_vline(xintercept = log1p(n_gene_cut), color = "red") + 
  geom_hline(yintercept = log1p(total_counts_cut), color = "red") 

plot_grid(plotlist= list(dp1,dp2,dp3), nrow = 1)

#get lists of bad cells
bad_cells <- list(adata$obs_names$values[adata$obs$total_counts < total_counts_cut],
                  adata$obs_names$values[adata$obs$n_genes_by_counts < n_gene_cut],
                 # adata$obs_names$values[adata$obs$pct_counts_in_top_50_genes > top_50_cut],
                  adata$obs_names$values[adata$obs$pct_counts_MT] > pct_mito_cut)
adata$obs$qc_cuts <- !adata$obs_names$values %in% Reduce(union, bad_cells)

adata <- adata[adata$obs$qc_cuts]

#run scrublet doublet detection on a per channel basis #####
adata <- scrublet_wrapper(adata, channel = "channel")
#do a grep for nuissance genes ####
hkGeneREGEX='^(DNAJ[ABC]|EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN)'
coreExcludeGenes = unique(c(grep('\\.[0-9]+',adata$var$Symbol,value=TRUE), #Poorly characterised
                            grep('MALAT1',adata$var$Symbol,value=TRUE), #Contamination or highly expressed poorly characterised
                            grep('^MT-',adata$var$Symbol,value=TRUE), #Mitochondria
                            grep("XIST", adata$var$Symbol, value = TRUE), #F gender
                            grep(hkGeneREGEX,adata$var$Symbol,value=TRUE) #Housekeeping genes
))
adata$var$core_exclude <- !adata$var$Symbol %in% coreExcludeGenes


#write out preprocessed data #####
adata$write_h5ad("data/preprocessed_adata.h5ad")

