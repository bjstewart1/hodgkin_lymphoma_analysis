#cellphoneDB
#import functions #####
source("code/functions.R")
source("code/cpdb_functions.R")
ad <- import("anndata")
sc <- import("scanpy")
#write out objects for cellphonedb input #####
cpdb_adata <- sc$read_h5ad("data/annotated_adata.h5ad")
#subset to variable genes 
cpdb_adata <- py_eval("r.cpdb_adata[:, r.cpdb_adata.var.highly_variable]")
#now we rationalise the celltype labels
cpdb_adata$obs$broad_celltype <- as.character(cpdb_adata$obs$cell_type)
cpdb_adata$obs$broad_celltype[cpdb_adata$obs$broad_celltype %in% c("NK_1", "NK_2")] <- "NK_cell"
cpdb_adata$obs$broad_celltype[cpdb_adata$obs$broad_celltype %in% c("naive_b", "memory_b", "proliferating_b")] <- "b_cell"
cpdb_adata$obs$broad_celltype[cpdb_adata$obs$broad_celltype %in% c("CD4_Th0", "CD4_Th17", "CD4_Th1")] <- "CD4_T_cell"
cpdb_adata$obs$broad_celltype[cpdb_adata$obs$broad_celltype %in% c("CD8_Tem_terminal", "CD8_Tem", "CD8_Tem_proliferating")] <- "CD8_Tem"
cpdb_adata$obs$broad_celltype[cpdb_adata$obs$broad_celltype %in% c("nonfenestrated_endothelium", "lymphatic_endothelium", "fenestrated_endothelium")] <- "endothelium"
cpdb_adata$obs$broad_celltype[cpdb_adata$obs$broad_celltype %in% c("CD4_Texh_proliferating", "CD4_Texh")] <- "CD4_Texh"

table(cpdb_adata$obs$broad_celltype)

cpdb_adata$write_h5ad("data/cpdb_adata.h5ad")
#write out the metadata
cpdb_meta_data <- cbind(cpdb_adata$obs_names$values, as.character(cpdb_adata$obs$broad_celltype))   #####  cluster is the user’s specific cluster column
#subset just to cells we are interested in
write.table(cpdb_meta_data, 'data/cpdb_meta.txt', sep='\t', quote=F, row.names=F)

#run cellphone db
system("bash code/run_cpdb.sh")

#plot the cpdb results #####
cells_use <- c("cDC1","cDC2", "classical_monocyte", "macrophage","CD4_Texh", "CD4_Tfh")
cellphone_db_plot("data/HLN_cpdb", senders = cells_use, receivers = cells_use, cluster_interactions = TRUE, 
                  filter_secreted  = FALSE, remove_auto = FALSE, scaling = 'min_max', 
                  filter_integrins = FALSE,  
                  annotations_drop = "guidetopharmacology.org")


