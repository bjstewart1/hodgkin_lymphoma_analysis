#microarray profiling of Reed Sternberg cells.
#stiedl data for HRSC vs germinal centre ####
gse <- GEOquery::getGEO("GSE39133")
gse <- gse$GSE39133_series_matrix.txt.gz
library(Biobase)
fData(gse)[, "Gene Symbol"] <- gsub(" ", "", unlist(lapply(strsplit(fData(gse)[, "Gene Symbol"], "///"), function(x){x[1]})))
pData(gse)
#trim out genes without symbols
gse <- gse[!is.na(fData(gse)[, "Gene Symbol"]), ]
#remove duplicates
gse <- gse[!duplicated(fData(gse)[, "Gene Symbol"]), ]
saveRDS(gse, "data/GSE39133.RDS")
#do a limma analysis for differential expression ####
library(limma)
group <- pData(gse)$characteristics_ch1.1
design <- model.matrix(~0+group)
colnames(design) <- c("HRS", "GC") #we compare HRS cells vs germinal centres (GC)
contr.matrix <- makeContrasts(
  HRSvsGC = HRS-GC,
  levels = colnames(design))
#do this with limma 
fit <- lmFit(exprs(gse), design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit)
summary(decideTests(efit))
tfit <- treat(fit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
top <- topTreat(tfit, coef=1, n=Inf)
top$gene <- fData(gse)[rownames(top), "Gene Symbol"]
top <- top[order(top$logFC), ]
library(ggplot2)
top$significant <- FALSE
top$significant[top$logFC > 3 & top$adj.P.Val < 0.01] <- "HRS"
top$significant[top$logFC < -3 & top$adj.P.Val < 0.01] <- "GC"
top <- top[order(top$significant,decreasing = FALSE), ] #order of plot.
ggplot(top, aes(x = logFC, y = -log(P.Value), color = significant, size = significant)) + ggrastr::geom_point_rast(pch = 19) + 
  theme_classic() + theme(axis.text = element_text(color = "black")) + scale_color_manual(values = c("FALSE" = "grey50", "HRS" = "black", "GC" = "black")) + 
  scale_size_manual(values = c("FALSE" = 0.5, "HRS" = 1, "GC" = 1))

saveRDS(top, "data/HRS_DE.RDS")

#get the HRS signature ####
HRS_signature <- top[top$significant %in% "HRS", ]
saveRDS(HRS_signature, "data/HRSC_signature.RDS")

#wrange the expression matrix and differential expression into shape for progeny and dorothea ####
library(dplyr)
expression_matrix <- data.frame(exprs(gse), gene = fData(gse)[, "Gene Symbol"], row.names = NULL) %>% dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>%  tibble::column_to_rownames(var = "gene") %>% as.matrix()
top <- data.frame(top, ID = top$gene, row.names = NULL)
top_mat <- top %>% dplyr::select(ID, t) %>% dplyr::filter(!is.na(t)) %>% tibble::column_to_rownames(var = "ID") %>% as.matrix()

#run dorothea ####
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% dplyr::filter(confidence %in% c("A", "B","C"))
tf_activities_stat <- dorothea::run_viper(top_mat, regulons, options =  list(minsize = 5, eset.filter = FALSE,  cores = 1, verbose = FALSE, nes = TRUE))
#get the top 10 TFs
tf_activities_stat_top10 <- tf_activities_stat %>% as.data.frame() %>% tibble::rownames_to_column(var = "GeneID") %>%dplyr::rename(NES = "t") %>% dplyr::top_n(10, wt = abs(NES)) %>%  dplyr::arrange(NES) %>%  dplyr::mutate(GeneID = factor(GeneID))
#plot these
ggplot(tf_activities_stat_top10,aes(x = reorder(GeneID, NES), y = NES)) +  geom_point(size = 5, aes(fill = NES), pch=21, stat = "identity") + 
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "grey90",midpoint = 0, limits = c(-max(abs(tf_activities_stat_top10$NES)), max(abs(tf_activities_stat_top10$NES)))) + 
  theme_bw() + theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size =10, face= "bold"), 
        axis.text.y = element_text(size =10, color = 'black', face= "italic")) + 
  xlab("Transcription Factors") + coord_flip()

#now look to see what the likely signalling upstream of this is
tf_activities <- tf_activities_stat %>% as.data.frame() %>%  tibble::rownames_to_column(var = "TF") 

#now create a network view of HRS cells and germinal centre cells #####
library(OmnipathR)
library(dplyr)
TFS <- as.character(tf_activities_stat_top10$GeneID) #get the top TF activities
HRS_genes <- HRS_signature$gene
B_cell_genes <- top[top$significant %in% "GC", 'gene']

#now make an edgebundled graph #####
dorothea_levels = c("A", "B", "C")
interactions <- import_all_interactions( organism = 9606, dorothea_levels =dorothea_levels)
interactions <- interactions[interactions$target_genesymbol %in% c(HRS_genes,TFS),  ]
interactions <- interactions[interactions$source_genesymbol %in% c(TFS, HRS_genes), ]
interaction_gr <- interaction_graph(interactions = interactions)

el1 <- data.frame("from" = "root", "to" = c("TF", "HRSC"))
el2 <- data.frame("from" = "TF", "to" = as.character(V(interaction_gr)$name)[V(interaction_gr)$name %in% TFS])
el3 <- data.frame("from" = "HRSC", "to" = as.character(V(interaction_gr)$name)[V(interaction_gr)$name %in% c(HRS_genes) ])

ig_gr <- data.frame(igraph::as_edgelist(interaction_gr))
colnames(ig_gr) <- c("from", "to")
ig_gr$from <- as.character(ig_gr$from)
ig_gr$to <- as.character(ig_gr$to)
master_el <- do.call(rbind, list(el1, el2, el3))
gr <- graph_from_edgelist(as.matrix(master_el))
from = match( ig_gr$from, V(gr)$name)
to = match(ig_gr$to, V(gr)$name)

#draw the graph
ggraph(gr, layout = "dendrogram", circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to =  to), 
                   alpha=0.2, color = 'grey30', tension = 0.5) + geom_node_point(pch =19, aes(filter = leaf)) + 
  geom_node_point(pch = 21, fill = "black", aes(filter = leaf)) + coord_fixed()  + 
  geom_node_label(aes(label = name, filter = leaf),   nudge_x = 0.1 )
