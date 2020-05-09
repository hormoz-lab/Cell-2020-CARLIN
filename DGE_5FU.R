library(Seurat)
library(cowplot)

DGE_5FU <- function(analysis_path, results_path) {
  
  integrated <- readRDS(paste(analysis_path, '5FU/V3Combined/Transcriptome/Seurat/Combined.rds', sep='/'))
  
  HSC_type_labels    = read.table(file=paste(results_path, '5FU/Seurat/DGE_HSC_Type_Labels.csv', sep='/'), sep=",", header=FALSE, stringsAsFactors=FALSE)
  HSC_cluster_labels = read.table(file=paste(results_path, '5FU/Seurat/DGE_HSC_Cluster_Labels.csv', sep='/'), sep=",", header=FALSE, stringsAsFactors=FALSE)
  
  meta = integrated@meta.data
  meta["HSC_type"] <- HSC_type_labels
  meta["HSC_cluster"] <- HSC_cluster_labels
  meta <- subset(meta, select=c("HSC_type", "HSC_cluster"))
  integrated <- AddMetaData(integrated, meta)
  
  DefaultAssay(integrated) <- 'RNA'
  
  Idents(object = integrated) <- 'HSC_type'
  DGE_HSC_type <- FindMarkers(object=integrated, logfc.threshold = 0.2, ident.1 = 'Parent', ident.2 = 'Childless')
  write.table(file=paste(results_path, '5FU/Seurat/DGE_HSC_Type.csv', sep='/'), x = DGE_HSC_type, sep=",")
  
  Idents(object = integrated) <- 'HSC_cluster'
  DGE_HSC_cluster <- FindMarkers(object=integrated, logfc.threshold = 0.2, ident.1 = 'Parent', ident.2 = 'Childless')
  write.table(file=paste(results_path, '5FU/Seurat/DGE_HSC_Cluster.csv', sep='/'), x = DGE_HSC_cluster, sep=",")
  
  VlnPlot(object=integrated, features=c("Plac8", "Cd34", "Mllt3"), idents=c('Parent', 'Childless'), ncol=6, pt.size=0.1)
  VlnPlot(object=integrated, features=c("Cdk6", "Serpinb1a", "Sox4", "Mpo"), idents=c('Parent', 'Childless'), ncol=6, pt.size=0.1)
  VlnPlot(object=integrated, features=c("Cd74", "Apoe", "Pdzk1ip1", "Hacd4"), idents=c('Parent', 'Childless'), ncol=6, pt.size=0.1)
  
  DGE_HSC_cluster <- DGE_HSC_cluster[DGE_HSC_cluster$p_val_adj<0.05,]
  DGE_HSC_cluster <- DGE_HSC_cluster[order(-DGE_HSC_cluster$avg_logFC),]
  DoHeatmap(integrated, cells=WhichCells(integrated, idents=c('Parent', 'Childless')), features=rownames(DGE_HSC_cluster[c(1:12, nrow(DGE_HSC_cluster)+(-11:0)),]))

  proliferation_genes = scan("proliferation_genes.txt", what="", sep="\n")
  integrated <- AddModuleScore(integrated, features=list(proliferation_genes), name="Proliferation.Score")
  
  DGE_HSC_prolif_labels <- data.frame(c(rep("Other", length(integrated$Proliferation.Score1))), stringsAsFactors = FALSE)
  q <- quantile(integrated$Proliferation.Score1[integrated$HSC_cluster!="Other"])
  DGE_HSC_prolif_labels[integrated$HSC_cluster!="Other" & integrated$Proliferation.Score1 <= q[2],] <- "Low"
  DGE_HSC_prolif_labels[integrated$HSC_cluster!="Other" & integrated$Proliferation.Score1 > q[2] & integrated$Proliferation.Score1 < q[4],] <- "Mid"
  DGE_HSC_prolif_labels[integrated$HSC_cluster!="Other" & integrated$Proliferation.Score1 >= q[4],] <- "High"
  
  meta <- integrated@meta.data
  meta["HSC_prolif"] <- DGE_HSC_prolif_labels
  meta <- subset(meta, select=c("HSC_prolif"))
  integrated <- AddMetaData(integrated, meta)
  
  Idents(object=integrated) <- 'HSC_prolif'
  DGE_HSC_prolif <- FindMarkers(object=integrated, logfc.threshold = 0.2, ident.1='High', ident.2='Low')
  
  write.table(file=paste(results_path, '5FU/Seurat/DGE_HSC_Prolif_Labels.csv', sep='/'), x = DGE_HSC_prolif_labels, sep=",", row.names=FALSE, col.names=FALSE)
  write.table(file=paste(results_path, '5FU/Seurat/DGE_HSC_Prolif.csv', sep='/'), x = DGE_HSC_prolif, sep=",")
 
}