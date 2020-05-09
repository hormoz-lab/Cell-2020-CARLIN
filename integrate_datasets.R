library(Seurat)
library(stringr)

integrate_datasets <- function(dataset_list, marker_genes, savepath) {
  
  print(paste('Integrating:', savepath, sep=' '))
  
  # Following workflow here: 
  # https://github.com/satijalab/seurat/issues/2023
  # https://github.com/satijalab/seurat/issues/2050
  
  integration_features <- SelectIntegrationFeatures(object.list = dataset_list)
  dataset_list <- PrepSCTIntegration(object.list = dataset_list, anchor.features = union(integration_features, marker_genes))
  integration_anchors <- FindIntegrationAnchors(object.list = dataset_list, normalization.method = "SCT", anchor.features = union(integration_features, marker_genes))
  
  dir.create(savepath, showWarnings = FALSE, recursive = TRUE)
  saveRDS(integration_anchors, file=paste(savepath, 'Anchors.rds', sep='/'))
  
  rm(dataset_list, pos=1)
  
  integrated <- IntegrateData(anchorset = integration_anchors, features.to.integrate = union(integration_anchors@anchor.features, marker_genes), normalization.method = "SCT")
  rm(integration_anchors)
  
  DefaultAssay(integrated) <- "integrated"
  integrated <- RunPCA(object = integrated, npcs = 30, verbose = FALSE)
  integrated <- RunUMAP(object = integrated, reduction = "pca", dims = 1:30)
  integrated <- FindNeighbors(object = integrated, reduction = "pca", dims = 1:30)
  integrated <- FindClusters(integrated, algorithm = 2, resolution = 1.5)
 
  DefaultAssay(integrated) <- "RNA"
  integrated <- NormalizeData(object = integrated, normalization.method = "LogNormalize", scale.factor = 10000)
  integrated <- ScaleData(object = integrated, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
  #integrated <- ScaleData(object = integrated)
  
  saveRDS(integrated, file=paste(savepath, 'Combined.rds', sep='/'))

  return(integrated)
  
}
