library(Seurat)
library(stringr)

preprocess_dataset <- function(filepath, experiment, sample_name, umi_cutoff, savepath) {
  
  print(paste('Preprocessing:', filepath, sep=' '))
  
  out.data <- Read10X_h5(file = filepath)
  out <- CreateSeuratObject(counts = out.data, project = experiment)
  rm(out.data)
  out[["percent.mt"]] <- PercentageFeatureSet(out, pattern = "^mt-|^MT-")
  expr1 <- FetchData(object = out, vars = "nCount_RNA")
  expr2 <- FetchData(object = out, vars = "percent.mt")
  out <- out[, which(x = expr1 >= umi_cutoff & expr2 <= 15)]
  out$sample <- sample_name
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  s.genes <- str_to_title(s.genes)
  g2m.genes <- str_to_title(g2m.genes)
   
  out <- CellCycleScoring(object = out, s.features = s.genes, g2m.features = g2m.genes)
  out <- SCTransform(object = out, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE, return.only.var.genes = FALSE)

  dir.create(savepath, showWarnings = FALSE)
  saveRDS(out, file=paste(savepath, 'Processed.rds', sep='/'))
  
  return(out)
  
}
