#####----------------------------------------------------------------------#####
##### 01
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))
library(scales)
outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20250102"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/release"
samplesheet <- read.csv(file.path(path.to.main.src, "SampleSheet_all_seurat_objects.csv"))
samplesheet <- samplesheet %>% rowwise() %>%
  mutate(dataset_name = ifelse(reIntegration == "yes", sprintf("%s_reIntegration", dataset_name), dataset_name))

for (row_i in seq(1, nrow(samplesheet))){
  PROJECT <- samplesheet[row_i, ]$PROJECT
  dataset_name <- samplesheet[row_i, ]$dataset_name
  path.to.s.obj <- samplesheet[row_i, ]$path
  re.integration <- samplesheet[row_i, ]$reIntegration
  
  if (re.integration %in% c("yes", "")){
    reduction.name <- "cca_UMAP"
    cluster.name <- "cca.cluster.0.5"
  } else {
    reduction.name <- "SCT_UMAP"
    cluster.name <- "seurat_clusters"
  }
    
  print(sprintf("Working on PROJECT %s, dataset name %s", PROJECT, dataset_name))
  
  path.to.s.obj <- str_replace(path.to.s.obj, ".rds", ".addedInfo.rds")
  
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  path.to.08.output <- file.path(path.to.main.output, "08_output", dataset_name)
  dir.create(path.to.08.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.seurat2anndata <- file.path(path.to.08.output, "seurat2anndata")
  dir.create(path.to.seurat2anndata, showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.seurat2anndata, sprintf('colordf_%s.csv', sprintf("%s_%s", PROJECT, dataset_name)))) == FALSE){
    s.obj <- readRDS(path.to.s.obj)
    
    s.obj$barcode <- colnames(s.obj)
    
    s.obj$UMAP_1 <- s.obj@reductions[[reduction.name]]@cell.embeddings[,1]
    s.obj$UMAP_2 <- s.obj@reductions[[reduction.name]]@cell.embeddings[,2]
    
    write.csv(s.obj@reductions[[reduction.name]]@cell.embeddings, 
              file=file.path(path.to.seurat2anndata, 
                             sprintf('pca_%s.csv', sprintf("%s_%s", PROJECT, dataset_name))), 
              quote=F, 
              row.names=F)
    
    write.csv(s.obj@meta.data, file=file.path(path.to.seurat2anndata, sprintf('metadata_%s.csv', sprintf("%s_%s", PROJECT, dataset_name))), quote=F, row.names=F)
    
    # write expression counts matrix
    library(Matrix)
    counts_matrix <- GetAssayData(s.obj, assay='SCT', slot='data')
    writeMM(counts_matrix, file=file.path(path.to.seurat2anndata, sprintf('counts_%s.mtx', sprintf("%s_%s", PROJECT, dataset_name))))
    
    # write gene names
    write.table( data.frame('gene'=rownames(counts_matrix)),file=file.path(path.to.seurat2anndata, sprintf('gene_names_%s.csv', sprintf("%s_%s", PROJECT, dataset_name))),
                 quote=F,row.names=F,col.names=F)
    
    coldf <- data.frame(cluster = unique(s.obj@meta.data[[cluster.name]]),
                        color = hue_pal()(length(unique(s.obj@meta.data[[cluster.name]]))))
    write.csv(coldf, file.path(path.to.seurat2anndata, sprintf('colordf_%s.csv', sprintf("%s_%s", PROJECT, dataset_name))))
  } else {
    print(sprintf("File %s exists", file.path(path.to.seurat2anndata, sprintf('colordf_%s.csv', sprintf("%s_%s", PROJECT, dataset_name)))))
  }
}