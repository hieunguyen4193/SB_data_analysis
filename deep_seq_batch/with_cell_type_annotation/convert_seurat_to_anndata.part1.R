##### clean up #####
gc()
rm(list = ls())

library(scales)
my_random_seed <- 42
set.seed(my_random_seed)

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)
if ("monocle" %in% installed.packages() == FALSE){
  BiocManager::install("monocle", update = FALSE)
}
library(monocle)

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/with_cell_type_annotation"
samplesheet <- read.csv(file.path(path.to.main.src, "SampleSheet_for_DGE_CellChat.CellAnnotated.csv"))
samplesheet <- samplesheet %>% rowwise() %>%
  mutate(full.dataset.name = sprintf("%s_%s", PROJECT, dataset_name))
samplesheet <- samplesheet[!duplicated(samplesheet$full.dataset.name), ]

for (full.name in unique(samplesheet$full.dataset.name)){
  print(sprintf("Working on dataset %s", full.name))
  PROJECT <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$PROJECT
  dataset_name <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$dataset_name
  path.to.input.s.obj <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$path
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis_with_cell_annotations")
  
  path.to.seurat2anndata <- file.path(path.to.main.output, "09_output", "seurat2anndata")
  dir.create(path.to.seurat2anndata, showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.seurat2anndata, sprintf("%s.csv", full.name))) == FALSE){
    s.obj <- readRDS(path.to.input.s.obj)
    
    s.obj$barcode <- colnames(s.obj)
    
    s.obj$UMAP_1 <- s.obj@reductions$integrated.cca@cell.embeddings[,1]
    s.obj$UMAP_2 <- s.obj@reductions$integrated.cca@cell.embeddings[,2]
    
    write.csv(s.obj@reductions$integrated.cca@cell.embeddings, 
              file=file.path(path.to.seurat2anndata, sprintf('pca_%s.csv', full.name)), 
              quote=F, 
              row.names=F)
    
    write.csv(s.obj@meta.data, 
              file = file.path(path.to.seurat2anndata, sprintf('metadata_%s.csv', full.name)), 
              quote = F, 
              row.names = F)
    
    # write expression counts matrix
    library(Matrix)
    counts_matrix <- GetAssayData(s.obj, assay='SCT', slot='data')
    writeMM(counts_matrix, 
            file = file.path(path.to.seurat2anndata, sprintf('counts_%s.mtx', full.name)))
    
    # write gene names
    write.table( data.frame(
      'gene' = rownames(counts_matrix)),
      file = file.path(path.to.seurat2anndata, sprintf('gene_names_%s.csv', full.name)),
      quote = F,
      row.names = F,
      col.names = F)
    
    coldf <- data.frame(cluster = unique(s.obj$cell.annotation),
                        color = hue_pal()(length(unique(s.obj$cell.annotation))))
    write.csv(coldf, file.path(path.to.seurat2anndata, sprintf('colordf_%s.csv', full.name)))
    write.csv(data.frame(status = c("finished generating conversion")), 
              file.path(path.to.seurat2anndata, sprintf("%s.csv", full.name)))
  }

}
