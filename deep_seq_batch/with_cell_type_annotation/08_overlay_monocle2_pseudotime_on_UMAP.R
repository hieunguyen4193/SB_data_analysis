##### clean up #####
gc()
rm(list = ls())

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
  # full.name <- unique(samplesheet$full.dataset.name)[[1]]
  
  print(sprintf("Working on dataset %s", full.name))
  PROJECT <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$PROJECT
  dataset_name <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$dataset_name
  path.to.input.s.obj <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$path
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  path.to.monocle2.input <- file.path(path.to.main.output, "08_output", "monocle2_input")
  path.to.08.output <- file.path(path.to.main.output, "08_output", "monocle_output", sprintf("%s.%s.monocle2", PROJECT, dataset_name))
  dir.create(path.to.monocle2.input, showWarnings = FALSE, recursive = TRUE)
  
  monocledf <- read.csv(file.path(path.to.08.output, "monocledf.csv"))
  monocle.reversedf <- read.csv(file.path(path.to.08.output, "monocledf.rev.csv")) %>%
    subset(select = -c(X))
  colnames(monocle.reversedf) <- c("barcode", "rev.state", "rev.pseudotime")
  s.obj <- readRDS(path.to.input.s.obj)
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
  meta.data <- merge(meta.data, monocledf, by.x = "barcode", by.y = "barcode")
  meta.data <- merge(meta.data, monocle.reversedf, by.x = "barcode", by.y = "barcode") 
  meta.data <- meta.data %>% column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data),]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$pseudotime, col.name = "pseudotime")
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$rev.pseudotime, col.name = "rev.pseudotime")
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$state, col.name = "state")
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$rev.state, col.name = "rev.state")
  
  Idents(s.obj) <- "cell_annotation"
  
  umap.pseudotime <- FeaturePlot(object = s.obj, reduction = "cca_UMAP", features = c("pseudotime"), label = TRUE)
  umap.rev.pseudotime <- FeaturePlot(object = s.obj, reduction = "cca_UMAP", features = c("rev.pseudotime"), label = TRUE)
  
  umap.state <- DimPlot(object = s.obj, reduction = "cca_UMAP", group.by = "state")
  umap.rev.state <- DimPlot(object = s.obj, reduction = "cca_UMAP", group.by = "rev.state")
  
  ggsave(plot = umap.pseudotime, filename = sprintf("%s.UMAP_pseudotime.svg", full.name), path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
  ggsave(plot = umap.rev.pseudotime, filename = sprintf("%s.UMAP_rev_pseudotime.svg", full.name), path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
  ggsave(plot = umap.state, filename = sprintf("%s.UMAP_state.svg", full.name), path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
  ggsave(plot = umap.rev.state, filename = sprintf("%s.UMAP_rev_state.svg", full.name), path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
}