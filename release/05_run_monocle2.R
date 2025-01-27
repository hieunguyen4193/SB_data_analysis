##### clean up #####
gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

# to solve unable to access index for repository https://mran.microsoft.com/snapshot/2020-07-16/src/contrib
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)})

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)
if ("monocle" %in% installed.packages() == FALSE){
  BiocManager::install("monocle", update = FALSE)
}
library(monocle)
library(dplyr)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"
source(file.path(path.to.main.src, "monocle2_helper_functions.R"))

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20250102"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/release"
samplesheet <- read.csv(file.path(path.to.main.src, "SampleSheet_all_seurat_objects.csv"))
samplesheet <- samplesheet %>% rowwise() %>%
  mutate(dataset_name = ifelse(reIntegration == "yes", sprintf("%s_reIntegration", dataset_name), dataset_name))

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

row_i <- 1
# for (row_i in seq(1, nrow(samplesheet))){
  PROJECT <- samplesheet[row_i, ]$PROJECT
  dataset_name <- samplesheet[row_i, ]$dataset_name
  re.integration <- samplesheet[row_i, ]$reIntegration
  if (dataset_name == "full"){
    if (re.integration %in% c("yes", "")){
      to.run.clusters <- c("cca.cluster.0.5", "cell.annotation")
      reduction.name <- "cca_UMAP"
    } else {
      to.run.clusters <- c("seurat_clusters", "cell.annotation")
      reduction.name <- "SCT_UMAP"
    }
  } else {
    if (re.integration %in% c("yes", "")){
      to.run.clusters <- c("cca.cluster.0.5")      
      reduction.name <- "cca_UMAP"
    } else {
      to.run.clusters <- c("seurat_clusters")
      reduction.name <- "SCT_UMAP"
    }
  }
  
  if (PROJECT == "SBharadwaj_20240318_Sample_3_6"){
    to.run.clusters <- setdiff(to.run.clusters, "cell.annotation")
  }
  
  all.cases <- sample.list[[PROJECT]]
  for (cluster.name in to.run.clusters){
    path.to.s.obj <- samplesheet[row_i, ]$path
    path.to.s.obj <- str_replace(path.to.s.obj, ".rds", ".addedInfo.rds")
    
    path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
    path.to.05.output <- file.path(path.to.main.output, "05_output")
    
    path.to.save.output <- file.path(path.to.05.output, dataset_name, cluster.name, "monocle2_input")
    dir.create(file.path(path.to.05.output, "monocle2_output"), showWarnings = FALSE, recursive = TRUE)
    
    monocle.obj <- readRDS(file.path(path.to.save.output, "monocle2_obj.rds"))
    monocle.obj <- run_monocle2_from_presave_obj(monocle.obj, 
                                                 file.path(path.to.05.output, 
                                                           "monocle2_output"))
    
    ##### plot cell trajectory, color by seurat clusters
    p <- plot_cell_trajectory(monocle.obj, color_by = "cca.cluster.0.5")
    ggsave(plot = p, 
           filename = sprintf("cell_trajectory_%s.seurat_clsuters.svg", full.name), 
           path = path.to.08.output, 
           device = "svg", 
           dpi = 300, 
           width = 14, 
           height = 10)
    
    ##### plot cell trajectory, color by monocle2 states
    p <- plot_cell_trajectory(monocle.obj, color_by = "State")
    ggsave(plot = p, 
           filename = sprintf("cell_trajectory_%s.State.svg", full.name), 
           path = path.to.08.output, 
           device = "svg", 
           dpi = 300, 
           width = 14, 
           height = 10)
    
    ##### plot cell trajectory, color by pseudotime
    p <- plot_cell_trajectory(monocle.obj, color_by = "Pseudotime")
    ggsave(plot = p, filename = sprintf("cell_trajectory_%s.pseudotime.svg", full.name), 
           path = path.to.08.output, 
           device = "svg", 
           dpi = 300, 
           width = 14, 
           height = 10)
    
    ##### plot cell trajectory, color by pseudotime
    monocle.obj.reverse <- orderCells(monocle.obj, reverse = TRUE)
    p <- plot_cell_trajectory(monocle.obj.reverse, color_by = "Pseudotime")
    ggsave(plot = p, 
           filename = sprintf("cell_trajectory_%s.rev_Pseudotime.svg", full.name), 
           path = path.to.08.output, 
           device = "svg", 
           dpi = 300, 
           width = 14, 
           height = 10)
    
    ##### save monocle data to csv file
    monocledf <- data.frame(
      barcode = colnames(monocle.obj),
      state = monocle.obj$State,
      pseudotime = monocle.obj$Pseudotime
    )
    monocle.reversedf <- data.frame(
      barcode = colnames(monocle.obj.reverse),
      state = monocle.obj.reverse$State,
      pseudotime = monocle.obj.reverse$Pseudotime
    )
    write.csv(monocledf, file.path(path.to.08.output, "monocledf.csv"))
    write.csv(monocle.reversedf, file.path(path.to.08.output, "monocledf.rev.csv"))
  }
# }