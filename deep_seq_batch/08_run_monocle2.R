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
outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"
source(file.path(path.to.main.src, "monocle2_helper_functions.R"))

samplesheet <- read.csv(file.path(path.to.main.src, "SampleSheet_for_DGE_and_CellChat.csv"))
samplesheet <- samplesheet %>% rowwise() %>%
  mutate(full.dataset.name = sprintf("%s_%s", PROJECT, dataset_name))
samplesheet <- samplesheet[!duplicated(samplesheet$full.dataset.name), ]

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

for (full.name in unique(samplesheet$full.dataset.name)){
  print(sprintf("Working on dataset %s", full.name))
  PROJECT <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$PROJECT
  
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  dataset_name <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$dataset_name
  path.to.monocle2.input <- file.path(path.to.main.output, "08_output", "monocle2_input")
  path.to.monocle.obj <- file.path(path.to.monocle2.input, 
                                   sprintf("%s.%s.monocle2.rds", 
                                           PROJECT, 
                                           dataset_name))
  path.to.08.output <- file.path(path.to.main.output, "08_output", 
                                 "monocle_output", sprintf("%s.%s.monocle2",
                                                           PROJECT,
                                                           dataset_name))
  dir.create(path.to.08.output, showWarnings = FALSE, recursive =TRUE)
  
  if (file.exists(file.path(path.to.08.output, "monocledf.rev.csv")) == FALSE){
    
    print("reading monocle object saved from the seurat data object")
    monocle.obj <- readRDS(path.to.monocle.obj)
    
    if (file.exists(file.path(path.to.08.output, "monocle_obj.rds")) == FALSE){
      print(sprintf("running monocle2 on dataset %s", full.name))
      monocle.obj <- readRDS(path.to.monocle.obj)
      monocle.obj <- run_monocle2_from_presave_obj(monocle.obj, path.to.08.output)
    } else {
      print("monocle result exists, reading in...")
      monocle.obj <- readRDS(file.path(path.to.08.output, "monocle_obj.rds"))
    }
    
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
  } else {
    print(sprintf("Analysis %s finished", file.path(path.to.08.output)))
  }
}