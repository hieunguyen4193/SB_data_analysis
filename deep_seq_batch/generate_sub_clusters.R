# gc()
# rm(list = ls())
# use the sub-clustering indices from
# /home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/sub_clustering_indices.R

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"
scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))
source(file.path(path.to.project.src, "sub_clustering_indices.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering_SeuratV5.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
# outdir <- "/home/hieunguyen/CRC1382/outdir"
outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"
# PROJECT <- "SBharadwaj_20240318_Sample_3_6"
# PROJECT <- "SBharadwaj_20240318_Sample_1_4_7_8_2_5"
PROJECT <- "SBharadwaj_20240318_Sample_1_4_7_8"
  
sub.cluster.idx <- "v0.1"

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(path.to.main.input, "data_analysis")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.save.subclusters <- file.path(path.to.main.output, "sub_clusters", sub.cluster.idx)
dir.create(path.to.save.subclusters, showWarnings = FALSE, recursive = TRUE)

s.obj.raw <- readRDS(file.path(path.to.01.output,sprintf("Project_%s.rds", PROJECT)))
s.obj <- subset(s.obj.raw, cca.cluster.0.5 %in% sub_clusters[[PROJECT]][[sub.cluster.idx]])

##### re-integration

num.PCA <- 30
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
num.dim.integration <- 30
num.dim.cluster <- 30
cluster.resolution <- 0.5
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt")

DefaultAssay(s.obj) <- "RNA"
s.obj <- JoinLayers(s.obj)
s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                     save.RDS.s8 = TRUE,
                                                     path.to.output = path.to.save.subclusters,
                                                     use.sctransform = TRUE,
                                                     num.PCA = num.PCA,
                                                     num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                     num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                     cluster.resolution = cluster.resolution,
                                                     vars.to.regress = vars.to.regress,
                                                     PROJECT = PROJECT)