# gc()
# rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"
scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))
source(file.path(path.to.project.src, "sub_clustering_indices.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering_SeuratV5.R"))
source(file.path(path.to.project.src, "sub_clustering_indices.continued_from_06.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
num.PCA <- 30
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
num.dim.integration <- 30
num.dim.cluster <- 30
cluster.resolution <- 0.5
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt")
my_random_seed <- 42
outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"
PROJECT <- "SBharadwaj_20240318_Sample_1_4_7_8_2_5"
sub.cluster.idx <- "v0.1"
# cont.sub.cluster.idx <- "Monocyte_Macrophages"

for (cont.sub.cluster.idx in names(sub_clusters[[PROJECT]][[sub.cluster.idx]])){
  print(sprintf("Working on sub cluster: %s", cont.sub.cluster.idx))
  path.to.main.input <- file.path(outdir, PROJECT)
  path.to.main.output <- file.path(path.to.main.input, "data_analysis")
  path.to.06.output <- file.path(path.to.main.output, "06_output", sub.cluster.idx)
  
  path.to.07.output <- file.path(path.to.main.output, "07_output", sub.cluster.idx, cont.sub.cluster.idx)
  dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)
  
  s.obj.raw <- readRDS(file.path(path.to.06.output, sprintf("Project_%s_%s.rds", PROJECT, sub.cluster.idx)))
  s.obj <- subset(s.obj.raw, cca.cluster.0.5 %in% sub_clusters[[PROJECT]][[sub.cluster.idx]][[cont.sub.cluster.idx]])
  
  if (file.exists(file.path(path.to.07.output, "s8_output", sprintf("%s_%s.noIntegration.rds", PROJECT, cont.sub.cluster.idx))) == FALSE){
    s.obj.no.integrated <- s.obj
    s.obj.no.integrated <- DietSeurat(s.obj.no.integrated)
    DefaultAssay(s.obj.no.integrated) <- "RNA"
    s.obj.no.integrated <- SCTransform(s.obj.no.integrated, vars.to.regress = vars.to.regress, verbose = FALSE)
    pca_reduction_name <- "SCT_PCA"
    umap_reduction_name <- "SCT_UMAP"
    
    s.obj.no.integrated <- RunPCA(s.obj.no.integrated, npcs = num.PCA, verbose = FALSE, reduction.name=pca_reduction_name)
    s.obj.no.integrated <- RunUMAP(s.obj.no.integrated, reduction = pca_reduction_name, 
                                   dims = 1:num.PC.used.in.UMAP, reduction.name=umap_reduction_name,
                                   seed.use = my_random_seed, umap.method = "uwot")
    # clustering 
    s.obj.no.integrated <- FindNeighbors(s.obj.no.integrated, reduction = pca_reduction_name, dims = 1:num.PC.used.in.Clustering)
    s.obj.no.integrated <- FindClusters(s.obj.no.integrated, resolution = cluster.resolution, random.seed = 0)
    
    dir.create(file.path(path.to.07.output, "s8_output"), showWarnings = FALSE, recursive = TRUE)
    
    saveRDS(s.obj.no.integrated, file.path(path.to.07.output, "s8_output", sprintf("%s_%s.noIntegration.rds", PROJECT, cont.sub.cluster.idx)))
  } 
  
  ##### re-integration
  if (file.exists(file.path(path.to.07.output, "s8_output", sprintf("%s_%s.output.s8.rds", PROJECT, cont.sub.cluster.idx))) == FALSE){
    DefaultAssay(s.obj) <- "RNA"
    s.obj <- JoinLayers(s.obj)
    s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj,
                                                         save.RDS.s8 = TRUE,
                                                         path.to.output = path.to.07.output,
                                                         use.sctransform = TRUE,
                                                         num.PCA = num.PCA,
                                                         num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                         num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                         cluster.resolution = cluster.resolution,
                                                         vars.to.regress = vars.to.regress,
                                                         integration.methods = "CCA",
                                                         PROJECT = sprintf("%s_%s", PROJECT, cont.sub.cluster.idx))
  } else {
    print("integrated data exists")
    # s.obj.integrated <- readRDS(file.path(path.to.07.output, "s8_output", sprintf("%s_%s.output.s8.rds", PROJECT, cont.sub.cluster.idx)))
  }
}
  

 