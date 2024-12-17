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
num.PCA <- 30
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
num.dim.integration <- 30
num.dim.cluster <- 30
cluster.resolution <- 0.5
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt")
my_random_seed <- 42

# outdir <- "/home/hieunguyen/CRC1382/outdir"
outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"
# PROJECT <- "SBharadwaj_20240318_Sample_3_6"
# PROJECT <- "SBharadwaj_20240318_Sample_1_4_7_8_2_5"
# PROJECT <- "SBharadwaj_20240318_Sample_1_4_7_8"

for (PROJECT in c("SBharadwaj_20240318_Sample_3_6",
                  "SBharadwaj_20240318_Sample_1_4_7_8_2_5",
                  "SBharadwaj_20240318_Sample_1_4_7_8")){
  sub.cluster.idx <- "v0.1"
  path.to.main.input <- file.path(outdir, PROJECT)
  path.to.main.output <- file.path(path.to.main.input, "data_analysis")
  path.to.01.output <- file.path(path.to.main.output, "01_output")
  path.to.save.subclusters <- file.path(path.to.main.output, "sub_clusters", sub.cluster.idx)
  dir.create(path.to.save.subclusters, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.save.subclusters, "s8_output"), showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.save.subclusters, "s8_output", sprintf("%s.noIntegration.rds", PROJECT))) == FALSE){
    print(sprintf("WORKING ON THE PROJECT %s", PROJECT))
    s.obj.raw <- readRDS(file.path(path.to.01.output,sprintf("Project_%s.rds", PROJECT)))
    s.obj <- subset(s.obj.raw, cca.cluster.0.5 %in% sub_clusters[[PROJECT]][[sub.cluster.idx]])

    print("generate non-integration data ...")
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
    
    saveRDS(s.obj.no.integrated, file.path(path.to.save.subclusters, "s8_output", sprintf("%s.noIntegration.rds", PROJECT)))
  } else {
    print(sprintf("File %s exists", file.path(path.to.save.subclusters, "s8_output", sprintf("%s.noIntegration.rds", PROJECT))))
  }
  
  ##### re-integration
  if (file.exists(file.path(path.to.save.subclusters, "s8_output", sprintf("%s.output.s8.rds", PROJECT))) == FALSE){
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
  }
}
  


