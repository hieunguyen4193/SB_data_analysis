# gc()
# rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

options(future.globals.maxSize = 8000 * 1024^2)

# __________GEX DATA ANALYSIS PIPELINE__________
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
path2src <- file.path(path.to.pipeline.src, "processes_src")

source(file.path(path2src, "import_libraries.R"))
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

path.to.storage <- "/media/hieunguyen/HD0/storage"

#####----------------------------------------------------------------------#####
##### HELPER FUNCTIONS
#####----------------------------------------------------------------------#####
run_pipeline_for_some_combinations <- function(stage_lst, PROJECT){
  
  path2input <- file.path(path.to.storage, "SBharadwaj_20240318_merged_dupl_samples")
  
  MINCELLS  <- 5
  MINGENES  <- 50
  
  save.RDS <- list(s1 = TRUE,
                   s2 = TRUE,
                   s3 = TRUE,
                   s4 = TRUE,
                   s5 = TRUE,
                   s6 = TRUE,
                   s7 = FALSE,
                   s8 = TRUE,
                   s8a = TRUE,
                   s9 = FALSE)
  
  sw <- list(s1 = "on",
             s2 = "on",
             s3 = "on",
             s4 = "on",
             s5 = "on",
             s6 = "on",
             s7 = "off",
             s8 = "on",
             s8a = "off",
             s9 = "off")
  
  rerun <- list(s1 = FALSE, 
                s2 = FALSE,
                s3 = FALSE,
                s4 = FALSE,
                s5 = FALSE,
                s6 = FALSE,
                s7 = FALSE,
                s8 = FALSE,
                s8a = FALSE,
                s9 = FALSE)
  
  
  filter.thresholds <- list(nFeatureRNAfloor = NULL,
                            nFeatureRNAceiling = NULL,
                            nCountRNAfloor = NULL, 
                            nCountRNAceiling = NULL,
                            pct_mitofloor = NULL, 
                            pct_mitoceiling = 10,
                            pct_ribofloor = NULL, 
                            pct_riboceiling = NULL,
                            ambientRNA_thres = 0.5)
  
  remove_doublet <- FALSE
  path.to.10X.doublet.estimation <- file.path(path.to.storage, "DoubletEstimation10X.csv")
  
  num.PCA <- 30
  num.PC.used.in.UMAP <- 30
  num.PC.used.in.Clustering <- 30
  num.dim.integration <- 30
  num.dim.cluster <- 30
  cluster.resolution <- 0.5
  use.sctransform <- TRUE
  vars.to.regress <- c("percent.mt")
  outdir <- "/home/hieunguyen/CRC1382/outdir"
  path.to.output <- file.path(outdir, "SBharadwaj_20240318_merged_dupl_samples", PROJECT)
  dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
  
  filtered.barcodes <- NULL
  
  s.obj <- run_pipeline_GEX(path2src = path2src,
                            path2input = path2input,
                            path.to.logfile.dir = file.path(path.to.output, "logs"),
                            stage_lst = stage_lst,
                            path.to.10X.doublet.estimation = path.to.10X.doublet.estimation,
                            MINCELLS = MINCELLS,
                            MINGENES = MINGENES,
                            PROJECT = PROJECT,
                            remove_doublet = remove_doublet,
                            save.RDS = save.RDS,
                            path.to.output = path.to.output,
                            rerun = rerun, 
                            DE.test = "wilcox",
                            num.PCA = num.PCA,
                            num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                            num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                            use.sctransform = use.sctransform,
                            filtered.barcodes = filtered.barcodes,
                            filter.thresholds = filter.thresholds,
                            input.method = "normal",
                            path.to.anno.contigs = path.to.anno.contigs,
                            path.to.count.clonaltype = path.to.count.clonaltype,
                            cluster.resolution = cluster.resolution,
                            num.dim.integration = num.dim.integration,
                            k.filter = NA,
                            with.TSNE = TRUE,
                            with.VDJ = FALSE,
                            sw = sw,
                            vars.to.regress = vars.to.regress)
  
  #### ALWAYS REMEMBER TO SAVE SESSIONINFO !!!!!!
  writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))
}


#####----------------------------------------------------------------------#####
##### PREPARATION
#####----------------------------------------------------------------------#####
input.stage_lst <- hash()

input.stage_lst[["Gut_CD45"]] <- c(
  ctrl_gut_CD45 =   "ctrl_gut_CD45",
  treated_gut_CD45 = "treated_gut_CD45"
)

input.stage_lst[["Gut_myeloid"]] <- c(
  ctrl_gut_myeloid =   "ctrl_gut_myeloid",
  treated_gut_myeloid = "treated_gut_myeloid"
)

input.stage_lst[["Liver_myeloid"]] <- c(
  ctrl_liver_myeloid =   "ctrl_liver_myeloid",
  treated_liver_myeloid = "treated_liver_myeloid"
)

for (input.name in names(input.stage_lst)){
  PROJECT <- input.name
  print(sprintf("Working on project %s", PROJECT))
  run_pipeline_for_some_combinations(stage_lst = input.stage_lst[[input.name]], PROJECT =  PROJECT)
}