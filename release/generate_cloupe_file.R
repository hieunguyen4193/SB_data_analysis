#####----------------------------------------------------------------------#####
##### 01
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20250102"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/release"
samplesheet <- read.csv(file.path(path.to.main.src, "SampleSheet_all_seurat_objects.csv"))
samplesheet <- samplesheet %>% rowwise() %>%
  mutate(dataset_name = ifelse(reIntegration == "yes", sprintf("%s_reIntegration", dataset_name), dataset_name))

sample.list <- list(
  SBharadwaj_20240318_Sample_1_4_7_8_2_5 = list(
    name = c("ctrl_gut_CD45_1",
             "treated_gut_CD45_1",
             "treated_gut_CD45_2",
             "ctrl_gut_CD45_2",
             "ctrl_gut_myeloid",
             "treated_gut_myeloid"),
    condition = c(c("ctrl_gut_CD45",
                    "treated_gut_CD45",
                    "ctrl_gut_myeloid",
                    "treated_gut_myeloid")),
    global.condition = c("ctrl", "treated")
  ),
  `SBharadwaj_20240318_Sample_1_4_7_8` = list(
    `condition` = c("treated_gut_CD45", "ctrl_gut_CD45"),
    name = c("ctrl_gut_CD45_1",
             "treated_gut_CD45_1",
             "treated_gut_CD45_2",
             "ctrl_gut_CD45_2")
  ),
  `SBharadwaj_20240318_Sample_3_6` = list(
    `condition` = c("treated_liver_myeloid", "ctrl_liver_myeloid")
  )
)

path.to.save.cloupe.file <- file.path(outdir, "cloupe_files")
dir.create(path.to.save.cloupe.file, showWarnings = FALSE, recursive = TRUE)

#### CONVERT seurat object to cloupe file
if ("loupeR" %in% installed.packages() == FALSE){
  install.packages("hdf5r")
  install.packages("/media/hieunguyen/HD01/storage/offline_pkgs/loupeR_Linux.tar.gz", repos = NULL, type = "source")
}

library(loupeR)

loupeR::setup()

for (row_i in seq(1, nrow(samplesheet))){
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
  path.to.s.obj <- samplesheet[row_i, ]$path
  path.to.s.obj <- str_replace(path.to.s.obj, ".rds", ".addedInfo.rds")
  save.cloupe.name <- sprintf("%s_%s", PROJECT, dataset_name)
  
  print(sprintf("Generating cloupe files for %s", save.cloupe.name))

  s.obj <- readRDS(path.to.s.obj)
  create_loupe_from_seurat(
    s.obj,
    output_dir = file.path(path.to.save.cloupe.file),
    output_name = save.cloupe.name,
    dedup_clusters = FALSE,
    executable_path = NULL,
    force = TRUE
  )
}

