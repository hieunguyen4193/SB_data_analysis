##### clean up #####
# gc()
# rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
outdir <- "/home/hieunguyen/CRC1382/outdir"

# all.PROJECTS <- c("SBharadwaj_20240318_Sample_1_4",
#                   "SBharadwaj_20240318_Sample_2_5",
#                   "SBharadwaj_20240318_Sample_1_4_7_8",
#                   "SBharadwaj_20240318_Sample_3_6",
#                   "SBharadwaj_20240318_Sample_1_7",
#                   "SBharadwaj_20240318_Sample_4_8",
#                   "SBharadwaj_20240318_Sample_2_3_5_6",
#                   "SBharadwaj_20240318_Sample_7_8")

all.PROJECTS <- c("SBharadwaj_20240318_all_samples")

for (PROJECT in all.PROJECTS){
  
  path.to.main.input <- file.path(outdir, "SeuratV5", PROJECT)
  s.obj.raw <- readRDS(file.path(path.to.main.input, "s1_output", sprintf("%s.output.s1.rds", PROJECT) ))
  s.obj <- readRDS(file.path(path.to.main.input, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
  s.obj <- subset(s.obj, Doublet_classifications == "Singlet")
  
  path.to.main.output <- file.path(path.to.main.input, "data_analysis")
  dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.01.output <- file.path(path.to.main.output, "01_output")
  path.to.save.cloupe <- file.path(path.to.01.output, "cloupe_files")
  dir.create(path.to.save.cloupe, showWarnings = FALSE, recursive = TRUE)

  print(sprintf("Working on file %s", PROJECT))
  s.obj <- readRDS(file.path(path.to.01.output,sprintf("Project_%s.rds", PROJECT)))
  ##### CONVERT seurat object to cloupe file
  if ("loupeR" %in% installed.packages() == FALSE){
    install.packages("hdf5r")
    install.packages("/media/hieunguyen/HD0/storage/offline_pkgs/loupeR_Linux.tar.gz", repos = NULL, type = "source")
  }

  library(loupeR)

  loupeR::setup()
  create_loupe_from_seurat(
      s.obj,
      output_dir = file.path(path.to.save.cloupe),
      output_name = sprintf("PROJECT_%s_cloupe_converted_from_seurat", PROJECT),
      dedup_clusters = FALSE,
      executable_path = NULL,
      force = FALSE
      )
}