#####----------------------------------------------------------------------#####
##### ADD CELL ANNOTATIONS TO SEURAT OBJECTS
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "cell_type_annotation.R"))

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

# all.PROJECTS <- c("SBharadwaj_20240318_Sample_2_5",
#                   "SBharadwaj_20240318_Sample_1_4", 
#                   "SBharadwaj_20240318_Sample_3_6",
#                   "SBharadwaj_20240318_Sample_1_4_7_8",
#                   "SBharadwaj_20240318_Sample_4_8",
#                   "SBharadwaj_20240318_Sample_1_7",
#                   "SBharadwaj_20240318_Sample_7_8",
#                   "SBharadwaj_20240318_Sample_2_3_5_6",
#                   "SBharadwaj_20240318_Sample_1_4_7_8_2_5")

all.PROJECTS <- c("SBharadwaj_20240318_Sample_1_4_7_8", 
                  "SBharadwaj_20240318_Sample_1_4_7_8_2_5")

for (PROJECT in all.PROJECTS){
  path.to.main.input <- file.path(outdir, PROJECT)
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  path.to.01.output <- file.path(path.to.main.output, "01_output")
  
  s.obj <- readRDS(file.path(path.to.01.output, sprintf("Project_%s.rds", PROJECT)))
  
  ##### add Cell annotation
  cell.annotations <- all.annotations[[PROJECT]]
  
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(cell.annotation = cell.annotations[[cca.cluster.0.5]]) %>%
    column_to_rownames("barcode")
  
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, col.name = "cell.annotation", metadata = meta.data$cell.annotation)
  
  ##### add CONDITIONS
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>% rowwise() %>%
    mutate(condition = paste(str_split(name, "_")[[1]][1:3], collapse = "_")) %>%
    mutate(sampletype = paste(str_split(name, "_")[[1]][2:3], collapse = "_")) %>%
    column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, col.name = "condition", metadata = meta.data$condition)
  s.obj <- AddMetaData(object = s.obj, col.name = "sampletype", metadata = meta.data$sampletype)
  
  saveRDS(s.obj, file.path(path.to.01.output, sprintf("Project_%s.addedCellAnnotation.rds", PROJECT)))
}
