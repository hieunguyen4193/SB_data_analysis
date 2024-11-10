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
  
  s.obj <- readRDS(file.path(path.to.main.input, 
                             "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
  s.obj <- subset(s.obj, Doublet_classifications == "Singlet")
  
  cell.annotations <- all.annotations[[PROJECT]]
  
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(cell.annotation = cell.annotations[[cca.cluster.0.5]]) %>%
    column_to_rownames("barcode")
  
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, col.name = "cell.annotation", metadata = meta.data$cell.annotation)
  saveRDS(s.obj, 
          file.path(path.to.main.input, 
                    "s8_output", sprintf("%s.addedCellAnnotation.output.s8.rds", PROJECT)))
}
