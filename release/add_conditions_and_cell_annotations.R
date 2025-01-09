#####----------------------------------------------------------------------#####
##### ADD CELL ANNOTATIONS TO SEURAT OBJECTS
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/release"

source(file.path(path.to.project.src, "cell_type_annotation.R"))

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))
source(file.path(path.to.project.src, "cell_type_annotation.R"))

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20250102"
path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/release"
samplesheet <- read.csv(file.path(path.to.main.src, "SampleSheet_all_seurat_objects.csv"))

rerun <- TRUE
for (row_i in seq(1, nrow(samplesheet))){
  PROJECT <- samplesheet[row_i, ]$PROJECT
  path.to.s.obj <- samplesheet[row_i, ]$path
  dataset_name <- samplesheet[row_i, ]$dataset_name
  
  if (file.exists(str_replace(path.to.s.obj, ".rds", ".addedInfo.rds")) == FALSE | rerun == TRUE){
    print(sprintf("Adding information to object %s", path.to.s.obj))
    s.obj <- readRDS(path.to.s.obj)
    
    ##### add conditions
    meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>% rowwise() %>%
      mutate(condition = paste(str_split(name, "_")[[1]][1:3], collapse = "_")) %>%
      mutate(sampletype = paste(str_split(name, "_")[[1]][2:3], collapse = "_")) %>% 
      mutate(global.condition = str_split(condition, "_")[[1]][[1]]) %>%
      column_to_rownames("barcode")
    
    meta.data <- meta.data[row.names(s.obj@meta.data), ]
    
    s.obj <- AddMetaData(object = s.obj, col.name = "condition", metadata = meta.data$condition)
    s.obj <- AddMetaData(object = s.obj, col.name = "global.condition", metadata = meta.data$global.condition)
    s.obj <- AddMetaData(object = s.obj, col.name = "sampletype", metadata = meta.data$sampletype)
    
    annotated.PROJECTS <- names(all.annotations)
    
    if (dataset_name == "full" & PROJECT %in% annotated.PROJECTS){
      cell.annotations <- all.annotations[[PROJECT]]
      
      meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
        rowwise() %>%
        mutate(cell.annotation = cell.annotations[[cca.cluster.0.5]]) %>%
        column_to_rownames("barcode")
      
      meta.data <- meta.data[row.names(s.obj@meta.data), ]
      s.obj <- AddMetaData(object = s.obj, col.name = "cell.annotation", metadata = meta.data$cell.annotation)
    }
    saveRDS(s.obj, str_replace(path.to.s.obj, ".rds", ".addedInfo.rds"))
  } else {
    print(sprintf("File %s exists", str_replace(path.to.s.obj, ".rds", ".addedInfo.rds")))
  }
}

