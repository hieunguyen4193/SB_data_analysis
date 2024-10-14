gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

samplesheet <- read.csv("/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/SampleSheet_for_DGE_and_CellChat.raw.csv")
all.s.obj <- unique(samplesheet$path)
for (file in all.s.obj){
  if (file.exists(str_replace(file, ".rds", ".addedConditions.rds")) == FALSE){
    s.obj <- readRDS(file)
    meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>% rowwise() %>%
      mutate(condition = paste(str_split(name, "_")[[1]][1:3], collapse = "_")) %>%
      mutate(sampletype = paste(str_split(name, "_")[[1]][2:3], collapse = "_")) %>%
      column_to_rownames("barcode")
    meta.data <- meta.data[row.names(s.obj@meta.data), ]
    s.obj <- AddMetaData(object = s.obj, col.name = "condition", metadata = meta.data$condition)
    s.obj <- AddMetaData(object = s.obj, col.name = "sampletype", metadata = meta.data$sampletype)
    all.conditions <- unique(s.obj$condition)
    print(all.conditions)
    saveRDS(s.obj, str_replace(file, ".rds", ".addedConditions.rds"))
    print(sprintf("finished saving %s", file))  
  } else {
    print(sprintf("%s exists", str_replace(file, ".rds", ".addedConditions.rds")))
  }
}
