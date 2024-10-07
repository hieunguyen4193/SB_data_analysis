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

all.PROJECTS <- c(
  "SBharadwaj_20240318_Sample_2_5",
  "SBharadwaj_20240318_Sample_1_4", 
  "SBharadwaj_20240318_Sample_3_6",
  "SBharadwaj_20240318_Sample_1_4_7_8",
  "SBharadwaj_20240318_Sample_4_8",
  "SBharadwaj_20240318_Sample_1_7",
  "SBharadwaj_20240318_Sample_7_8",
  "SBharadwaj_20240318_Sample_2_3_5_6"
)

if ("CellChat" %in% installed.packages() == FALSE){
  devtools::install_github("jinworks/CellChat")
}
library(CellChat)

for (PROJECT in all.PROJECTS){
  print("#####-----------------------------------------------------------------#####")
  print(sprintf("working on %s", PROJECT))
  print("#####-----------------------------------------------------------------#####")
  
  path.to.main.input <- file.path(outdir, PROJECT)
  path.to.main.output <- file.path(path.to.main.input, "data_analysis")
  dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.01.output <- file.path(path.to.main.output, "01_output")
  path.to.02.output <- file.path(path.to.main.output, "02_output")
  path.to.03.output <- file.path(path.to.main.output, "03_output")
  dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
  
  s.obj <- readRDS(file.path(path.to.01.output,sprintf("Project_%s.rds", PROJECT)))
  
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>% rowwise() %>%
    mutate(condition = str_split(name, "_")[[1]][[1]]) %>%
    mutate(organ = str_split(name, "_")[[1]][[2]]) %>%
    mutate(celltype = str_split(name, "_")[[1]][[3]])  %>% 
    column_to_rownames("barcode")
  
  meta.data <- meta.data[row.names(s.obj@meta.data),]
  
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$condition, col.name = "condition")
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$organ, col.name = "organ")
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$celltype, col.name = "celltype")
  
  cells <- hash()
  for (input.gene in c("Ciita", "Foxp3", "Klrg1", "Gata3")){
    tmp <- GetAssayData(s.obj, assay = "RNA", slot = "data")[input.gene, ]  
    cells[[input.gene]] <- tmp[tmp > 0] %>% names()
  }
  
  group1.1 <- intersect(cells$Foxp3, cells$Klrg1)
  group1.2 <- intersect(cells$Foxp3, cells$Gata3)
  
  group1.cells <- unique(c(group1.1, group1.2))
  
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(group = ifelse(barcode %in% group1.cells, "Group1", ifelse(barcode %in% cells$Ciita, "Group2", "others"))) %>%
    column_to_rownames("barcode")
  
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$group, col.name = "group")
  
  if (file.exists(file.path(path.to.03.output, sprintf("%s.cellchat.rds", PROJECT))) == FALSE){
    data.input <- s.obj[["RNA"]]@data 
    Idents(s.obj) <- "group"
    labels <- Idents(s.obj)
    meta <- data.frame(labels = labels, row.names = names(labels)) 
    cellchat <- createCellChat(object = s.obj, group.by = "ident", assay = "RNA")
    
    
    ##### Select the CellChat database: human or mouse
    CellChatDB <- CellChatDB.mouse 
    showDatabaseCategory(CellChatDB)
    
    CellChatDB.use <- subsetDB(CellChatDB) # use all, except non protein signaling
    cellchat@DB <- CellChatDB.use
    
    ##### data in cellchat object pre-processing
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 4) # do parallel
    if ("presto" %in% installed.packages() == FALSE){
      devtools::install_github('immunogenomics/presto')
    }
    
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    ##### Compute the communication prob
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    ##### infer cell-cell communication at a signaling pathway
    cellchat <- computeCommunProbPathway(cellchat)
    
    ##### aggregate cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    
    saveRDS(object = cellchat, file.path(path.to.03.output, sprintf("%s.cellchat.rds", PROJECT)))  
  } 
}

