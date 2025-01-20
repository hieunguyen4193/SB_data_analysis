#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

#####----------------------------------------------------------------------#####
##### install packages for monocle3
#####----------------------------------------------------------------------#####
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
if ("monocle3" %in% installed.packages() == FALSE){
  devtools::install_github('cole-trapnell-lab/monocle3', force = TRUE, upgrade = FALSE)
} else if (packageVersion("monocle3") != "1.3.7"){
  devtools::install_github('cole-trapnell-lab/monocle3', force = TRUE, upgrade = FALSE)
}
if ("tradeSeq" %in% installed.packages() == FALSE){
  install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/tradeSeq_1.20.0.tar.gz", type = "source", repos = NULL)  
}
if ("monocle" %in% installed.packages()){
  remove.packages("monocle")  
}
library(monocle3)

#####----------------------------------------------------------------------#####
##### MAIN RUN
#####----------------------------------------------------------------------#####

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

path.to.rmd <- file.path(path.to.main.src, "04_run_monocle3.Rmd")

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
  
  if (PROJECT == "SBharadwaj_20240318_Sample_3_6"){
    to.run.clusters <- setdiff(to.run.clusters, "cell.annotation")
  }
  
  all.cases <- sample.list[[PROJECT]]
  for (cluster.name in to.run.clusters){
    path.to.s.obj <- samplesheet[row_i, ]$path
    path.to.s.obj <- str_replace(path.to.s.obj, ".rds", ".addedInfo.rds")
    
    path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
    path.to.05.output <- file.path(path.to.main.output, "05_output")
    
    path.to.save.output <- file.path(path.to.05.output, dataset_name, cluster.name, "monocle2_input")
    dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
    
    if (file.exists(file.path(path.to.save.output, "monocle2_obj.rds")) == FALSE){
      data <- GetAssayData(s.obj, slot = "data", assay = "SCT")
      
      pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)
      
      fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
      fd <- new('AnnotatedDataFrame', data = fd)
      
      library(monocle)
      monocle.obj <- newCellDataSet(data,
                                    phenoData = pd,
                                    featureData = fd,
                                    lowerDetectionLimit = 0.5,
                                    expressionFamily = negbinomial.size())
      saveRDS(object = monocle.obj, file.path(path.to.save.output, "monocle2_obj.rds"))
    } else {
      print(sprintf("file %s exists", file.path(path.to.save.output, "monocle2_obj.rds")))
    }
  }
}