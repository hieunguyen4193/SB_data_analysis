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

pair.samples <- list(
  `SBharadwaj_20240318_Sample_1_4_7_8_2_5` = list(
    `1` = c("treated_gut_CD45", "ctrl_gut_CD45"),
    `2` = c("treated_gut_myeloid", "ctrl_gut_myeloid"),
    `global.condition` = c("treated", "ctrl")
  ),
  `SBharadwaj_20240318_Sample_1_4_7_8` = list(
    `1` = c("treated_gut_CD45", "ctrl_gut_CD45")
  ),
  `SBharadwaj_20240318_Sample_3_6` = list(
    `1` = c("treated_liver_myeloid", "ctrl_liver_myeloid")
  )
)

path.to.rmd <- file.path(path.to.main.src, "03_CellChat_diff_analysis_2_samples.Rmd")

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
  
  path.to.s.obj <- samplesheet[row_i, ]$path
  
  path.to.s.obj <- str_replace(path.to.s.obj, ".rds", ".addedInfo.rds")
  
  path.to.save.html <- file.path(outdir, PROJECT, "html_output", "03_output_CellChat_diff", dataset_name)
  dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)
  
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  path.to.03.output <- file.path(path.to.main.output, "03_output")
  
  filter10cells <- "Filter10"
  
  all.pairs <- pair.samples[[PROJECT]]
  for (j in names(all.pairs)){
    if (j != "global.condition"){
      condition.name <- "condition"
    } else {
      condition.name <- "global.condition"
    }
    sample1 <- all.pairs[[j]][[1]]
    sample2 <- all.pairs[[j]][[2]]
    
    for (cluster.name in to.run.clusters){
      path.to.save.output <- file.path(path.to.03.output, dataset_name, cluster.name, condition.name, sprintf("%s_vs_%s", sample1, sample2))
      dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
      
      print(sprintf("Working on DGE analysis, PROJECT %s, dataset %s", PROJECT, dataset_name))
      print(sprintf("cluster name: %s, condition.name: %s", cluster.name, condition.name))
      print(sprintf("Perform CellChat diff. analysis %s vs %s", sample1, sample2))
      
      path.to.cellchat1 <- file.path(path.to.03.output, dataset_name, cluster.name, condition.name, sample1, 
                                     sprintf("CellChat_object.%s.%s.rds", sample1, filter10cells))
      path.to.cellchat2 <- file.path(path.to.03.output, dataset_name, cluster.name, condition.name, sample2, 
                                     sprintf("CellChat_object.%s.%s.rds", sample2, filter10cells))
      
      print(sprintf("Reading cell chat object for %s from %s", sample1, path.to.cellchat1))
      print(sprintf("Reading cell chat object for %s from %s", sample2, path.to.cellchat2))
      
      print(sprintf("Reading seurat object from %s", path.to.s.obj))
      print(sprintf("Saving output to %s", path.to.save.output))
      
      rmarkdown::render(input = path.to.rmd,
                        params = list(
                          path.to.s.obj = path.to.s.obj,
                          sample1 = sample1,
                          path.to.cellchat1 = path.to.cellchat1,
                          sample2 = sample2,
                          path.to.cellchat2 = path.to.cellchat2,
                          path.to.save.output = path.to.save.output,
                          filter10cells = filter10cells,
                          condition.name = condition.name,
                          cluster.name = cluster.name,
                          reduction.name = reduction.name
                        ),
                        output_dir = path.to.save.html,
                        output_file = sprintf("03_CellChat_%s_%s_Sample_%s_vs_%s.html",
                                              condition.name, cluster.name, sample1, sample2))
    }
  }
}
