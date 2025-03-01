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
# write.csv(samplesheet, file.path(path.to.main.src, "SampleSheet_all_seurat_objects.modified.csv"))

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

path.to.rmd <- file.path(path.to.main.src, "02_DGE_analysis.Rmd")

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
  
  path.to.save.html <- file.path(outdir, PROJECT, "html_output", "02_output", dataset_name)
  dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)
  
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  path.to.02.output <- file.path(path.to.main.output, "02_output")
  
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
      path.to.save.output <- file.path(path.to.02.output, dataset_name, cluster.name, condition.name, sprintf("%s_vs_%s", sample1, sample2))
      dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
      
      print(sprintf("Working on DGE analysis, PROJECT %s, dataset %s", PROJECT, dataset_name))
      print(sprintf("cluster name: %s, condition.name: %s", cluster.name, condition.name))
      
      if (file.exists(file.path(path.to.save.html, sprintf("02_DGE_analysis_%s_vs_%s_%s_%s.html",
                                                           sample1,
                                                           sample2, 
                                                           cluster.name,
                                                           condition.name))) == FALSE){
        rmarkdown::render(input = path.to.rmd, 
                          params = list(
                            sample1 = sample1,
                            sample2 = sample2,
                            path.to.s.obj = path.to.s.obj,
                            path.to.save.output = path.to.save.output,
                            cluster.name = cluster.name,
                            condition.name = condition.name,
                            reduction.name = reduction.name
                          ),
                          output_file = sprintf("02_DGE_analysis_%s_vs_%s_%s_%s.html",
                                                sample1,
                                                sample2, 
                                                cluster.name,
                                                condition.name),
                          output_dir = path.to.save.html)
      } else {
        print(sprintf("File %s exists",file.path(path.to.save.html, sprintf("02_DGE_analysis_%s_vs_%s_%s_%s.html",
                                                                            sample1,
                                                                            sample2, 
                                                                            cluster.name,
                                                                            condition.name))))
      }
    }
  }
}


