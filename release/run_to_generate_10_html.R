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

library(argparse)
parser <- ArgumentParser()

parser$add_argument("-i", "--row_i", action="store",
                    help="Full name of the input project/dataset name")

args <- parser$parse_args()

row_i <- args$row_i

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

path.to.rmd <- file.path(path.to.main.src, "10_run_monocle3_for_single_sample.Rmd")

# for (row_i in seq(1, nrow(samplesheet))){
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
  
  path.to.save.html <- file.path(outdir, PROJECT, "html_output", "10_output", dataset_name)
  dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)
  
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  path.to.10.output <- file.path(path.to.main.output, "10_output")
  
  all.cases <- sample.list[[PROJECT]]
  
  for (condition.name in names(all.cases)){
    for (cluster.name in to.run.clusters){
      for (sample.id in all.cases[[condition.name]]){
        path.to.save.output <- file.path(path.to.10.output, dataset_name, cluster.name, condition.name, sample.id)
        dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
        
        print(sprintf("Working on monocle3 for single sample analysis, PROJECT %s, dataset %s", PROJECT, dataset_name))
        print(sprintf("cluster name: %s, condition.name: %s", cluster.name, condition.name))
        print(sprintf("Sample: %s", sample.id))
        
        if (file.exists(file.path(path.to.save.html, sprintf("10_monocle3_single_sample_%s_%s_%s.html", 
                                                             condition.name, 
                                                             cluster.name, 
                                                             sample.id))) == FALSE){
          excluded.clusters <- NULL
          use_partition <- FALSE
          
          input.params <- list(
            path.to.s.obj = path.to.s.obj,
            path.to.save.output = path.to.save.output,
            use_partition = use_partition,
            excluded.clusters = excluded.clusters,
            reduction.name = reduction.name, 
            cluster.name = cluster.name,
            sample.id = sample.id,
            condition.name = condition.name
          )
          
          print("----------------------------------")
          print("List of input params:")
          for (p in names(input.params)){
            print(sprintf("%s = %s", p, input.params[[p]]))
          }
          print("----------------------------------")
          
          rmarkdown::render(input = path.to.rmd,
                            params = input.params,
                            output_dir = path.to.save.html,
                            output_file = sprintf("10_monocle3_single_sample_%s_%s_%s.html", 
                                                  condition.name, cluster.name, sample.id))
        } else {
          print(sprintf("File %s exists", file.path(path.to.save.html, sprintf("10_monocle3_single_sample_%s_%s_%s.html", 
                                                                               condition.name, 
                                                                               cluster.name, 
                                                                               sample.id))))
        }
      }
    }
  }
# }
