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

path.to.rmd <- file.path(path.to.main.src, "03_CellChat_general_analysis.Rmd")

row_i <- 1

PROJECT <- samplesheet[row_i, ]$PROJECT
dataset_name <- samplesheet[row_i, ]$dataset_name
path.to.s.obj <- samplesheet[row_i, ]$path

path.to.s.obj <- str_replace(path.to.s.obj, ".rds", ".addedInfo.rds")

path.to.save.html <- file.path(outdir, PROJECT, "html_output", "03_output", dataset_name)
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.03.output <- file.path(path.to.main.output, "03_output")

all.cases <- sample.list[[PROJECT]]

for (condition.name in names(all.cases)){
  for (cluster.name in c("cca.cluster.0.5", "cell.annotation")){
    for (sample.id in all.cases[[condition.name]])
      path.to.save.output <- file.path(path.to.03.output, dataset_name, cluster.name, condition.name, sample.id)
      dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
      
      print(sprintf("Working on CellChat analysis, PROJECT %s, dataset %s", PROJECT, dataset_name))
      print(sprintf("cluster name: %s, condition.name: %s", cluster.name, condition.name))
      print(sprintf("Sample: %s", sample.id))
      
      rmarkdown::render(input = path.to.rmd,
                        params = list(
                          sample.id = sample.id,
                          path.to.s.obj = path.to.s.obj,
                          path.to.save.output = path.to.save.output,
                          filter10cells = "Filter10",
                          condition.name = condition.name,
                          cluster.name = cluster.name
                        ),
                        output_dir = path.to.save.html,
                        output_file = sprintf("03_CellChat_%s_%s_%s.html", 
                                              condition.name, cluster.name, sample.id))
  }
}
