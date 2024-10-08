#####----------------------------------------------------------------------#####
##### 01
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

samplesheet <- read.csv("/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/SampleSheet_for_DGE_and_CellChat.csv")

for (i in seq(1, nrow(samplesheet))){
  PROJECT <- samplesheet[i, "PROJECT"]
  dataset_name <- samplesheet[i, "dataset_name"]
  sample1 <- samplesheet[i, "sample1"]
  sample2 <- samplesheet[i, "sample2"]
  path.to.main.output <- file.path(outdir, "SeuratV5", PROJECT, "data_analysis" )
  path.to.02.output <- file.path(path.to.main.output, "02_output")
  dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)
  path.to.save.output <- file.path(path.to.02.output, dataset_name, sprintf("%s_vs_%s", sample1, sample2))
  dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
  path.to.s.obj <- samplesheet[i, "path"]
  path.to.html.outputs <- file.path(outdir, 
                                    "SeuratV5", 
                                    PROJECT, 
                                    "html_output", 
                                    "02_output",
                                    dataset_name, 
                                    sprintf("%s_vs_%s", sample1, sample2))
  dir.create(path.to.html.outputs, showWarnings = FALSE, recursive = TRUE)
  
  path.to.Rmd.file <- file.path(path.to.main.src, "02_DGE_analysis.Rmd")
  html_name <- basename(path.to.Rmd.file) %>% str_replace(".Rmd", sprintf("DGE_%s_vs_%s.html", sample1, sample2))
  if (file.exists(file.path(path.to.html.outputs, html_name)) == FALSE){
    rmarkdown::render(input = path.to.Rmd.file,
                      output_file = html_name,
                      output_dir = path.to.html.outputs,
                      params = list(sample1 = sample1,
                                    sample2 = sample2,
                                    path.to.s.obj = path.to.s.obj,
                                    path.to.save.output = path.to.save.output))   
  }  
}

