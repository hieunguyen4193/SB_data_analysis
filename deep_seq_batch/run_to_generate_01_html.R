#####----------------------------------------------------------------------#####
##### 01
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/home/hieunguyen/CRC1382/outdir"

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

all.PROJECTS <- c("SBharadwaj_20240318_Sample_3_6")

for (PROJECT in all.PROJECTS){
  path.to.html.outputs <- file.path(outdir, "SeuratV5", PROJECT, "html_output")
  dir.create(path.to.html.outputs, showWarnings = FALSE, recursive = TRUE)
  
  path.to.Rmd.file <- file.path(path.to.main.src, "01_preliminary_analysis_all_samples_integrated.Rmd")
  html_name <- basename(path.to.Rmd.file) %>% str_replace(".Rmd", sprintf(".%s.html", PROJECT))
  if (file.exists(file.path(path.to.html.outputs, html_name)) == FALSE){
    rmarkdown::render(input = path.to.Rmd.file,
                      output_file = html_name,
                      output_dir = path.to.html.outputs,
                      params = list(PROJECT = PROJECT))   
  }  
}

