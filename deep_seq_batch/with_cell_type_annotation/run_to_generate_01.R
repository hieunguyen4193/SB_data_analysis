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

all.PROJECTS <- c("SBharadwaj_20240318_Sample_1_4_7_8", 
                  "SBharadwaj_20240318_Sample_1_4_7_8_2_5")

for (PROJECT in all.PROJECTS){
  path.to.html.outputs <- file.path(outdir, PROJECT, "html_output_with_cell_annotations")
  dir.create(path.to.html.outputs, showWarnings = FALSE, recursive = TRUE)
  
  path.to.Rmd.file <- file.path(path.to.main.src, 
                                "with_cell_type_annotation", 
                                "01_preliminiary_analysis.Rmd")
  html_name <- basename(path.to.Rmd.file) %>% str_replace(".Rmd", sprintf(".%s.html", PROJECT))
  if (file.exists(file.path(path.to.html.outputs, html_name)) == FALSE){
    rmarkdown::render(input = path.to.Rmd.file,
                      output_file = html_name,
                      output_dir = path.to.html.outputs,
                      params = list(
                        PROJECT = PROJECT
                        ))   
  }  
}

