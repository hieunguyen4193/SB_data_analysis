#####----------------------------------------------------------------------#####
##### 01
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/release"

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20250102"
PROJECT <- "SBharadwaj_20240318_Sample_1_4_7_8_2_5"

path.to.html.outputs <- file.path(outdir, PROJECT, "html_output")
dir.create(path.to.html.outputs, showWarnings = FALSE, recursive = TRUE)

path.to.Rmd.file <- file.path(path.to.main.src, "CUSTOM_reclustering_DCs_Cd3g_genes.Rmd")

for (integration in c(TRUE, FALSE)){
  if (integration == FALSE){
    html_name <- str_replace(path.to.Rmd.file, ".Rmd", sprintf(".noReIntegration.html"))  
  } else {
    html_name <- str_replace(path.to.Rmd.file, ".Rmd", sprintf(".html"))
  }
  
  if (file.exists(file.path(path.to.html.outputs, html_name)) == FALSE){
    rmarkdown::render(input = path.to.Rmd.file,
                      output_file = html_name,
                      output_dir = path.to.html.outputs,
                      params = list(integration = integration))   
  }  
}

