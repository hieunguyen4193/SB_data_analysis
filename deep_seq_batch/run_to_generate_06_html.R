gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"
path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

# PROJECT <- "SBharadwaj_20240318_Sample_1_4_7_8_2_5"
# PROJECT <- "SBharadwaj_20240318_Sample_3_6"
# PROJECT <- "SBharadwaj_20240318_Sample_1_4_7_8"

for (PROJECT in c("SBharadwaj_20240318_Sample_3_6",
                  "SBharadwaj_20240318_Sample_1_4_7_8_2_5",
                  "SBharadwaj_20240318_Sample_1_4_7_8")){
  sub.cluster.idx <- "v0.1"
  
  path.to.html.outputs <- file.path(outdir, PROJECT, "html_output")
  dir.create(path.to.html.outputs, showWarnings = FALSE, recursive = TRUE)
  
  path.to.Rmd.file <- file.path(path.to.main.src, "06_sub_clustering_analysis.Rmd")
  html_name <- basename(path.to.Rmd.file) %>% str_replace(".Rmd", sprintf(".%s.%s.html", PROJECT, sub.cluster.idx))
  
  for (integration in c(TRUE, FALSE)){
    if (integration == FALSE){
      html_name <- str_replace(html_name, ".html", ".noIntegration.html")
    }
    if (file.exists(file.path(path.to.html.outputs, html_name)) == FALSE){
      rmarkdown::render(input = path.to.Rmd.file,
                        output_file = html_name,
                        output_dir = path.to.html.outputs,
                        params = list(PROJECT = PROJECT, 
                                      sub.cluster.idx = sub.cluster.idx,
                                      integration = integration))   
    } else {
      print(sprintf("file %s existed", PROJECT))
    }
  }
}

