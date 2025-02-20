gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/release"
scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering_SeuratV5.R"))



outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20250102"
path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/release"

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20250102"

all.PROJECTS <- c("SBharadwaj_20240318_Sample_1_4_7_8",
                  "SBharadwaj_20240318_Sample_1_4_7_8_2_5",
                  "SBharadwaj_20240318_Sample_3_6")

sub.cluster.idx <- "v0.1"

output.number <- "09_output"
# output.number <- "07_output"

if (output.number == "07_output"){
  source(file.path(path.to.project.src, "sub_clustering_indices.continued_from_06.R"))
  integration.cases <- c(FALSE, TRUE)
} else if (output.number == "09_output"){
  source(file.path(path.to.project.src, "sub_clustering_indices.continued_from_06.noReIntegration.R"))
  integration.cases <- c(FALSE)
}

for (PROJECT in all.PROJECTS){
  for (cont.sub.cluster.idx in names(sub_clusters[[PROJECT]][[sub.cluster.idx]])){
    for (integration in integration.cases){
      input.params <- list(
        PROJECT = PROJECT,
        sub.cluster.idx = sub.cluster.idx, 
        cont.sub.cluster.idx = cont.sub.cluster.idx,
        integration = integration,
        outdir = outdir,
        output.number = output.number
      )
      path.to.html.outputs <- file.path(outdir, PROJECT, "html_output", output.number, sub.cluster.idx, cont.sub.cluster.idx)
      dir.create(path.to.html.outputs, showWarnings = FALSE, recursive = TRUE)
      
      path.to.Rmd.file <- file.path(path.to.main.src, "07_continued_sub_clusters_from_06.Rmd")
      if (integration == FALSE){
        html_name <- basename(path.to.Rmd.file) %>% str_replace(".Rmd", sprintf(".%s.noIntegration.html", PROJECT))    
      } else {
        html_name <- basename(path.to.Rmd.file) %>% str_replace(".Rmd", sprintf(".%s.html", PROJECT))
      }
      
      html_name <- str_replace(html_name, ".html", sprintf(".%s.%s.html", sub.cluster.idx, cont.sub.cluster.idx))
      
      if (file.exists(file.path(path.to.html.outputs, html_name)) == FALSE){
        for (n in names(input.params)){
          print(sprintf("%s: %s", n, input.params[[n]]))
        }
        rmarkdown::render(input = path.to.Rmd.file,
                          output_file = html_name,
                          output_dir = path.to.html.outputs,
                          params = input.params)   
      } else {
        print(sprintf("File %s exists", file.path(path.to.html.outputs, html_name)))
      }
    }  
  }
}
