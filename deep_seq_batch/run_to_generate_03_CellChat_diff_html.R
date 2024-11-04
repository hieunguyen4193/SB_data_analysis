gc()
rm(list = ls())

##### rerun Cellchat after tmp fixing bugs
# devtools::install_github("https://github.com/hieunguyen4193/CellChat", upgrade = FALSE)
# reticulate::install_python(version = '3.8')
# reticulate::py_install(packages = 'umap-learn')
##### BUGS in cellchat: https://github.com/jinworks/CellChat/issues/159
##### temporary fix --> https://github.com/jinworks/CellChat/issues/202
# cellchat.sample1@data.smooth <- cellchat.sample1@data.project
# cellchat.sample2@data.smooth <- cellchat.sample2@data.project

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

samplesheet <- read.csv("/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/SampleSheet_for_DGE_and_CellChat.csv")

all.s.objs <- unique(samplesheet$path)

filter10cells <- "Filter10"

path.to.rmd.file <- file.path(path.to.main.src, "03_CellChat_diff_analysis_2_samples.Rmd")
for (i in seq(1, nrow(samplesheet))){
    tmp.samplesheet <- samplesheet[i, ]
    path.to.s.obj <- tmp.samplesheet$path
    dataset_name <- tmp.samplesheet$dataset_name
    PROJECT <- tmp.samplesheet$PROJECT
    sample1 <- tmp.samplesheet$sample1
    sample2 <- tmp.samplesheet$sample2
    if (grepl("_gut_CD45", sample1) == TRUE & grepl("_gut_CD45", sample2) == TRUE){
      chosen.group <- "condition"
    } else {
      chosen.group <- "name"
    }
    path.to.html.outputs <- file.path(outdir,
                                      PROJECT,
                                      "html_output",
                                      "03_output",
                                      dataset_name
                                      )
    html.filename <- sprintf("%s_vs_%s.CellChat.html", sample1, sample2)
    path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
    path.to.03.output <- file.path(path.to.main.output, "03_output")
    dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
    path.to.save.output <- file.path(path.to.03.output, dataset_name, sprintf("%s_vs_%s", sample1, sample2))
    path.to.cellchat1 <- file.path(path.to.03.output, dataset_name, sample1, sprintf("CellChat_object.%s.Filter10.rds", sample1))
    path.to.cellchat2 <- file.path(path.to.03.output, dataset_name, sample2, sprintf("CellChat_object.%s.Filter10.rds", sample2))
    dir.create(path.to.html.outputs, showWarnings = FALSE, recursive = TRUE)
    if (file.exists(file.path(path.to.html.outputs, html.filename)) == FALSE){
      rmarkdown::render(
        input = path.to.rmd.file,
        params = list(
          path.to.s.obj = path.to.s.obj,
          sample1 = sample1,
          path.to.cellchat1 = path.to.cellchat1,
          sample2 = sample2,
          path.to.cellchat2 = path.to.cellchat2,
          path.to.save.output = path.to.save.output,
          filter10cells = filter10cells,
          chosen.group = chosen.group
        ),
        output_file = html.filename,
        output_dir = path.to.html.outputs
      )
    }
  }
