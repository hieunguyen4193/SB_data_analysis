gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

reticulate::install_python(version = '3.8')
reticulate::py_install(packages = 'umap-learn')

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

samplesheet <- read.csv("/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/SampleSheet_for_DGE_and_CellChat.csv")

all.s.objs <- unique(samplesheet$path)
sample.list <- list(
  SBharadwaj_20240318_Sample_1_4_7_8_2_5 = c( "ctrl_gut_CD45_1",
                                              "treated_gut_CD45_1",
                                              "treated_gut_CD45_2",
                                              "ctrl_gut_CD45_2",
                                              "ctrl_gut_myeloid",
                                              "treated_gut_myeloid"),
  SBharadwaj_20240318_Sample_3_6 <- c("ctrl_liver_myeloid", 
                                      "treated_liver_myeloid")
)

# for (p in unique(samplesheet$PROJECT)){
#     sample.list[[p]] <- c(subset(samplesheet, samplesheet$PROJECT == p)$sample1, 
#                           subset(samplesheet, samplesheet$PROJECT == p)$sample2) %>% unique()
# }

filter10cells <- "Filter10"

path.to.rmd.file <- file.path(path.to.main.src, "03_CellChat_diff_analysis_2_samples.Rmd")
for (i in seq(1, nrow(samplesheet))){
    tmp.samplesheet <- samplesheet[i, ]
    path.to.s.obj <- tmp.samplesheet$path
    dataset_name <- tmp.samplesheet$dataset_name
    PROJECT <- tmp.samplesheet$PROJECT
    sample1 <- tmp.samplesheet$sample1
    sample2 <- tmp.samplesheet$sample2

    path.to.html.outputs <- file.path(outdir,
                                      "SeuratV5",
                                      PROJECT,
                                      "html_output",
                                      dataset_name,
                                      "03_output")
    html.filename <- sprintf("%s_vs_%s.CellChat.html", sample1, sample2)
    path.to.main.output <- file.path(outdir, "SeuratV5", PROJECT, "data_analysis")
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
          filter10cells = filter10cells
        ),
        output_file = html.filename,
        output_dir = path.to.html.outputs
      )
    }
  }
