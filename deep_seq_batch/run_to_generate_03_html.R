gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

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
                                              "treated_gut_myeloid",
                                              "ctrl_gut_CD45",
                                              "treated_gut_CD45"),
  SBharadwaj_20240318_Sample_3_6 = c("ctrl_liver_myeloid", 
                                      "treated_liver_myeloid")
)

# for (p in unique(samplesheet$PROJECT)){
#     sample.list[[p]] <- c(subset(samplesheet, samplesheet$PROJECT == p)$sample1, 
#                           subset(samplesheet, samplesheet$PROJECT == p)$sample2) %>% unique()
# }

filter10cells <- "Filter10"

path.to.rmd.file <- file.path(path.to.main.src, "03_CellChat_general_analysis.Rmd")
for (path.to.s.obj in all.s.objs){
  PROJECT <- str_split(path.to.s.obj, "/data_analysis")[[1]][[1]] %>% basename()
  dataset_name <- subset(samplesheet, samplesheet$path == path.to.s.obj)$dataset_name %>% unique()
  if (length(dataset_name) == 1){
    dataset_name <- dataset_name[[1]]
  } else {
    print("error, more than 1 dataset_name")
  }
  dir.create(file.path(outdir, "SeuratV5", PROJECT, "html_output", dataset_name, "03_output"), showWarnings = FALSE, recursive = TRUE)  
  sample.in.data <- sample.list[[PROJECT]]
  for (sample.id in sample.in.data){
    if (sample.id %in% c("ctrl_gut_CD45", "treated_gut_CD45")){
      chosen.group <- "condition"
    } else {
      chosen.group <- "name"
    }
    path.to.html.outputs <- file.path(outdir,
                                      "SeuratV5",
                                      PROJECT,
                                      "html_output",
                                      "03_output",
                                      dataset_name
                                      )
    html.filename <- sprintf("%s.CellChat.html", sample.id)
    path.to.main.output <- file.path(outdir, "SeuratV5", PROJECT, "data_analysis" )
    path.to.03.output <- file.path(path.to.main.output, "03_output")
    dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
    path.to.save.output <- file.path(path.to.03.output, dataset_name, sample.id)
    
    if (file.exists(file.path(path.to.html.outputs, html.filename)) == FALSE){
      input.params <- list(
        sample.id = sample.id,
        path.to.s.obj = path.to.s.obj,
        path.to.save.output = path.to.save.output,
        filter10cells = filter10cells,
        chosen.group = chosen.group
      )
      for (i in names(input.params)){
        print(sprintf("%s: %s", i, input.params[[i]]))
      }
      rmarkdown::render(
        input = path.to.rmd.file,
        params = input.params,
        output_file = html.filename,
        output_dir = path.to.html.outputs
      )
    }
  }
}
