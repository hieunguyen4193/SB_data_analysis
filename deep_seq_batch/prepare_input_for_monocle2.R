##### clean up #####
gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)
if ("monocle" %in% installed.packages() == FALSE){
  BiocManager::install("monocle", update = FALSE)
}
library(monocle)

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"
samplesheet <- read.csv(file.path(path.to.main.src, "SampleSheet_for_DGE_and_CellChat.csv"))
samplesheet <- samplesheet %>% rowwise() %>%
  mutate(full.dataset.name = sprintf("%s_%s", PROJECT, dataset_name))
samplesheet <- samplesheet[!duplicated(samplesheet$full.dataset.name), ]

for (full.name in unique(samplesheet$full.dataset.name)){
  print(sprintf("Working on dataset %s", full.name))
  PROJECT <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$PROJECT
  dataset_name <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$dataset_name
  path.to.input.s.obj <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$path
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  path.to.monocle2.input <- file.path(path.to.main.output, "08_output", "monocle2_input")
  dir.create(path.to.monocle2.input, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(path.to.monocle2.input, 
                            sprintf("%s.%s.monocle2.rds", PROJECT, dataset_name))) == FALSE){
    print(sprintf("File %s does not exists, generating ...",
                  file.path(path.to.monocle2.input, 
                            sprintf("%s.%s.monocle2.rds", PROJECT, dataset_name))))
    s.obj <- readRDS(path.to.input.s.obj)
    
    data <- GetAssayData(s.obj, slot = "data", assay = "SCT")
    
    pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)
    
    fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fd)
    
    library(monocle)
    monocle.obj <- newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  lowerDetectionLimit = 0.5,
                                  expressionFamily = negbinomial.size())
    saveRDS(monocle.obj, file.path(path.to.monocle2.input, 
                                   sprintf("%s.%s.monocle2.rds", PROJECT, dataset_name)))
  } else {
    print(sprintf("File %s exists", file.path(path.to.monocle2.input, 
                                              sprintf("%s.%s.monocle2.rds", PROJECT, dataset_name))))
  }

}

