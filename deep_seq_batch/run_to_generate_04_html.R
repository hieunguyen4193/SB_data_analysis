##### clean up #####
gc()
rm(list = ls())

library(argparse)
my_random_seed <- 42
set.seed(my_random_seed)

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)
library(dplyr)
outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

samplesheet <- read.csv(file.path(path.to.main.src, "SampleSheet_for_DGE_and_CellChat.csv"))
samplesheet <- samplesheet %>% rowwise() %>%
  mutate(full.dataset.name = sprintf("%s_%s", PROJECT, dataset_name))
samplesheet <- samplesheet[!duplicated(samplesheet$full.dataset.name), ]

# write.table(samplesheet, file.path(path.to.main.src, "SampleSheet_for_DGE_and_CellChat_monocle3.csv"), row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

if ("monocle3" %in% installed.packages() == FALSE){
  devtools::install_github('cole-trapnell-lab/monocle3', force = TRUE, upgrade = FALSE)
} else if (packageVersion("monocle3") != "1.3.7"){
  devtools::install_github('cole-trapnell-lab/monocle3', force = TRUE, upgrade = FALSE)
}

if ("tradeSeq" %in% installed.packages() == FALSE){
  install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/tradeSeq_1.20.0.tar.gz", type = "source", repos = NULL)  
}

if ("monocle" %in% installed.packages()){
  remove.packages("monocle")  
}
library(monocle3)
path.to.rmd <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/04_run_monocle3.Rmd"

parser <- ArgumentParser()

parser$add_argument("-i", "--full_name", action="store",
                    help="Full name of the input project/dataset name")

args <- parser$parse_args()

full.name <- args$full_name

# for (full.name in unique(samplesheet$full.dataset.name)){
  
print(sprintf("Working on dataset %s", full.name))
PROJECT <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$PROJECT

path.to.save.html <- file.path(outdir, PROJECT, "html_output", "04_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
dataset_name <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$dataset_name
path.to.s.obj <- subset(samplesheet, samplesheet$full.dataset.name == full.name)$path
path.to.04.output <- file.path(path.to.main.output, "04_output", PROJECT, dataset_name)
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)
excluded.clusters <- NULL
use_partition <- FALSE
save.html.name <- sprintf("04_monocle3_usePartition_%s_%s_%s.html", use_partition, PROJECT, dataset_name)
if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
  rmarkdown::render(
    input = path.to.rmd,
    params = list(
      path.to.s.obj = path.to.s.obj,
      path.to.save.output = path.to.04.output,
      use_partition = use_partition,
      excluded.clusters = excluded.clusters
    ),
    output_file = save.html.name,
    output_dir = path.to.save.html
  )
}

# }

