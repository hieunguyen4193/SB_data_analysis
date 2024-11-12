gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### Install new packages
#####----------------------------------------------------------------------#####
list.of.packages <- c("Seurat",
                      "SingleCellExperiment",
                      "optparse", 
                      "comprehenr", 
                      "tidyverse", 
                      "ggplot2", 
                      "SoupX",
                      "comprehenr",
                      "DoubletFinder",
                      "vroom",
                      "hash",
                      "DT",
                      "janitor",
                      "knitr",
                      "circlize",
                      "formattable",
                      "htmlwidgets",
                      "plotly",
                      "stringr"
)

bioc.packages <- c("celda", 
                   "BiocSingular", 
                   "PCAtools", 
                   "SingleCellExperiment",
                   "scRepertoire", 
                   "sctransform")

# Check if packages are installed 

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.bioc.packages <- bioc.packages[!(bioc.packages %in% installed.packages()[,"Package"])]

# Install new packages 
if(length(new.packages)) install.packages(new.packages)
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages)

# Import all packages 
package_loading_Status <- lapply(list.of.packages, 
                                 require, 
                                 character.only = TRUE)

package_loading_Status_bioc <- lapply(bioc.packages, 
                                      require, 
                                      character.only = TRUE)

#####----------------------------------------------------------------------#####
##### MAIN FUNCTIONS
#####----------------------------------------------------------------------#####

path.to.s.obj <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318/SBharadwaj_20240318_Sample_1_4_7_8_2_5/data_analysis/01_output/Project_SBharadwaj_20240318_Sample_1_4_7_8_2_5.addedCellAnnotation.rds"
path.to.save.output <- "/home/hieunguyen/CRC1382/outdir/tmp"
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(path.to.s.obj)