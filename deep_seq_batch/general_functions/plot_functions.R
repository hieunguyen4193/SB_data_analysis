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

# path to the RDS file
path.to.s.obj <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318/SBharadwaj_20240318_Sample_1_4_7_8_2_5/data_analysis/01_output/Project_SBharadwaj_20240318_Sample_1_4_7_8_2_5.addedCellAnnotation.rds"

# path to save output svg file
path.to.save.output <- "/home/hieunguyen/CRC1382/outdir/tmp"
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(path.to.s.obj)

# Idents(s.obj) <- "cca.cluster.0.5" # <<<< CHANGE HERE to use cluster number.
Idents(s.obj) <- "cell.annotation" # <<<<< CHANGE HERE to switch to cell annotation
# update 12.11.2024: cell annotation available only in samples 1 4 7 8 2 5 and 1 4 7 8. 

# set default assay =  SCT
DefaultAssay(s.obj) < "SCT"

# genes to show on the UMAP, violin plot
plot.genes <- c("Cox20",
                "Hnrnpu",
                "Efcab2",
                "Smyd3" )
  
##### UMAP plot, all samples
umap.plot <- DimPlot(object = s.obj, reduction = "cca_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, pt.size = 1)

##### UMAP plot, splitted by conditions
umap.plot.conditions <- DimPlot(object = s.obj, reduction = "cca_UMAP", label = TRUE, label.box = FALSE, repel = TRUE, split.by = "condition", ncol = 2, 
                                pt.size = 1, 
                                label.size = 12)

##### UMAP plot, splitted by samples
umap.plot.sample <- DimPlot(object = s.obj, reduction = "cca_UMAP", label = TRUE, label.box = FALSE, repel = TRUE, split.by = "name", ncol = 2, pt.size = 1)

##### feature plot
feature.plot.sample <- FeaturePlot(object = s.obj, reduction = "cca_UMAP", label = TRUE, repel = TRUE, ncol = 2, pt.size = 1, features = plot.genes)

##### violin plot
violin.plot <- VlnPlot(object = s.obj, features = plot.genes, ncol = 2, pt.size = 1)  # set pt.size = 0 to remove points on violin plots

##### violin plot, multiple samples
violin.plot.condition1 <- VlnPlot(object = s.obj, features = plot.genes, ncol = 2, pt.size = 1, split.by = "condition")  # set pt.size = 0 to remove points on violin plots
violin.plot.condition2 <- VlnPlot(object = s.obj, features = plot.genes, ncol = 2, pt.size = 1, group.by = "condition")  # set pt.size = 0 to remove points on violin plots
violin.plot.condition3 <- VlnPlot(object = s.obj, features = plot.genes, ncol = 2, pt.size = 1, group.by = "name")

##### save function
choose_your_filename <- "test.svg"
ggsave(plot = umap.plot.sample, filename = choose_your_filename, path = path.to.save.output, device = "svg", width = 14, height = 10, dpi = 300)


unique(s.obj$name)
unique(s.obj$condition)
