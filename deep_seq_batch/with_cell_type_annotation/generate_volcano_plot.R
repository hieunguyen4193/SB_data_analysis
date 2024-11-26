#####----------------------------------------------------------------------#####
##### generate volcano plot
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())
library(ggpubr)
path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path.to.save.output <- file.path(outdir, "volcano_plots")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"
##### read the sample sheet containing paths to seurat object with cell annotations.
samplesheet <- read.csv(file.path(path.to.project.src, 
                        "/with_cell_type_annotation/SampleSheet_for_DGE_CellChat.CellAnnotated.csv"))

PROJECT <- "SBharadwaj_20240318_Sample_1_4_7_8"

sample1 <- "treated_gut_CD45"
sample2 <- "ctrl_gut_CD45"
dataset.name <- "full"

path.to.save.output <- file.path(path.to.save.output, 
                                 PROJECT,
                                 dataset.name, 
                                 sprintf("%s_vs_%s", sample1, sample2))
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
# get cell annotations
source(file.path(path.to.project.src, "cell_type_annotation.R"))
cell.annotations <- all.annotations[[PROJECT]]
path.to.dge.obj <- file.path(outdir, 
                             PROJECT, 
                             "data_analysis",
                             "02_output", 
                             dataset.name,
                             sprintf("%s_vs_%s", sample1, sample2), 
                             sprintf("sample_%s_vs_%s.raw.single_cell.rds", 
                                     sample1, 
                                     sample2))
dge.obj <- readRDS(path.to.dge.obj)

new.dge.obj <- list()
for (i in names(dge.obj)){
  new.dge.obj[[cell.annotations[[str_replace(i, "cluster_", "")]]]] <- dge.obj[[i]]
}

cutoff.adjp <- 0.05

# i <- names(new.dge.obj)[[1]]
for (i in names(new.dge.obj)){
  input.df <- new.dge.obj[[i]]
  top10.up <- input.df %>% subset(p_val_adj <= cutoff.adjp) %>% arrange(desc(avg_log2FC)) %>% head(10) %>% pull(Gene)
  top10.down <- input.df %>% subset(p_val_adj <= cutoff.adjp) %>% arrange(avg_log2FC) %>% head(10) %>% pull(Gene)
  
  input.df <- input.df %>% mutate(abs.log2FoldChange = abs(avg_log2FC))
  input.df <- input.df %>% rowwise() %>%
    mutate(sig = ifelse(p_val_adj <= cutoff.adjp, "sig", "non.sig")) %>%
    mutate(show.gene.name = ifelse(Gene %in% c(top10.up, top10.down), Gene, NA))
  
  volcano.plot <- ggplot(data=input.df, 
                         aes(x=avg_log2FC, y=-log10(p_val_adj), col=sig, label=show.gene.name)) + 
    geom_point(size = 2) +
    theme_pubr() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    geom_vline(xintercept=0, linetype = "dotted") +
    # geom_vline(xintercept=c(-1, 1), col="#9a9fa6", linetype='dotted') +
    # geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6", linetype='dotted') +
    geom_text_repel(size = 8) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 25),
          axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25)) +
    xlab("logFC") + ylab("-log10(p-value)") + 
    xlim(c(-max(input.df$abs.log2FoldChange), max(input.df$abs.log2FoldChange))) 
    
  ggsave(plot = volcano.plot, filename = sprintf("%s_%s_vs_%s.volcano_plot.svg", i, sample1, sample2), 
         path = path.to.save.output, dpi = 300, width = 14, height = 10)
}

