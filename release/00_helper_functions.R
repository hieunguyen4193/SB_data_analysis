`%ni%` = Negate(`%in%`)

#####----------------------------------------------------------------------#####
# function to create an interactive data table
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}


run_qc <- function(s.obj, object.name = "newQC"){
  all.QC <- list()
  
  # Number of cells obtained per experiment/sample
  all.QC$cell.counts.plot <- ggplot(s.obj@meta.data, 
                                    aes(x=name , fill=name)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    ggtitle("Number of cells in each dataset")
  
  # distribution of number of UMI in each sample
  all.QC$nCountRNA.distribution <- ggplot(s.obj@meta.data,
                                          aes(color=name, x=nCount_RNA, fill = name)) + 
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    ylab("Cell density") +
    geom_vline(xintercept = 500, color = "red") +
    ggtitle("Distribution of read depths in each sample")
  
  # distribution of number of features (genes detected) in each sample
  all.QC$nFeature_RNA.distribution <- ggplot(s.obj@meta.data,
                                             aes(color=name, x=nFeature_RNA, fill = name, y = ..scaled..)) + 
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    ylab("Cell density") +
    geom_vline(xintercept = 1000, color = "red") +
    xlim(1000, 10000) +
    ggtitle("Distribution of number of detected genes in each sample")
  
  
  # scatter plot showing the relation between cell read-depth and number of genes detected.
  ## with Mitochondrial percentage
  all.QC$nCount.vs.nFeature.MT <- ggplot(s.obj@meta.data, 
                                         aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm, formula = y ~ x) + # apply a linear regression to show the relation, if existed.
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    facet_wrap(~name) +
    ggtitle("Scatter plot: nCount_RNA vs. nFeature_RNA, cmap % Mitochondrial genes")
  
  ## with Ribosome percentage
  all.QC$nCount.vs.nFeature.Ribo <- ggplot(s.obj@meta.data, 
                                           aes(x=nCount_RNA, y=nFeature_RNA, color=percent.ribo)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm, formula = y ~ x) + # apply a linear regression to show the relation, if existed.
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    facet_wrap(~name) +
    ggtitle("Scatter plot: nCount_RNA vs. nFeature_RNA, cmap % Ribosome genes")
  
  # Complexity: 
  
  ## We can see the samples where we sequenced each cell less have a higher overall complexity, 
  # that is because we have not started saturating the sequencing for any given gene for these samples. 
  # Outlier cells in these samples might be cells that have a less complex RNA species than other cells. 
  # Sometimes we can detect contamination with low complexity cell types like red blood cells via this metric. 
  # Generally, we expect the novelty score to be above 0.80. 
  
  ## More expanations: they are looking for cells that have a low number of genes with a high number of UMI counts. 
  # This likely means that you only captured transcripts from a low number of genes, and simply sequenced transcripts 
  # from those lower number of genes over and over again. This could be because of the cell type 
  # (such as a red blood cell having little to no RNA as they mentioned), or some other strange artifact.
  
  # Compute the complexity and add it to the s.obj@meta.data
  s.obj@meta.data <- s.obj@meta.data %>% 
    mutate(log10GenesPerUMI = log10(nFeature_RNA) / log10(nCount_RNA))
  
  all.QC$complexity <- ggplot(s.obj@meta.data,
                              aes(x=log10GenesPerUMI, color = name, fill=name)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    geom_vline(xintercept = 0.8) +
    ggtitle("Complexity: Log10(nCount_RNA) / log10(nFeature_RNA)")
  
  # add new slot for all.QC into the existed SEURAT OBJECT. 
  s.obj@misc[[object.name]] <- all.QC
}

##### PRE-DEFINED CLUSTERS TO FOCUS ON IN CELL CHAT, DEPENDING ON THE EXPRESSION OF 
##### THE GENE FOXP3, CIITA, KLRG3, GATA3

predefined.clusters <- list(
  SBharadwaj_20240318_Sample_7_8 = list(
    Ciita = c(3, 5, 13, 14),
    Foxp3 = c(7, 11),
    Klrg1 = c(4, 7, 11),
    Gata3 = c(4, 7, 11)
  ),
  SBharadwaj_20240318_Sample_1_4 = list(
    Ciita = c(1, 3, 5, 12, 15),
    Foxp3 = c(7, 11),
    Klrg1 = c(4, 7, 11),
    Gata3 = c(2, 4, 7, 11)
  ),
  SBharadwaj_20240318_Sample_1_4_7_8 = list(
    Ciita = c(2, 3, 5, 13, 14, 16),
    Foxp3 = c(6, 12),
    Klrg1 = c(4, 6, 12),
    Gata3 = c(4, 6, 12)
  ),
  SBharadwaj_20240318_Sample_1_7 = list(
    Ciita = c(1, 2, 5, 13),
    Foxp3 = c(6, 11),
    Klrg1 = c(4, 6, 11),
    Gata3 = c(3, 4, 6, 11)
  ),
  SBharadwaj_20240318_Sample_2_3_5_6 = list(
    Ciita = c(0, 5, 6, 7, 9, 16, 17, 19),
    Foxp3 = c(),
    Klrg1 = c(8),
    Gata3 = c(12)
  ),
  SBharadwaj_20240318_Sample_2_5 = list(
    Ciita = c(0, 1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13),
    Foxp3 = c(),
    Klrg1 = c(),
    Gata3 = c(6)
  ),
  SBharadwaj_20240318_Sample_3_6 = list(
    Ciita = c(5, 23),
    Foxp3 = c(),
    Klrg1 = c(8),
    Gata3 = c()
  ),
  SBharadwaj_20240318_Sample_4_8 = list(
    Ciita = c(2, 3, 5, 13, 14, 15, 16),
    Foxp3 = c(7, 10),
    Klrg1 = c(4, 7, 10),
    Gata3 = c(1, 4, 7, 8, 10)
  )
)

excluded_cellchat_runs <- c("SBharadwaj_20240318_Sample_3_6",
                            "SBharadwaj_20240318_Sample_2_5",
                            "SBharadwaj_20240318_Sample_2_3_5_6"
                            )
#####----------------------------------------------------------------------#####
##### function to subset a cellchat object
#####----------------------------------------------------------------------#####
subset_cellchat <- function(object, selected.clusters){
  net <- object@net
  for (net.j in names(net)) {
    values <- net[[net.j]]
    if (net.j %in% c("prob","pval")) {
      values.new <- values[selected.clusters, selected.clusters, , drop = FALSE]
      net[[net.j]] <- values.new
    }
    if (net.j %in% c("count","sum","weight")) {
      values.new <- values[selected.clusters, selected.clusters, drop = FALSE]
      net[[net.j]] <- values.new
    }
  }
  net.subset <- net  
  return(net.subset)
}


