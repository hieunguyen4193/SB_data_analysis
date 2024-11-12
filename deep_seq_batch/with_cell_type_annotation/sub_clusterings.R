##### read lists of selected clusters for each case. 
source("/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/sub_clustering_indices.R")

##### read cell type annotation
source("/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/cell_type_annotation.R")

for (PROJECT in c(
  "SBharadwaj_20240318_Sample_1_4_7_8_2_5",
  "SBharadwaj_20240318_Sample_1_4_7_8"
)){
  groups <- names(sub_clusters[[PROJECT]])
  for (g in groups){
    tmp <- sub_clusters[[PROJECT]][[g]]
    tmp <- unlist(lapply(sub_clusters[[PROJECT]][[g]], function(x){
      all.annotations[[PROJECT]][[sprintf("%s", x)]]
    }))
    sub_clusters[[PROJECT]][[g]] <- tmp
  }
}