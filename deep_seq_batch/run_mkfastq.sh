path_to_cellranger="/home/hieunguyen/CRC1382/src_2023/src_pipeline/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH
sample_sheet="/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch/240318_Bharadwaj_Izcue_MolMed_scRNAseq.modified.csv";
inputdir="/home/hieunguyen/CRC1382/outdir/SBharadwaj_20240318/BCL"
cellranger mkfastq --run ${inputdir} --sample-sheet ${sample_sheet};


