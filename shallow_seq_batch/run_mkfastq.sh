path_to_cellranger="/home/hieunguyen/CRC1382/src_2023/src_pipeline/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH
sample_sheet="/home/hieunguyen/CRC1382/src_2023/SBharadwaj/SampleSheet.csv";
inputdir="/home/hieunguyen/CRC1382/outdir/tmp"
cellranger mkfastq --run ${inputdir} --sample-sheet ${sample_sheet};

