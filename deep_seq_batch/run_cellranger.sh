path_to_cellranger="/home/hieunguyen/CRC1382/src_2023/src_pipeline/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH

path_to_fastq="/home/hieunguyen/CRC1382/outdir/H252MAFX7";
path_to_outputdir="/home/hieunguyen/CRC1382/outdir/SBharadwaj_output";
mkdir -p $path_to_outputdir;

for sample in ctrl_gut_CD45  ctrl_gut_myeloid  treated_gut_CD45  treated_gut_myeloid;do \
cellranger count --id=SB_${sample} \
                   --transcriptome=/home/hieunguyen/CRC1382/storage/build-mm10/mm10 \
                   --fastqs=${path_to_fastq}/${sample} \
                   --sample=${sample} \
                   --localcores=20;

done