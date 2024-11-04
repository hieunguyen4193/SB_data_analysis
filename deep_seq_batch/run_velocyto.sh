export PATH=/home/hieunguyen/samtools/bin/:$PATH
path_to_masked_regions_gtf="/media/hieunguyen/HD01/storage/REF/mm10_rmsk.gtf";
path_to_modified_gtf="/media/hieunguyen/HD01/storage/REF/build-mm10/mm10-2020-A_build/gencode.vM23.primary_assembly.annotation.gtf.filtered";
for sample in ctrl_gut_CD45_1 treated_gut_myeloid ctrl_gut_CD45_2 treated_liver_myeloid ctrl_gut_myeloid treated_gut_CD45_1 ctrl_liver_myeloid treated_gut_CD45_2;do \
samtools_threads=12;
path_to_cellranger_output=/media/hieunguyen/HD01/tmp/240318_Bharadwaj_Izcue_CellRanger_output/${sample};
velocyto run10x -m ${path_to_masked_regions_gtf} ${path_to_cellranger_output} ${path_to_modified_gtf} -@ ${samtools_threads};done
