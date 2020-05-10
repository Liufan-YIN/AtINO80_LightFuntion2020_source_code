#!/bin/sh
set -e
path=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/BigWig
computeMatrix scale-regions -S ${path}/WT1_H2AZ_merge_RPKM.bigwig ${path}/WT1_H3_merge_RPKM.bigwig ${path}/WT1_input_merge_RPKM.bigwig  \
-R ~/genome/tair10/ylf_genome/tair10_col6.bed  \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero  \
-out  WT1_matrix.mat.gz
computeMatrix reference-point  \
--referencePoint TSS  -b 1000 -a 2000  \
-R ~/genome/tair10/ylf_genome/tair10_col6.bed  \
-S ${path}/WT1_H2AZ_merge_RPKM.bigwig ${path}/WT1_H3_merge_RPKM.bigwig ${path}/WT1_input_merge_RPKM.bigwig  \
--binSize 10  --sortRegions keep --missingDataAsZero  \
-out  WT1_TSS_matrix.mat.gz
computeMatrix reference-point  \
--referencePoint TES  -b 2000 -a 1000  \
-R ~/genome/tair10/ylf_genome/tair10_col6.bed  \
-S ${path}/WT1_H2AZ_merge_RPKM.bigwig ${path}/WT1_H3_merge_RPKM.bigwig ${path}/WT1_input_merge_RPKM.bigwig  \
--binSize 10  --sortRegions keep --missingDataAsZero  \
-out  WT1_TES_matrix.mat.gz
computeMatrix scale-regions -S ${path}/ino80_H2AZ_merge_RPKM.bigwig ${path}/ino80_H3_merge_RPKM.bigwig ${path}/ino80_input_merge_RPKM.bigwig  \
-R ~/genome/tair10/ylf_genome/tair10_col6.bed  \
-b 1000 -a 1000 -m 3000 --binSize 10 --sortRegions keep --missingDataAsZero  \
-out ino80_matrix.mat.gz
computeMatrix reference-point  \
--referencePoint TSS  -b 1000 -a 2000  \
-R ~/genome/tair10/ylf_genome/tair10_col6.bed  \
-S ${path}/ino80_H2AZ_merge_RPKM.bigwig ${path}/ino80_H3_merge_RPKM.bigwig ${path}/ino80_input_merge_RPKM.bigwig  \
--binSize 10  --sortRegions keep --missingDataAsZero  \
-out  ino80_TSS_matrix.mat.gz
computeMatrix reference-point  \
--referencePoint TES  -b 2000 -a 1000  \
-R ~/genome/tair10/ylf_genome/tair10_col6.bed  \
-S ${path}/ino80_H2AZ_merge_RPKM.bigwig ${path}/ino80_H3_merge_RPKM.bigwig ${path}/ino80_input_merge_RPKM.bigwig  \
--binSize 10  --sortRegions keep --missingDataAsZero  \
-out  ino80_TES_matrix.mat.gz
gzip -d *gz
sed '1d'   WT1_matrix.mat     > WT1_matrix_nohead.mat 
sed '1d'   WT1_TSS_matrix.mat > WT1_TSS_matrix_nohead.mat 
sed '1d'   WT1_TES_matrix.mat > WT1_TES_matrix_nohead.mat 
sed '1d'   ino80_matrix.mat   > ino80_matrix_nohead.mat
sed '1d' ino80_TSS_matrix.mat > ino80_TSS_matrix_nohead.mat 
sed '1d' ino80_TES_matrix.mat > ino80_TES_matrix_nohead.mat 
echo "done"   