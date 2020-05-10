#!/bin/sh 
set -e    
ls -1 *repeat1*.bam | while read id
do
echo "${id%_repeat*}"   
samtools merge  ${id%_repeat*}_merge.bam  ${id%_repeat*}_repeat1_q20_s_rm.bam ${id%_repeat*}_repeat2_q20_s_rm.bam
echo "merge done"
samtools index ${id%_repeat*}_merge.bam
echo "index done"
bamCoverage -p 4 -b ${id%_repeat*}_merge.bam -o ${id%_repeat*}_merge_RPKM.bigwig --binSize 10 --normalizeUsing RPKM
echo "bigwig done"
bedtools bamtobed -i ${id%_repeat*}_merge.bam > ${id%_repeat*}_merge.bed
echo "bed done"
awk '{print "chr"$0}'  ${id%_repeat*}_merge.bed > ${id%_repeat*}_merge_chr.bed ;
echo "chr done"
done
echo "all done"