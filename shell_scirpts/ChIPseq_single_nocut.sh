#!/bin/sh
set -e
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP127/SRP127400/SRR6412321/SRR6412321.sra
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP127/SRP127400/SRR6412322/SRR6412322.sra
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP127/SRP127400/SRR6412323/SRR6412323.sra
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP127/SRP127400/SRR6412324/SRR6412324.sra
echo "download done"
fastq-dump  SRR6412321.sra
fastq-dump  SRR6412322.sra
fastq-dump  SRR6412323.sra
fastq-dump  SRR6412324.sra
fastq-dump  SRR6412325.sra
fastq-dump  SRR6412326.sra
fastq-dump  SRR6412327.sra
fastq-dump  SRR6412328.sra
echo "unzip done"
mkdir QC
fastqc -o QC -f fastq  *.fastq
echo "qc done"
#quality is good and don't need to cut cutadapter

ls *fastq|awk -F "." '{print $1}'|while read id
do
echo "${id}"
/software/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -p 8 -x ~/genome/tair10/Bowtie2Index/genome -U  ${id}.fastq  -S  ${id}.sam
echo "${id} sam bowtie done"
samtools view -q 20 ${id}.sam -bS > ${id}_q20.bam
echo "${id} q20 done"
samtools sort ${id}_q20.bam  -o ${id}_q20_sort.bam
echo "${id} sort done"
samtools rmdup -s ${id}_q20_sort.bam   ${id}_q20_sort_rmdup.bam
echo "${id} rmdup done"
samtools index ${id}_q20_sort_rmdup.bam
echo "${id} index done"
bamCoverage -b ${id}_q20_sort_rmdup.bam -o ${id}_q20_sort_rmdup_RPKM.bigwig --binSize 10 --normalizeUsing RPKM
echo "bigwig done"
bedtools bamtobed -i ${id}_q20_sort_rmdup.bam > ${id}_q20_sort_rmdup.bed
echo "bed done"
awk '{print "chr"$0}'  ${id}_q20_sort_rmdup.bed > ${id}_q20_sort_rmdup_chr.bed 
echo "bed chr done"
done
echo "all done"
