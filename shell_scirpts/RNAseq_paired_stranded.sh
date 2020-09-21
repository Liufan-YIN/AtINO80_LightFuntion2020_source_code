#!/bin/sh
set -e
gzip -d *.fastq.gz 
echo "unzip done"
mkdir QC 
fastqc -o QC -f fastq  *.fastq
echo "qc done"

ls -1 *R1*fastq | while read id
do
echo "${id%R1*}"
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  -u -20 -U -20 -m 80  -o trim_${id%R1*}R1.fastq -p trim_${id%R1*}R2.fastq  ${id%R1*}R1_001.fastq  ${id%R1*}R2_001.fastq
done
echo "trim done"
mkdir trimQC
fastqc -o trimQC -f fastq  trim*.fastq
echo "trim qc done"
echo "done"

ls -1 trim*_R1.fastq | while read id
do
echo "${id%_R1*}"
hisat2 -p 8 --rna-strandness RF -x ~/genome/tair10/hisat2_index/tair10_tran -1 ${id%_R1*}_R1.fastq  -2  ${id%_R1*}_R2.fastq  -S  ${id%_R1*}.sam
echo "mapping done"
samtools view -q 20  -bS ${id%_R1*}.sam   > ${id%_R1*}_q20.bam
echo "q20 bam done"
samtools sort ${id%_R1*}_q20.bam -o ${id%_R1*}_q20_s.bam
echo "sort done"
samtools index ${id%_R1*}_q20_s.bam
echo "index done"
bamCoverage -b ${id%_R1*}_q20_s.bam -o ${id%_R1*}.bigwig --binSize 10 --normalizeUsing RPKM
echo "bigwig done" 
featureCounts -T 8 -p -s 2 -t exon -g gene_id -F GTF  -a  ~/genome/tair10/genes.gtf -o  ${id%_R1*}.count  ${id%_R1*}_q20_s.bam
echo "count done"
sed '1d' ${id%_R1*}.count |  awk '{print $1"\t"$7}' > ${id%_R1*}.rawcount
done
echo "all done"

