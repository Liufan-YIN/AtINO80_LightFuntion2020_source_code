#/bin/sh
set -e
gzip -d *gz
echo "jieya done"
mkdir QC 
fastqc -o QC -f fastq  *.fastq
echo "qc done"
#quality is good and don't need to cut cutadapter

ls -1 *R1.fastq |while read id
do
echo "${id}"
/software/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -p 8 -x ~/genome/tair10/Bowtie2Index/genome -1 ${id%_R1*}_R1.fastq -2  ${id%_R1*}_R2.fastq   -S  ${id%_R1*}.sam
echo "map done"
samtools view -q 20 ${id%_R1*}.sam -bS > ${id%_R1*}_q20.bam
echo "unique done"
samtools sort ${id%_R1*}_q20.bam  -o ${id%_R1*}_q20_s.bam
echo "sort done"
samtools rmdup ${id%_R1*}_q20_s.bam   ${id%_R1*}_q20_s_rm.bam
echo "rmdup done"
samtools index ${id%_R1*}_q20_s_rm.bam
echo "index done"
samtools view -F 4 -c {id%_R1*}_q20_s.bam
samtools view -F 4 -c ${id%_R1*}_q20_s_rm.bam
echo "reads done"
bamCoverage -b ${id%_R1*}_q20_s_rm.bam -o ${id%_R1*}_q20_s_rm_RPKM.bigwig --binSize 10 --normalizeUsing RPKM
echo "bigwig done"
bedtools bamtobed -i ${id%_R1*}_q20_s_rm.bam > ${id%_R1*}_q20_s_rm.bed
echo "bed done"
awk '{print "chr"$0}'  ${id%_R1*}_q20_s_rm.bed > ${id%_R1*}_q20_s_rm_chr.bed 
echo "bed chr done"
done
echo "all done"