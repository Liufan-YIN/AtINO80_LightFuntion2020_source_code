#!/bin/sh                                                                                           
set -e                                                                                              
sh SICER.sh ChIPseq/BED SRR6412321_q20_sort_rmdup_chr.bed SRR6412323_q20_sort_rmdup_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh SICER.sh ChIPseq/BED SRR6412322_q20_sort_rmdup_chr.bed SRR6412324_q20_sort_rmdup_chr.bed . tair10 1 200 150 0.9 200 0.05  
#identify significantly enriched regions with FDR<=0.05 
awk '{if($8<=0.05) print $0}'  SRR6412321_q20_sort_rmdup_chr-W200-G200-islands-summary    > sig_SRR6412321.bed       
awk '{if($8<=0.05) print $0}'  SRR6412322_q20_sort_rmdup_chr-W200-G200-islands-summary    > sig_SRR6412322.bed 
#obtain files having seqnames without prefix "chr" for annotation
ls -1 sig*.bed |while read id
do
echo "${id%.bed}" 
awk -F "hr" '{print  $2}' ${id} >  ${id%.bed}_nochr.bed
done  
wc -l sig*_nochr.bed
echo "done"