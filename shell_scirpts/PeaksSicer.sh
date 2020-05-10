#!/bin/sh                                                                                           
set -e                                                                                              
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed    WT1_H2AZ_repeat1_q20_s_rm_chr.bed    WT1_input_repeat1_q20_s_rm_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed      WT1_H3_repeat1_q20_s_rm_chr.bed    WT1_input_repeat1_q20_s_rm_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed  ino80_H2AZ_repeat1_q20_s_rm_chr.bed  ino80_input_repeat1_q20_s_rm_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed    ino80_H3_repeat1_q20_s_rm_chr.bed  ino80_input_repeat1_q20_s_rm_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed    WT1_H2AZ_repeat2_q20_s_rm_chr.bed    WT1_input_repeat2_q20_s_rm_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed      WT1_H3_repeat2_q20_s_rm_chr.bed    WT1_input_repeat2_q20_s_rm_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed  ino80_H2AZ_repeat2_q20_s_rm_chr.bed  ino80_input_repeat2_q20_s_rm_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed    ino80_H3_repeat2_q20_s_rm_chr.bed  ino80_input_repeat2_q20_s_rm_chr.bed . tair10 1 200 150 0.9 200 0.05 
echo "done"  
awk '{if($8<=0.05) print $0}'    WT1_H2AZ_repeat1_q20_s_rm_chr-W200-G200-islands-summary    > sig_WT1_H2AZ_repeat1.bed       
awk '{if($8<=0.05) print $0}'  ino80_H2AZ_repeat1_q20_s_rm_chr-W200-G200-islands-summary    > sig_ino80_H2AZ_repeat1.bed 
awk '{if($8<=0.05) print $0}'    WT1_H2AZ_repeat2_q20_s_rm_chr-W200-G200-islands-summary    > sig_WT1_H2AZ_repeat2.bed       
awk '{if($8<=0.05) print $0}'  ino80_H2AZ_repeat2_q20_s_rm_chr-W200-G200-islands-summary    > sig_ino80_H2AZ_repeat2.bed 
awk '{if($8<=0.05) print $0}'      WT1_H3_repeat1_q20_s_rm_chr-W200-G200-islands-summary    > sig_WT1_H3_repeat1.bed 
awk '{if($8<=0.05) print $0}'    ino80_H3_repeat1_q20_s_rm_chr-W200-G200-islands-summary    > sig_ino80_H3_repeat1.bed 
awk '{if($8<=0.05) print $0}'      WT1_H3_repeat2_q20_s_rm_chr-W200-G200-islands-summary    > sig_WT1_H3_repeat2.bed 
awk '{if($8<=0.05) print $0}'    ino80_H3_repeat2_q20_s_rm_chr-W200-G200-islands-summary    > sig_ino80_H3_repeat2.bed 
ls -1 sig*.bed |while read id
do
echo "${id%.bed}" 
awk -F "hr" '{print  $2}' ${id} >  ${id%.bed}_nochr.bed
done  
wc -l sig*_nochr.bed
echo "done"