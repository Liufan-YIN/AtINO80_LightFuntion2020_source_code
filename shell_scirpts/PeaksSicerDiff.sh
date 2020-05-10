#!/bin/sh                                                                                           
set -e                                                                                              
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed      WT1_H2AZ_merge_chr.bed      WT1_input_merge_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed        WT1_H3_merge_chr.bed      WT1_input_merge_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed    ino80_H2AZ_merge_chr.bed    ino80_input_merge_chr.bed . tair10 1 200 150 0.9 200 0.05 
sh /software/SICER_V1.1fuben/SICER/SICER.sh /mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bed    ino80_H3_merge_chr.bed      ino80_input_merge_chr.bed . tair10 1 200 150 0.9 200 0.05 
echo " sicer done"
python /software/SICER_V1.1fuben/SICER/src/find_union_islands.py -s tair10 -a ./WT1_H2AZ_repeats_overlap_peaks.bed -b ./ino80_H2AZ_repeats_overlap_peaks.bed  -o merge_H2AZ.bed 
python /software/SICER_V1.1fuben/SICER/src/compare_two_libraries_on_islands.py -s tair10 -a  ino80_H2AZ_merge_chr-1-removed.bed  -b  WT1_H2AZ_merge_chr-1-removed.bed -d  ./merge_H2AZ.bed  -f 150 -o ./FCH2AZ_merge.summary                                                                  
python /software/SICER_V1.1fuben/SICER/src/find_union_islands.py -s tair10 -a ./WT1_H3_repeats_overlap_peaks.bed   -b ./ino80_H3_repeats_overlap_peaks.bed    -o merge_H3.bed               
python /software/SICER_V1.1fuben/SICER/src/compare_two_libraries_on_islands.py -s tair10 -a  ino80_H3_merge_chr-1-removed.bed    -b  WT1_H3_merge_chr-1-removed.bed   -d  ./merge_H3.bed    -f 150 -o ./FCH3_merge.summary     
echo "sicer diff done"
