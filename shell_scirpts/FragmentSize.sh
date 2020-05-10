#!/bin/sh
set -e
WTH2AZ_IP1=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/WT1_H2AZ_repeat1_q20_s_rm.bam              
WTH3_IP1=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/WT1_H3_repeat1_q20_s_rm.bam                  
WT_in1=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/WT1_input_repeat1_q20_s_rm.bam                 
muH2AZ_IP1=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/ino80_H2AZ_repeat1_q20_s_rm.bam            
muH3_IP1=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/ino80_H3_repeat1_q20_s_rm.bam                
mu_in1=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/ino80_input_repeat1_q20_s_rm.bam               
WTH2AZ_IP2=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/WT1_H2AZ_repeat2_q20_s_rm.bam              
WTH3_IP2=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/WT1_H3_repeat2_q20_s_rm.bam                
WT_in2=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/WT1_input_repeat2_q20_s_rm.bam              
muH2AZ_IP2=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/ino80_H2AZ_repeat2_q20_s_rm.bam            
muH3_IP2=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/ino80_H3_repeat2_q20_s_rm.bam                
mu_in2=/mnt/data4/ylf/atino80/WL_ara_ino80_H2AZ/reChIPseq/Bam/ino80_input_repeat2_q20_s_rm.bam   
#H2AZ
bamPEFragmentSize  --bamfiles  ${WTH2AZ_IP1} ${WTH2AZ_IP2}  ${muH2AZ_IP1} ${muH2AZ_IP2}  \
--histogram H2AZ_bam_fragsize_2.pdf  --numberOfProcessors 4  \
--samplesLabel WTH2AZ_IP1 WTH2AZ_IP2  muH2AZ_IP1 muH2AZ_IP2 --maxFragmentLength 500  
echo "1done"
#H3
bamPEFragmentSize  --bamfiles  ${WTH3_IP1} ${WTH3_IP2}  ${muH3_IP1} ${muH3_IP2}  \
--histogram H3_bam_fragsize_2.pdf  --numberOfProcessors 4  \
--samplesLabel  WTH3_IP1 WTH3_IP2 muH3_IP1 muH3_IP2  --maxFragmentLength 500 
echo "2done"
#input
bamPEFragmentSize  --bamfiles  ${WT_in1} ${WT_in2}  ${mu_in1} ${mu_in2} \
--histogram input_bam_fragsize_2.pdf  --numberOfProcessors 4  \
--samplesLabel WT_in1 WT_in2 mu_in1 mu_in2   --maxFragmentLength 500 
echo "3done"       