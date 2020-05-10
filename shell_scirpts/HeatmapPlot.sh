#add the header below to each matrix firstly
#@{"verbose":false,"scale":1,"skip zeros":false,"nan after end":false,"sort using":"mean","unscaled 5 prime":[0,0],"body":[0,0],"sample_labels":["TSS","TES"],"downstream":[2000,1000],"unscaled 3 prime":[0,0],"group_labels":["genes"],"bin size":[10,10],"upstream":[1000,2000],"group_boundaries":[0,33083],"sample_boundaries":[0,300,600],"max threshold":null,"ref point":["TSS","TES"],"min threshold":null,"sort regions":"keep","proc number":1,"bin avg type":"mean","missing data as zero":true}

#!/bin/sh
set -e
gzip *mat
plotHeatmap -m    WT_H2AZ_TSS_TES.matrix.mat.gz -out    WT_H2AZ_TSS_TES_order.pdf --sortRegions keep --colorList black,red -min 0 -max 40
plotHeatmap -m ino80_H2AZ_TSS_TES.matrix.mat.gz -out ino80_H2AZ_TSS_TES_order.pdf --sortRegions keep --colorList black,red -min 0 -max 40
plotHeatmap -m      WT_H3_TSS_TES.matrix.mat.gz -out      WT_H3_TSS_TES_order.pdf --sortRegions keep --colorList black,red -min 0 -max 40
plotHeatmap -m   ino80_H3_TSS_TES.matrix.mat.gz -out   ino80_H3_TSS_TES_order.pdf --sortRegions keep --colorList black,red -min 0 -max 40
echo "done"