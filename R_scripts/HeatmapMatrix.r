setwd("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/Matrix/merge/heatmap/figure_matrix")
name <- "WT"
Sample1 <- read.delim("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/Matrix/merge/heatmap/WT1_TSS_matrix_nohead.mat",header=F)
Sample2 <- read.delim("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/Matrix/merge/heatmap/WT1_TES_matrix_nohead.mat",header=F)
genebed <- Sample1[,1:6]
row.names(genebed) <- genebed[,4]
Sample1H2AZ <- Sample1[,7:306] - Sample1[,607:906];row.names(Sample1H2AZ) <- genebed[,4]
Sample2H2AZ <- Sample2[,7:306] - Sample2[,607:906];row.names(Sample2H2AZ) <- genebed[,4]
H2AZ <- cbind(genebed,Sample1H2AZ,Sample2H2AZ)
H2AZ <- H2AZ[which(rowSums(H2AZ[,7:306]==0) < 200),]
H2AZ <- H2AZ[order(rowSums(H2AZ[,107:306])),]
write.table(H2AZ,paste(name,"_H2AZ_TSS_TES.matrix.mat",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
WT_H2AZ <- H2AZ
Sample1H3 <- Sample1[,307:606]  - Sample1[,607:906];row.names(Sample1H3) <- genebed[,4]
Sample2H3 <- Sample2[,307:606]  - Sample2[,607:906];row.names(Sample2H3) <- genebed[,4]
H3 <- cbind(genebed,Sample1H3,Sample2H3)
H3 <- subset(H3,is.element(row.names(H3),row.names(WT_H2AZ))==T)
H3 <- H3[row.names(WT_H2AZ),]
write.table(H3,paste(name,"_H3_TSS_TES.matrix.mat",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
name <- "ino80"
Sample1 <- read.delim("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/Matrix/merge/heatmap/ino80_TSS_matrix_nohead.mat",header=F)
Sample2 <- read.delim("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/Matrix/merge/heatmap/ino80_TES_matrix_nohead.mat",header=F)
genebed <- Sample1[,1:6]
row.names(genebed) <- genebed[,4]
Sample1H2AZ <- Sample1[,7:306] - Sample1[,607:906];row.names(Sample1H2AZ) <- genebed[,4]
Sample2H2AZ <- Sample2[,7:306] - Sample2[,607:906];row.names(Sample2H2AZ) <- genebed[,4]
H2AZ <- cbind(genebed,Sample1H2AZ,Sample2H2AZ)
H2AZ <- subset(H2AZ,is.element(row.names(H2AZ),row.names(WT_H2AZ))==T)
H2AZ <- H2AZ[row.names(WT_H2AZ),]
write.table(H2AZ,paste(name,"_H2AZ_TSS_TES.matrix.mat",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
Sample1H3 <- Sample1[,307:606] - Sample1[,607:906];row.names(Sample1H3) <- genebed[,4]
Sample2H3 <- Sample2[,307:606] - Sample2[,607:906];row.names(Sample2H3) <- genebed[,4]
H3 <- cbind(genebed,Sample1H3,Sample2H3)
H3 <- subset(H3,is.element(row.names(H3),row.names(WT_H2AZ))==T)
H3 <- H3[row.names(WT_H2AZ),]
write.table(H3,paste(name,"_H3_TSS_TES.matrix.mat",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
nrow(H2AZ)
#[1] 33083