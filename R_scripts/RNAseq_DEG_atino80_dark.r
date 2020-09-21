library(Hmisc)                                                                                      
library(DESeq2)                                                                                     
library(pheatmap)                                                                                   
WT_1 <- read.delim("trim_D_WT_1.rawcount")
WT_2 <- read.delim("trim_D_WT_2.rawcount")
WT_3 <- read.delim("trim_D_WT_3.rawcount")
mu_1 <- read.delim("trim_D_in80_1.rawcount")
mu_2 <- read.delim("trim_D_in80_2.rawcount")
mu_3 <- read.delim("trim_D_in80_3.rawcount")
raw_count <- data.frame(gene_id=WT_1[,1],WT_1=WT_1[,2],WT_2=WT_2[,2],WT_3=WT_3[,2],mu_1=mu_1[,2],mu_2=mu_2[,2],mu_3=mu_3[,2]) 
count_data=raw_count[,2:ncol(raw_count)]                                                            
row.names(count_data)=raw_count[,1]                                                                 
colnames(count_data) = c("WT_1","WT_2","WT_3","mu_1","mu_2","mu_3")                                 
condition=factor(c("WT","WT","WT","mu","mu","mu"))                                                  
dds=DESeqDataSetFromMatrix(count_data,DataFrame(condition),design=~condition);
dds_filter=dds[rowSums(counts(dds)) > 1,]
dds_DES=DESeq(dds_filter);
samplelist=c("mu")
name="Darkino80"
for (n in 1:length(samplelist) ) 
{
sample=samplelist[n]
print(sample)
res=results(dds_DES,contrast=c("condition",sample,"WT"))
diff=subset(res,abs(res$log2FoldChange)>=log2(1.5) & res$pvalue <= 0.05)
up=subset(res,res$log2FoldChange >=log2(1.5) & res$pvalue <=0.05)
down=subset(res,res$log2FoldChange <= -log2(1.5) & res$pvalue <=0.05)
write.table(res,paste(name,"vsWT.xls",sep=""),sep="\t",row.names=T)
write.table(diff,paste(name,"diff.xls",sep="_"),sep="\t",row.names=T)
write.table(up,paste(name,"up.xls",sep="_"),sep="\t",row.names=T)
write.table(down,paste(name,"down.xls",sep="_"),sep="\t",row.names=T)
print(c(sample,nrow(res),nrow(diff),nrow(up),nrow(down)))
} 
                                                                                                  
genebed <- read.delim("tair10_col6.bed",header=F)
row.names(genebed) <- genebed[,4]
genebed <- genebed[as.character(row.names(data.frame(res))),]
gene<- GRanges(seqnames=genebed[,1],ranges=IRanges(start=genebed[,2],end=genebed[,3]),names=genebed[,4],strand=genebed[,6])
gene_sp <- split(gene, factor(elementMetadata(gene)$names,levels=as.character(row.names(data.frame(res)))))
rowRanges(dds_filter) <- gene_sp
FPKM <- data.frame(fpkm(dds_filter))
FPKM$WT_mean <- data.frame(rowMeans(FPKM[,1:3]))[,1]
FPKM$mu_mean <- data.frame(rowMeans(FPKM[,4:6]))[,1]
FPKM$log2FC <- data.frame(res)$log2FoldChange
FPKM$pvalue <- data.frame(res)$pvalue
FPKM$padj   <- data.frame(res)$padj
write.table(FPKM,"Darkino80_FPKM_FC_pvalue.xls",sep="\t")  