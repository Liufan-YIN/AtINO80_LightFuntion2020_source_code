Plot <- function(Path,Afilename,Bfilename,Asamplename,Bsamplename,mainname) {
  suppressMessages(library(VennDiagram))  
  suppressMessages(library("ggplot2"))
  suppressMessages(library(clusterProfiler))
  suppressMessages(library(org.At.tair.db))
  #path-the cutrrent path,Afilename-target genes, Bfilename-diffregulated genes
  setwd(Path)
  Afile <- read.delim(Afilename,header=T)
  Bfile <- read.delim(Bfilename,header=T)
  #venn
  print("VennPlot")
  common <- length(intersect(Afile[,1],row.names(Bfile))); A <- nrow(Afile); B <- nrow(Bfile); M <- matrix(c(common,A-common,B-common,33602-A-B+common),nrow=2);
  result <- fisher.test(M,alternative="greater");
  print(c("Target",A,"diff",B,"common",common))
  print(c("pvalue",result$p.value))
  print(result$estimate)
  venn.diagram(list(Asamplename=Afile[,1],Bsamplename=row.names(Bfile)),resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5),fill=c("red","green"),filename = paste(mainname,"venn.tiff")) 
  #scatterplot
  print("scatterplot")
  SetA   <- Bfile
  Target <- Afile
  exp <- subset(SetA,is.element(row.names(SetA),Target[,1])==T)[,c(2,5)]
  loc_up <- intersect(which(exp$log2FoldChange>=log2(1.5)),which(exp$pvalue<=0.05))
  loc_down <- intersect(which(exp$log2FoldChange <= -log2(1.5)),which(exp$pvalue<=0.05))
  sig <- rep("no",times=nrow(exp))
  sig[loc_up] <- "up"             
  sig[loc_down] <- "down"  
  sig <- factor(sig,levels=c("up","down")) 
  exp$type <- sig
  exp$log10FDR <- -log(exp$pvalue,10)
  exp$log10FDR2 <- exp$log10FDR
  exp$log10FDR2[exp$log10FDR2 > 50] <- 50
  print(c("sig",length(sig),"up",length(loc_up),"down",length(loc_down))); 
  pdf("volcanoplot_overlap.pdf",width=9,height=6)
  yline <- c(-log2(1.5),log2(1.5))
  xline <- -log(0.05,10)   
  p <- qplot(x=exp$log10FDR2,y=exp$log2FoldChange,xlab="Expression Score: -log10[pvalue]",ylab="log2[expression fold-change]",size=I(1),colour=sig,shape=sig) + xlim(0,50)+ ylim(-10,10) + scale_color_manual(values=c("up"="red","down"="green")) + scale_shape_manual(values=c("up"=19,"down"=19)) + geom_hline(yintercept=yline,lty=2,size=I(0.5),colour="grey")+geom_vline(xintercept=xline,lty=2,size=I(0.5),colour="grey") + theme_bw() + theme(panel.background=element_rect(colour="black",size=1,fill="white"),panel.grid=element_blank())
  print(p)
  dev.off()  
  exp_up <- subset(exp,exp$type=="up")
  exp_down <- subset(exp,exp$type=="down")
  write.table(exp,paste(Asamplename,Bsamplename,"overlap.xls",sep="_"),sep="\t")
  write.table(exp_up,paste(Asamplename,Bsamplename,"overlap_up.xls",sep="_"),sep="\t")
  write.table(exp_down,paste(Asamplename,Bsamplename,"overlap_down.xls",sep="_"),sep="\t")
  #GO
  print("GO")
  setwd(paste(Path,"GO",sep=""))
  ego <- enrichGO(gene=row.names(exp_up),OrgDb=org.At.tair.db,keyType='TAIR',ont='BP',pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1);
  ego <- setReadable(ego,org.At.tair.db,keyType='TAIR');
  write.table(data.frame(ego),"overlap_up_GO_BP.xls",sep="\t",row.names=F) 
  ego <- enrichGO(gene=row.names(exp_down),OrgDb=org.At.tair.db,keyType='TAIR',ont='BP',pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1);
  ego <- setReadable(ego,org.At.tair.db,keyType='TAIR');
  write.table(data.frame(ego),"overlap_down_GO_BP.xls",sep="\t",row.names=F)    
}
