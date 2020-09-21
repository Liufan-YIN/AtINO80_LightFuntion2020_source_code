VennPlot1 <- function(Afilename,Bfilename,Asamplename,Bsamplename,mainname) {
  #packages
  suppressMessages(library(VennDiagram))  
  suppressMessages(library("ggplot2"))
  suppressMessages(library(clusterProfiler))
  suppressMessages(library(org.At.tair.db))
  #Afilename-input file name for target genes; Bfilename-input file name for differentially regulated genes; Asamplename-sample name for Afile; Bsamplename-sample name for Bfile; mainname-prefix for output file
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
  #GO
  commonlist <-  intersect(Afile[,1],row.names(Bfile))
  write.table(data.frame(gene_id=commonlist),paste(mainname,"overlap.xls"),sep="\t",row.names=F)    
  ego <-  enrichGO(gene=commonlist,OrgDb=org.At.tair.db,keyType='TAIR',ont='BP',pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1);
  ego <- setReadable(ego,org.At.tair.db,keyType='TAIR');                                            
  write.table(data.frame(ego),paste(mainname,"overlap_BP.xls",sep=""),sep="\t",row.names=F)        
  }