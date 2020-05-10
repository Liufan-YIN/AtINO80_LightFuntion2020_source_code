VennPlot2 <- function(Path,Afilename,Bfilename,Asamplename,Bsamplename,mainname) {
  suppressMessages(library(VennDiagram))  
  suppressMessages(library("ggplot2"))
  #path-the cutrrent path,Afilename-target1 genes, Bfilename-target2 genes
  setwd(Path)
  Afile <- read.delim(Afilename,header=T)
  Bfile <- read.delim(Bfilename,header=T)
  #venn
  print("VennPlot")
  common <- length(intersect(Afile[,1],Bfile[,1])); A <- nrow(Afile); B <- nrow(Bfile); M <- matrix(c(common,A-common,B-common,33602-A-B+common),nrow=2);
  result <- fisher.test(M,alternative="greater");
  print(c("Target",A,"diff",B,"common",common))
  print(c("pvalue",result$p.value))
  print(result$estimate)
  venn.diagram(list(Asamplename=Afile[,1],Bsamplename=Bfile[,1]),resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5),fill=c("red","green"),filename = paste(mainname,"venn.tiff"))               
  }