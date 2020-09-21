VennPlot2 <- function(Afilename,Bfilename,Asamplename,Bsamplename,mainname) {
  #packages
  suppressMessages(library(VennDiagram))  
  suppressMessages(library("ggplot2"))
  #Afilename-input file name for targetA genes; Bfilename-input file name for targetB genes; Asamplename-sample name for Afile; Bsamplename-sample name for Bfile; mainname-prefix for output file
  Afile <- read.delim(Afilename,header=T)
  Bfile <- read.delim(Bfilename,header=T)
  #venn
  print("VennPlot")
  common <- length(intersect(Afile[,1],Bfile[,1])); A <- nrow(Afile); B <- nrow(Bfile); M <- matrix(c(common,A-common,B-common,33602-A-B+common),nrow=2);
  result <- fisher.test(M,alternative="greater");
  print(c("TargetA",A,"TargetB",B,"common",common))
  print(c("pvalue",result$p.value))
  print(result$estimate)
  venn.diagram(list(Asamplename=Afile[,1],Bsamplename=Bfile[,1]),resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5),fill=c("red","green"),filename = paste(mainname,"venn.tiff"))            
  }