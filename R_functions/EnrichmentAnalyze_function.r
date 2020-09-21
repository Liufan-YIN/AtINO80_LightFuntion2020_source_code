EnrichAna <- function(Afilename,Bfilename,Asamplename,Bsamplename,mainname) {
  suppressMessages(library(VennDiagram)) 
  #Afilename-input file name for target genes; Bfilename-input file name for differentially regulated genes; Asamplename-sample name for Afile; Bsamplename-sample name for Bfile; mainname-prefix for output file  
  Afile <- read.delim(Afilename,header=T)
  Bfile <- read.delim(Bfilename,header=T)
  common <- length(intersect(Afile[,1],row.names(Bfile))); A <- nrow(Afile); B <- nrow(Bfile); M <- matrix(c(common,A-common,B-common,33602-A-B+common),nrow=2);
  result <- fisher.test(M,alternative="greater");
  venn.diagram(list(Asamplename=Afile[,1],Bsamplename=row.names(Bfile)),resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5),fill=c("red","green"),filename = paste(mainname,"venn.tiff")) 
  print(result$p.value)
  print(result$estimate)
  print(common/A)
  genebed <- read.delim("tair10_col6.bed",header=F)
  Random <- sample(genebed[,4],B)
  common <- length(intersect(Afile[,1],Random)); A <- nrow(Afile); B <- length(Random); M <- matrix(c(common,A-common,B-common,33602-A-B+common),nrow=2); 
  print(common/A)  
}