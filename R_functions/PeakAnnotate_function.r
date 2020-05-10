peakanno=function(filename,samplename) {
  suppressMessages(library("ChIPpeakAnno"));
  suppressMessages(library("TxDb.Athaliana.BioMart.plantsmart28"));
  txdb=TxDb.Athaliana.BioMart.plantsmart28;
  annDatagene=genes(txdb);
  peaks <- read.delim(filename,header=T)
  peaks_ranges <- GRanges(seqnames=peaks[,1],ranges=IRanges(start=peaks[,2],end=peaks[,3],names=paste("peak",rep(1:nrow(peaks)),sep="")),strand="*",seqinfo=seqinfo(txdb))
  annoOvgene <<- annotatePeakInBatch(peaks_ranges, AnnotationData=annDatagene,output="overlapping",maxgap=0L)
  write.table(annoOvgene, paste(samplename,"_annoOvgene.xls",sep=""), sep="\t",row.names=FALSE)
  annogene <- data.frame(genename=unique(as.character(subset(annoOvgene,annoOvgene$feature != "NA")$feature)))
  write.table(annogene, paste(samplename,"_annoOvgene_genelist.xls",sep=""), sep="\t",row.names=FALSE)
  print(nrow(annogene))  
}