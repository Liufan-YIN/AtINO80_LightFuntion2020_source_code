Peakover <- function(Afilename,Bfilename,peak1name,peak2name,mainname){
 #packages
 suppressMessages(library("ChIPpeakAnno"))
 suppressMessages(library("GenomicRanges"))
 suppressMessages(library(BSgenome))
 suppressMessages(library(regioneR))
 suppressMessages(library(BSgenome.Athaliana.TAIR.TAIR9))
 #Afilename-filename of Afile; Bfilename-filename of Bfile; peak1name-name of peaks in Afile; peak2name-name of peaks in Bfile; mainname-prefiex of name in outputfile
 a <- data.frame(getGenome("TAIR9"))
 TAIR10 <- GRanges(seqnames=c(1,2,3,4,5,"Mt","Pt"),ranges=IRanges(start=a[,2],end=a[,3]),strand="*")
 Afile <- read.delim(Afilename,header=F)
 Bfile <- read.delim(Bfilename,header=T)
 Agranges <- GRanges(seqnames=Afile[,1],ranges=IRanges(start=Afile[,2],end=Afile[,3],names=paste("A",rep(1:nrow(Afile)),sep="")))
 Bgranges <- GRanges(seqnames=Bfile[,1],ranges=IRanges(start=Bfile[,2],end=Bfile[,3],names=paste("B",rep(1:nrow(Bfile)),sep="")))
 result <<- findOverlapsOfPeaks(Agranges,Bgranges,connectedPeaks="merge")
 pdf(paste(mainname,".pdf",sep=""))
 makeVennDiagram(result,NameOfPeaks=c(peak1name,peak2name),height=3000,width=3000,col="transparent",fill=c("red","green"),alpha=c(0.5,0.5),main=mainname)
 dev.off()
 write.table(data.frame(result$mergedPeaks)[,1:3],paste(mainname,".bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F) 
 permTest(A=Agranges,B=Bgranges,ntimes=1000,randomize.function=randomizeRegions,evaluate.function=numOverlaps,count.once=TRUE,genome=TAIR10,mc.set.seed=FALSE,mc.core=1) 
 }