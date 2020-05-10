DiffPeak1 <- function(Path,samplename,peakfile,Fold,pvalue) {
setwd(Path)
Peak <- read.delim(peakfile,header=T)
Peaksummary_down <- Peak[,c(1:3,11,13)]                              
Peaksummary_up <- Peak[,c(1:3,8,10)] 
colnames(Peaksummary_down)=c("seqnames","start","end","FC","FDR")                        
colnames(Peaksummary_up)=c("seqnames","start","end","FC","FDR")                                              
rownames(Peaksummary_down)=paste(rep("diffPeak", nrow(Peaksummary_down)), 1:nrow(Peaksummary_down), sep="_") 
rownames(Peaksummary_up)=paste(rep("diffPeak", nrow(Peaksummary_up)), 1:nrow(Peaksummary_up), sep="_")      
sigPeakdown <- subset(Peaksummary_down,Peaksummary_down$FC>=Fold & Peaksummary_down$FDR <=pvalue)
sigPeakup <- subset(Peaksummary_up,Peaksummary_up$FC>=Fold & Peaksummary_up$FDR <=pvalue)
print(c(nrow(sigPeakdown),nrow(sigPeakup))) 
write.table(sigPeakdown,paste(samplename,"_sigPeakdown.xls",sep=""),sep="\t",row.names=F,quote=F)
write.table(sigPeakup,paste(samplename,"_sigPeakup.xls",sep=""),sep="\t",row.names=F,quote=F)
}
