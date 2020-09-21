BubblePlot <- function(GOfile,SizeLimits,SizeRange,BarValue) {
#package
suppressMessages(library("ggplot2"))
#GOfile-input file of GO analysis results including term order, function terms, gene count number and p value; SizeLimits-dot size limits for  gene number; Barvalue-color bar for log10(p value)
allGO=read.table(GOfile,sep="\t",header=TRUE)
GOup <- allGO[,c(1,3,4)]
GOup$names <- "Up"
GOdown <- allGO[,c(1,5,6)]
GOdown$names <- "Down"
colnames(GOup) <- c("Order","Count","Pvalue","names")
colnames(GOdown) <- c("Order","Count","Pvalue","names")
allGO=rbind(GOdown,GOup)
pdf("GO_bubble.pdf",width=4,height=9) 
p = ggplot(allGO,aes(factor(names,levels=c("Up","Down")),factor(allGO[,1]) ))
p=p + geom_point(aes(size=Count,colour=-1*log10(Pvalue)))+ scale_size("Gene Number",limits=SizeLimits,range=SizeRange) + scale_colour_gradient2(low ="green", mid="darkblue", high ="red",midpoint=0,breaks=BarValue)
print(p)
dev.off()
}