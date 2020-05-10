BubblePlot <- function(Path,GOfile,SizeLimits,SizeRange,BarValue) {
suppressMessages(library("ggplot2"))
setwd(Path)
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