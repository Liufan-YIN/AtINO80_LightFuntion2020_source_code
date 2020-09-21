ProfilePlot <- function(colMeanDataA,colMeanDataB,DataAname,DataBname,mainname,filename,Ylim) {
#colMeanDataA-average density value for each bin in WT; colMeanDataB-average density value for each bin in mutant; DataAname-name of WT; DataBname-name of mutant; mainname-title of profile; filename-output file name,Ylim-ylim of profile
pdf(filename)
Dot <- c(seq(-1000,-10,length.out=100),seq(0,2990,length.out=300),seq(3000,3990,length.out=100))
Xlim <- c(-1000,4000)
sample1 <- colMeanDataA
sample2 <- colMeanDataB
lengendname <- c(DataAname,DataBname)
Color <- c("black","red")  
plot(x=Dot,y=sample1,tck=0.02,,bty="l",ty="l", lty=1, lwd=4,xaxt="n",col=Color[1],ylim=Ylim,xlim=Xlim,xlab='',ylab="occupancy",main=mainname)
lines(x=Dot,y=sample2,tck=0.02,bty="l",ty="l", lty=1, lwd=4,xaxt="n",col=Color[2],ylim=Ylim,xlim=Xlim,xlab='',ylab="occupancy")
legend("topright",lengendname,col=Color,lty=1,lwd=4,bty='n')
axis(side=c(1,4),at=c(-1000,0,3000,4000),
labels=c("-1kb","TSS","TTS","1kb"),font=2,tck=0.02)
axis(side=2,font=2,tck=0.02)
dev.off() 
print(t.test(sample1,sample2,alternative="greater")$p.value)    
print(t.test(sample1[101:400],sample2[101:400],alternative="greater")$p.value) 
}

############# read input files
WT <- read.delim("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/Matrix/merge/profile/WT1_matrix_nohead.mat",header=F)
mu <- read.delim("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/Matrix/merge/profile/ino80_matrix_nohead.mat",header=F)
row.names(WT) <- WT[,4]
row.names(mu) <- mu[,4]
genebed <- WT[,1:6]

############# plot H2AZ profile
WT_IP <-   WT[,7:506]- WT[,1007:1506];
mu_IP <-   mu[,7:506]- mu[,1007:1506];  
WT_IP_enrich <- WT_IP[rowSums(WT_IP)>0,]
mu_IP_enrich <- mu_IP[rowSums(mu_IP)>0,]
genelist <- unique(union(row.names(WT_IP_enrich),row.names(mu_IP_enrich)))
length(genelist)
#[1] 30376
WT_IP <- subset(WT_IP,is.element(row.names(WT_IP),genelist)==T) 
mu_IP <- subset(mu_IP,is.element(row.names(mu_IP),genelist)==T) 
ProfilePlot(colMeans(WT_IP)/10,colMeans(mu_IP)/10,"WT","atino80","@H2A.Z","H2AZ_IP-input.pdf",c(0,8))
#[1] 1.12731e-22
#[1] 1.78391e-28
Target <- read.delim("HY5_ino80_overlap.xls") 
WT_target <- subset(WT_IP,is.element(row.names(WT_IP),row.names(Target))==T)  
mu_target <- subset(mu_IP,is.element(row.names(mu_IP),row.names(Target))==T)
nrow(WT_target)#[1] 808
WT_target <- colMeans(WT_target)/10
mu_target <- colMeans(mu_target)/10
ProfilePlot(WT_target,mu_target,"WT","atino80","@H2AZ in overlap","H2AZ_in_overlap.pdf",c(-2,10))
#[1] 1.41215e-24
#[1] 3.15699e-39
Target <- read.delim("HY5_ino80_overlap_up.xls") 
WT_target <- subset(WT_IP,is.element(row.names(WT_IP),row.names(Target))==T)  
mu_target <- subset(mu_IP,is.element(row.names(mu_IP),row.names(Target))==T)
nrow(WT_target)#[1] 514
WT_target <- colMeans(WT_target)/10
mu_target <- colMeans(mu_target)/10
ProfilePlot(WT_target,mu_target,"WT","atino80","@H2AZ in overlap_up","H2AZ_in_overlap_up.pdf",c(-2,10))
#[1] 1.042385e-29
#[1] 1.092678e-47
Target <- read.delim("/Volumes/OpheliaData/atino80/atino80_0317/HY5_regulate_genes/HY5_INO80_1/HY5_ino80_overlap_down.xls") 
WT_target <- subset(WT_IP,is.element(row.names(WT_IP),row.names(Target))==T)  
mu_target <- subset(mu_IP,is.element(row.names(mu_IP),row.names(Target))==T)
nrow(WT_target)#[1] 294
WT_target <- colMeans(WT_target)/10
mu_target <- colMeans(mu_target)/10
ProfilePlot("/Volumes/OpheliaData/atino80/atino80_0317/ChIPseq2/Matrix/merge/profile/HY5_ino80_1",WT_target,mu_target,"WT","atino80","@H2AZ in overlap_down","H2AZ_in_overlap_down.pdf",c(-2,10))
#[1] 1.182877e-13
#1] 7.011148e-19

############# plot H3 profile
WT_IP <-   WT[,507:1006]- WT[,1007:1506];
mu_IP <-   mu[,507:1006]- mu[,1007:1506];  
WT_IP_enrich <- WT_IP[rowSums(WT_IP)>0,]
mu_IP_enrich <- mu_IP[rowSums(mu_IP)>0,]
genelist <- unique(union(row.names(WT_IP_enrich),row.names(mu_IP_enrich)))
length(genelist)
#[1] 32670
WT_IP <- subset(WT_IP,is.element(row.names(WT_IP),genelist)==T) 
mu_IP <- subset(mu_IP,is.element(row.names(mu_IP),genelist)==T) 
ProfilePlot(colMeans(WT_IP)/10,colMeans(mu_IP)/10,"WT","atino80","@H3","H3_IP-input.pdf",c(-2,6))
#[1] 5.220671e-27
#[1] 1.208631e-73
Target <- read.delim("HY5_ino80_overlap.xls") 
WT_target <- subset(WT_IP,is.element(row.names(WT_IP),row.names(Target))==T)  
mu_target <- subset(mu_IP,is.element(row.names(mu_IP),row.names(Target))==T)
nrow(WT_target)#[1] 826
WT_target <- colMeans(WT_target)/10
mu_target <- colMeans(mu_target)/10
ProfilePlot(WT_target,mu_target,"WT","atino80","@H3 in overlap","H3_in_overlap.pdf",c(-2,10))
#[1] 8.609018e-21
#[1] 1.339471e-53
Target <- read.delim("HY5_ino80_overlap_up.xls") 
WT_target <- subset(WT_IP,is.element(row.names(WT_IP),row.names(Target))==T)  
mu_target <- subset(mu_IP,is.element(row.names(mu_IP),row.names(Target))==T)
nrow(WT_target)#[1] 530
WT_target <- colMeans(WT_target)/10
mu_target <- colMeans(mu_target)/10
ProfilePlot(WT_target,mu_target,"WT","atino80","@H3 in overlap_up","H3_in_overlap_up.pdf",c(-2,10))
#[1] 4.653906e-21
#[1] 4.278609e-47
Target <- read.delim("HY5_ino80_overlap_down.xls") 
WT_target <- subset(WT_IP,is.element(row.names(WT_IP),row.names(Target))==T)  
mu_target <- subset(mu_IP,is.element(row.names(mu_IP),row.names(Target))==T)
nrow(WT_target)#[1] 296
WT_target <- colMeans(WT_target)/10
mu_target <- colMeans(mu_target)/10
ProfilePlot(WT_target,mu_target,"WT","atino80","@H3 in overlap_down","H3_in_overlap_down.pdf",c(-2,10))
#[1] 6.593643e-19
#[1] 1.629341e-67

############# plot H2AZ/H3 profile
Target <- read.delim("HY5_ino80_overlap.xls") 
WT_target <- subset(WT,is.element(row.names(WT),row.names(Target))==T)  
mu_target <- subset(mu,is.element(row.names(mu),row.names(Target))==T)
nrow(WT_target)#[1] 826
WT_target <- 10*log2(colMeans(WT_target[,7:506])/colMeans(WT_target[,507:1006]))
mu_target <- 10*log2(colMeans(mu_target[,7:506])/colMeans(mu_target[,507:1006]))
ProfilePlot(WT_target,mu_target,"WT","atino80","@H2AZ/H3 in overlap","H2AZvsH3_in_overlap.pdf",c(-2,10))
#[1] 7.697697e-05
#[1] 0.0005787423
Target <- read.delim("HY5_ino80_overlap_up.xls") 
WT_target <- subset(WT,is.element(row.names(WT),row.names(Target))==T)  
mu_target <- subset(mu,is.element(row.names(mu),row.names(Target))==T)
nrow(WT_target)#[1] 530
WT_target <- 10*log2(colMeans(WT_target[,7:506])/colMeans(WT_target[,507:1006]))
mu_target <- 10*log2(colMeans(mu_target[,7:506])/colMeans(mu_target[,507:1006]))
ProfilePlot(WT_target,mu_target,"WT","atino80","@H2AZ/H3 in overlap_up","H2AZvsH3_in_overlap_up.pdf",c(-2,10))
#[1] 5.59668e-08
#[1] 2.307303e-07
Target <- read.delim("HY5_ino80_overlap_down.xls") 
WT_target <- subset(WT,is.element(row.names(WT),row.names(Target))==T)  
mu_target <- subset(mu,is.element(row.names(mu),row.names(Target))==T)
nrow(WT_target)#[1] 296
WT_target <- 10*log2(colMeans(WT_target[,7:506])/colMeans(WT_target[,507:1006]))
mu_target <- 10*log2(colMeans(mu_target[,7:506])/colMeans(mu_target[,507:1006]))
ProfilePlot(WT_target,mu_target,"WT","atino80","@H2AZ/H3 in overlap_down","H2AZvsH3_in_overlap_down.pdf",c(-2,10))
#[1] 0.1701653
#[1] 0.4878595

