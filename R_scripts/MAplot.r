library("ggplot2")
FC <- read.delim("H2AZ_H3_FC.xls") 
MA_data <- data.frame(Avalue=FC[,6],Mvalue= -FC[,9],pvalue=FC[,10])
MA_data <- MA_data[order(MA_data$pvalue,decreasing=TRUE),]
row.names(MA_data) <- rep(1:nrow(MA_data))
Colour <- rep("grey",times=nrow(MA_data))
Colour[which(MA_data$pvalue <=0.05 & MA_data$Mvalue > 0)] <- "red"
Colour[which(MA_data$pvalue <=0.05 & MA_data$Mvalue < 0)] <- "dark green"
MA_data$Colour <- Colour
pdf(file="MAplot.pdf",width=9,height=6)
p = ggplot(MA_data,aes(Avalue,Mvalue))+xlim(0,12)+ylim(-2,2)
p=p + geom_point(shape=20,size=0.3,colour=Colour)+theme_bw()+theme(panel.background=element_rect(colour="black",size=1,fill="white"),panel.grid=element_blank())
yline2 <- 0    
p <- p + geom_hline(yintercept=yline2,lty=1,size=I(0.5),colour="blue")  
print(p)   
dev.off()