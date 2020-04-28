#!/usr/bin/env Rscript
#MA_SCE_plot.R
#the following will plot the levels of multiallelism across SCE segments
#input 1 = strain id, input 2 = nod id
args = commandArgs(trailingOnly=TRUE)
input_path <- paste("ma/",args[1],"/",args[2],"/",args[2],".sce_bed2",sep="")
output_path <- paste("ma/",args[1],"/",args[2],"/",args[2],".ma_plot.pdf",sep="")
sce_boundaries_path <- paste("bin/",args[1],"_chr_prop",sep="")
sce_boundaries <- read.table(sce_boundaries_path,sep=" ",header=FALSE)
ma_sce <- read.table(input_path,sep=" ",header=FALSE)
pdf(output_path)
plot(ma_sce$V7,ma_sce$V6,xlim=c(0,1),ylim=c(0,0.5),cex=0,xlab="Chromosome",ylab="Prop. multiallelism",axes = FALSE, main=args[2]) +axis(1,at=(sce_boundaries$V4),labels=as.vector(sce_boundaries$V1),las=2) +axis(side = 2,las=2)+ box()
segments(sce_boundaries$V5,0,sce_boundaries$V5,1, lwd=0.5, col="gray")
segments(ma_sce$V8,ma_sce$V7,ma_sce$V9,ma_sce$V7, lwd=2, col="red")
dev.off()
