#!/usr/bin/env Rscript
##plot_driver_times.R
args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(ggplot2)

input_path <- paste("ma/",args[1],"/",args[1],".timing",sep="")
driver_times <- data.table(input_path,sep="\t",header=FALSE)
colnames(driver_times) <- c("nod","ma_prop","first_gen","driver","braf","egfr","Hras","Kras")

data <- data.table(DRIVER=character(), V2=numeric(), V3=numeric(), V4=numeric())
for (DRIVER in c("braf","egfr","hras","kras")){
 A <- NROW(driver_times[driver_times[[DRIVER]] == 1 & driver_times$first_gen == 1,])
 B <- NROW(driver_times[driver_times[[DRIVER]] == 1 & driver_times$first_gen == 0,])
 C <- NROW(driver_times[driver_times[[DRIVER]] == 0 & driver_times$first_gen == 1,])
 D <- NROW(driver_times[driver_times[[DRIVER]] == 0 & driver_times$first_gen == 0,])
 x <- fisher.test(matrix(cbind(A,B,C,D),nrow=2))
 data <- rbind(data, data.table(DRIVER,log2(c(x$conf.int[1])),log2(c(x$conf.int[2])),log2(x$estimate)))
}

colnames(data) <- c("driver","lower","upper","OR")
data$driver <- factor(data$driver, levels=rev(data$driver))
output_path <- paste("ma/",args[1],"/",args[1],".timing_plot.pdf",sep="")
pdf(output_path)
ggplot(data=data, aes(x=driver, y=OR, ymin=lower, ymax=upper)) + geom_pointrange() + geom_hline(yintercept=1, lty=2) + coord_flip() + xlab("Driver") + ylab("Log2/Odds") + theme_bw() 
dev.off()
