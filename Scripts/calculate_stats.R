source("/home/mohit/Downloads/liblsb-0.2.2/R/utils.R")
source("/home/mohit/Downloads/liblsb-0.2.2/R/stats.R")
source("/home/mohit/Downloads/liblsb-0.2.2/R/aes.R")

args = commandArgs(trailingOnly=TRUE)

data.in <- read.table(args[1], header=TRUE)
#library(plyr)
#data.grouped <- ddply(data.in,c("nP","runid"),summarise,wtime=max(wtime))
data.grouped <- data.in
data.stats <- CalculateDataSummary(data.grouped, measurevar="wtime", groupvars=c("nP"))
write.table(data.stats, args[2], row.names=FALSE)
