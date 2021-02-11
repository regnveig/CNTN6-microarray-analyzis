library(tidyr)
library(ggplot2)

setwd("~/Projects/CNTN6_microarray/results/")
load(file="aLL_results.Robj")

data_path = "/mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data/"
my_design = read.table(paste(data_path,"2017-12-targets.txt",sep=""),header=T)


pdf("Different_normalizations.pdf")
for (i in 1:length(all_results)) {
   for (j in 1:length(all_results[[i]])){
     all_results[[i]][[j]] = all_results[[i]][[j]] / all_results[[i]][[j]][1]
  }
  temp = as.data.frame(all_results[[i]])
  temp2 = gather(temp)
  cellLines = paste(1:8,my_design$Name,sep=".")
  temp2$CellLine = factor(rep(cellLines,ncol(temp)),levels = cellLines )
  #print (temp2$CellLine)

  iTAF1nor25_positions = which(my_design$Name %in% c("iTAF1nor25"))
  points = as.data.frame(list(x = iTAF1nor25_positions, y = rep(0.6353053119,length(iTAF1nor25_positions))))
  #print (points)
  print(
      ggplot(data=temp2, aes(x=CellLine, y=value, fill=key)) +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_point(data=points, aes(x=x, y=y, fill=NULL)) 
	)
}
dev.off()
