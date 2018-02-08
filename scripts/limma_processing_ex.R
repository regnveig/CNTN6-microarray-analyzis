library(limma)
setwd("~/Projects/CNTN6_microarray/") #Optionally - set working dir here
source("scripts/limma_processing.R") #Path to limma_processing.R file, accodring to working dir

#data_path = "/mnt/storage/home/vsfishman/Projects/CNTN6_microarray//"
#my_design = read.table(paste(data_path,"2017-12-targets.txt",sep=""),header=T)
my_design = read.table("data/2018-01-10-targets.txt",header=T)
#my_design = my_design[my_design$Condition %in% c("WT","DEL"),] #Optionally - remove GLIA, DUP, etc.

#pdf(paste("results/",format(Sys.Date(),"%Y-%m-%d"),"allDensPlots.txt",sep=""))
#my_data = my_loadAndNormalize_data(data_path,my_design,averageDups = T, doDensityPlot = T)
#dev.off()

pdf("results/DensityPlots.pdf")
my_data = my_loadAndNormalize_data(my_design,averageDups = T, doDensityPlot = T)
save(my_data,file="results/2016-02-05_my_data.Robj")

my_data = my_loadAndNormalize_data(my_design,
                                   infile = "results/2016-02-05_my_data.Robj")

##################Plot density##########################
#pdf("results/DensityPlotBeforeCorr.pdf")
#plotDensities(my_data)
#dev.off()

#pdf("results/DensityPlotAfterCorrAndDupMerged.pdf")
#plotDensities(my_data)
#dev.off()
#########################################################

#save(my_data,file = "results/2018-01-10-my_data_all.noProbDups.normexp.quantile.Robj")

#sink(paste("results/",format(Sys.Date(),"%Y-%m-%d"),"results/2018-01-07-DEG.txt",sep=""))
#by_genotype_comparisons(my_data)
#sink()

#sink("results/2018-01-07-DEG-full-report.txt")
#by_genotype_comparisons(my_data,printDEgenes=T)
#sink()

#sink("results/2018-01-10-DEG-full-report.txt")
library(xlsx)
wb<-createWorkbook(type="xlsx")
b=doDEG(my_data,"Condition","WT-DEL",exclude=c("13.iNAF1nor6"),deGenesWb=wb)
b=doDEG(my_data,"Name","iTAF3del17C98 - (iTAF3del21+iTAF3del17)/2",exclude=c("13.iNAF1nor6"),deGenesWb=wb)
saveWorkbook(wb, paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-DEG.xlsx",sep=""))

sink("results/2018-01-10-DEG-summary.txt")
b=doDEG(my_data,"Name","iTAF3del17C92 - (iTAF3del21+iTAF3del17)/2",exclude=c("13.iNAF1nor6"),printDEgenes=FALSE)
b=doDEG(my_data,"Condition","WT-DEL",exclude=c("13.iNAF1nor6"),printDEgenes=FALSE)
sink()

pdf(paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-logExpAll_PCAplotNoDupProbes.pdf",sep=""))
my_plotPCA(my_data$E)
dev.off()

pdf(paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-logExpAll_HistAndHmNoDupProbes.pdf",sep=""))
my_correlations(my_data$E,top=0.3)
dev.off()

pdf(paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-fit_genotype_NoDupProbesPlotPCA.pdf",sep=""))
my_plotPCA(fit_genotype$coefficients)
dev.off()

pdf(paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-fit_genotype_HistAndHmNoDupProbes.pdf",sep=""))
my_correlations(fit_genotype$coefficients, top = 0.3)
dev.off()