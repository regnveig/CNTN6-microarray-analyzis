library(limma)
setwd("~/Projects/CNTN6_microarray/")
source("scripts/limma_processing.R")

data_path = "/mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data/"
my_design = read.table(paste(data_path,"2017-12-targets.txt",sep=""),header=T)

#my_data = my_loadAndNormalize_data(data_path,my_design,averageDups = T)
my_data = my_loadAndNormalize_data(data_path,my_design,
            infile = "results/018-01-07-my_data.noProbDups.normexp.quantile.Robj")
#save(my_data,file = "results/018-01-07-my_data.noProbDups.normexp.quantile.Robj")
#file = "results/2017-01-06-my_data.normexp.quantile.Robj"

#sink("results/2018-01-07-DEG.txt")
by_genotype_comparisons(my_data)
#sink()

sink("results/2018-01-07-DEG-full-report.txt")
by_genotype_comparisons(my_data,printDEgenes=T)
sink()

pdf("results/2018-01-fit_genotype_PCAplotNoDupProbes.pdf")
my_plotPCA(fit_genotype)
dev.off()

pdf("results/2018-01-fit_genotype_PCAplotNoDupProbes.pdf")
my_plotPCA(fit_genotype)
dev.off()

pdf("results/2018-01-fit_genotype_HistAndHmapNoDupProbes.pdf")
my_correlations(fit_genotype, top = 0.10)
dev.off()


temp = list()
temp$coefficients = my_data$E
colnames(temp$coefficients) = paste(1:8,my_design$Name)

pdf("results/unlog_unnormed_noBackExpVals_PCA.pdf")
pdf("results/unlog_unnormed_BackCorrectExpVals_PCA.pdf")

pdf("results/log_normed_ExpVals_PCA.pdf")
my_plotPCA(temp)
dev.off()
pdf("results/log_normed_ExpVals_hClast.pdf")
my_correlations(temp, top=0.3)
dev.off()
