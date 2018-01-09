#!/mnt/storage/home/vsfishman/Distr/R-3.4.1/bin/Rscript
cat("\nHi\n")
library(limma)

data_path = "/mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data/"
results_path = "/mnt/storage/home/vsfishman/Projects/CNTN6_microarray/results/"
my_design = read.table(paste(data_path,"2017-12-targets.txt",sep=""),header=T)

all_columns = list(list(G="gProcessedSignal",Gb="gBGUsed"),
		list(G="gMeanSignal",Gb="gBGMeanSignal"),
		list(G="gMedianSignal",Gb="gBGMedianSignal"))
all_bg_corrections = list("none", "subtract", "half", "minimum","normexp")
all_bwArr_normalizations = list("none", "scale", "quantile","cyclicloess")

count = 0

all_results = list()
for (column in all_columns) {
	test_data = read.maimages(my_design$FileName,path=data_path,source="agilent",green.only=T,
			 columns=column)
	column_cond_name = paste(column$G, column$Gb,sep=".")
	all_results[[column_cond_name]] = list()
	for (bg_correction in all_bg_corrections) {
		if (column$G == "gProcessedSignal" & bg_correction != "none") {next}
		test_data_bgCor = backgroundCorrect(test_data,method=bg_correction)
		for (normMethod in all_bwArr_normalizations) {
			test_data_norm = normalizeBetweenArrays(test_data_bgCor,method=normMethod)
			cond_name = paste(column$G, column$Gb, bg_correction,normMethod,sep=".")
			print(cond_name)
			all_results[[column_cond_name]][[cond_name]] = as.vector(test_data_norm[test_data_norm$gene$GeneName=="CNTN6",]$E)
		}
	}
}

print (all_results)
save(all_results,file="temp.Robj")
save(all_results,file=paste(results_path,"all_results.Robj",sep="/")
cat("\nDone\n")
