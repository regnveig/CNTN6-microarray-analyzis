library(limma)
library(xlsx)
library(stringr)
library(tidyr)

PlotColors = c("aliceblue", "antiquewhite2", "antiquewhite4", "blue", "black",
			   "azure3", "brown", "chartreuse", "burlywood4", "chocolate4",
			   "hotpink", "lightgoldenrodyellow", "limegreen", "orange3",
			   "gray72", "lightpink3", "gray17", "sienna1")

ReadDesign = function(FileName) { return(read.table(FileName, header=T)) }
ReadTargets = function(Design) { return(read.maimages(Design$FileName, source="agilent", green.only=T, names=paste(1:length(Design$Name), Design$Name, sep="."))) }
SaveData = function(Data, FileName) { save(Data,file=FileName) }
LoadData = function(FileName) { return(load(file=FileName)) }

# NORMALIZE

NormalizeData = function(Data, averageDups=T) {
	
	Data = backgroundCorrect(Data, method="normexp")
	plotDensities(Data, col=PlotColors[1:ncol(Data$E)], legend="bottomright"); title("Densities [before]");
	Data = normalizeBetweenArrays(Data, method="quantile")
	if (averageDups) { Data = avereps(Data, Data$genes[,"ProbeName"]) }
	plotDensities(Data, col=PlotColors[1:ncol(Data$E)], legend="bottomright"); title("Densities [after]");
	return(Data)
}

# BAYESS

BT = function(Bool) { if (Bool) { return("T") } else { return("F") } }

MakeeBayess = function(Fit, ResultXlsx, lfc=1) {
	eBayess = list()
	sheetInd = 0
	Set = 0
	SummarySheet = createSheet(ResultXlsx, "Summary")
	for (trend in c(F,T)) {
		for (robust in c(F,T)) {
			ind = 0
			Set = Set + 1
			name = paste("Trend=", BT(trend), ", Robust=", BT(robust), sep="")
			eBayess[[name]] = eBayes(Fit, trend=trend, robust=robust)
			sheetNames = colnames(eBayess[[name]]$contrasts)
			for (coeff in sheetNames) {
				sheetInd = sheetInd + 1
				ind = ind + 1
				sheetName = paste("(", as.character(sheetInd), ") ", coeff, ", t=", BT(trend), ", r=", BT(robust), sep="")
				sheet = createSheet(ResultXlsx,sheetName)
				Samples = str_split(c(coeff), "-")
				sheetTitle = c(
					paste("# Samples: ", paste(as.character(Samples[[1]]), collapse=", ")),
					paste("# Trend: ", as.character(trend)),
					paste("# Robust: ", as.character(robust)),
					"# DATA")
				addDataFrame(data.frame(matrix(unlist(sheetTitle), nrow=length(sheetTitle), byrow=TRUE)), sheet, startRow=1, startColumn=1, col.names=F, row.names=F)
				addDataFrame(data.frame(matrix(unlist(c(paste("# SUMMARY: ", name, sep=""))), nrow=1, byrow=TRUE)), SummarySheet, startRow=(Set - 1) * 6 + 1, startColumn=1, col.names=F, row.names=F)
				addDataFrame(spread(as.data.frame(summary(decideTests(eBayess[[name]], adjust.method="BH", p.value=0.05))), key=Var2, value=Freq), SummarySheet, startRow=(Set - 1) * 6 + 2, startColumn=1, row.names=F)
				addDataFrame(as.data.frame(topTable(eBayess[[name]], number=Inf, p.value=0.05, coef=coeff, lfc=lfc)), sheet, startRow=5, startColumn=1)
			}
		}
	}
}

# COMPARISON BY GT

ComparisonByGenotype = function(Data, Design, Contrasts, ResultXlsx) {
	DesignMatrix = model.matrix(~ 0 + Design$Name)
	colnames(DesignMatrix) = levels(Design$Name)
	FitGenotype = lmFit(Data, DesignMatrix)
	MadeContrasts = makeContrasts(contrasts=Contrasts, levels=DesignMatrix)
	ContrastsFit = contrasts.fit(FitGenotype, MadeContrasts)
	MakeeBayess(ContrastsFit, ResultXlsx)
}

# ----------------------------

# help(make.names)
# q()
Design = ReadDesign("/dev/datasets/FairWind/.cloud/core/CNTN6-microarray-analyzis/datasets/Tomsk_202101/Targets.list")
Data = ReadTargets(Design)

pdf("/dev/datasets/FairWind/.cloud/core/CNTN6-microarray-analyzis/datasets/Tomsk_202101/re.pdf")

Data = NormalizeData(Data)

wb<-createWorkbook(type="xlsx")
ComparisonByGenotype(Data, Design, c("s13_26-s14_10", "s13_27-s14_10", "s13_26-s13_27"), ResultXlsx=wb)
saveWorkbook(wb, "DEGexSamplesGLIA.xlsx")
dev.off()
