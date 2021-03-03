options(java.parameters = "-Xmx8000m")

library(conflicted)
library(dplyr)
library(gplots)
library(ggplot2)
library(gridExtra)
library(limma)
library(maptools); gpclibPermit()
library(rjson)
library(sp)
library(stringr)
library(tidyr)
library(xlsx)

PlotColors = c("aliceblue", "antiquewhite2", "antiquewhite4", "blue", "black",
			   "azure3", "brown", "chartreuse", "burlywood4", "chocolate4",
			   "hotpink", "lightgoldenrodyellow", "limegreen", "orange3",
			   "gray72", "lightpink3", "gray17", "sienna1")

# I/O

ReadDesign = function(FileName) { return(read.table(FileName, header=T)) }
ReadTargets = function(Design) { return(read.maimages(Design$FileName, source="agilent", green.only=T, names=paste(1:length(Design$Name), Design$Name, sep="."), verbose=F)) }
SaveData = function(Data, FileName) { save(Data,file=FileName) }
LoadData = function(FileName) { return(load(file=FileName)) }

# NORMALIZE

NormalizeData = function(Data, averageDups=T) {
	
	Data = backgroundCorrect(Data, method="normexp", verbose=F)
	DataBefore = Data
	Data = normalizeBetweenArrays(Data, method="quantile")
	if (averageDups) { Data = avereps(Data, Data$genes[,"ProbeName"]) }
	DataAfter = Data
	par(mfrow = c(2, 1))
	plotDensities(DataBefore, col=PlotColors[1:ncol(Data$E)], legend="bottomright"); title("Densities [before]");
	plotDensities(DataAfter, col=PlotColors[1:ncol(Data$E)], legend="bottomright"); title("Densities [after]");
	par(mfrow = c(1, 1))
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

# COMPARISON

ComparisonDEG = function(Data, Design, Level, Contrasts, ResultXlsx) {
	DesignMatrix = model.matrix(~ 0 + Design[[Level]])
	colnames(DesignMatrix) = levels(Design[[Level]])
	FitGenotype = lmFit(Data, DesignMatrix)
	MadeContrasts = makeContrasts(contrasts=Contrasts, levels=DesignMatrix)
	ContrastsFit = contrasts.fit(FitGenotype, MadeContrasts)
	MakeeBayess(ContrastsFit, ResultXlsx)
}

# PLOTTING

PlotPCA = function(Data) {
	TransposedData = aperm(Data)
	pca = prcomp(TransposedData, scale=T, center=T)
	xlab = paste("PC1", as.character(pca$sdev[1] ^ 2 / sum(pca$sdev ^ 2)))
	ylab = paste("PC2", as.character(pca$sdev[2] ^ 2 / sum(pca$sdev ^ 2)))
	plot(pca$x, type="p", col="blue", xlab=xlab, ylab=ylab, pch=25, lwd=2)
	text(pca$x, labels=rownames(pca$x), pos=3, offset=0.9)
	title("PCA");
	plot.new()
	grid.table(as.data.frame(t(summary(pca)$importance)))
	title("PCA Summary Importance");
}

PlotCorrelations = function(Data, top=0.15) {
	#This code is based on example https://rstudio-pubs-static.s3.amazonaws.com/98999_50e28d4bc1324523899f9b27949ba4fd.html
	Vars = apply(Data, 1, sd)
	thrhold = quantile(Vars, 1 - top) 
	Data = Data[Vars > thrhold,]
	pearsonCorr = as.dist(1 - cor(Data))
	hC = hclust(pearsonCorr)
	plot(hC, hang=-1)
	heatmap.2(as.matrix(Data),
		  dendrogram="column",
		  trace="none",
		  col=greenred(100),
		  main=paste("Heatmap of top", as.character(top * 100), "% variative genes"),
		  symm=F,
		  symkey=F,
		  symbreaks=F,
		  scale="none")
}

AggreGather = function(Data, FUN, ColName) {
	Data = aggregate(Data[, -ncol(Data)], list(Data[, "Line"]), FUN, drop=FALSE)
	Data$Line = Data$Group.1
	Data = gather(Data, key="Group.1", value="value", -Line, -Group.1, na.rm=FALSE)
	colnames(Data)[colnames(Data) == "Group.1"] = "Probe"
	colnames(Data)[colnames(Data) == "value"] = ColName
	return(Data)
}

PlotGeneExpression = function(Data, Design, GeneName, namesColumn) {
	GenePos = which(Data$genes[[namesColumn]] == GeneName)
	TempData = t(2 ^ Data$E[GenePos, , drop=F])
	ColNames = colnames(TempData)
	TempData = cbind(as.data.frame(TempData), as.matrix(Design$Name))
	colnames(TempData) = c(ColNames, "Line")
	Exprs = AggreGather(TempData, function(x) base::mean(x, na.rm=TRUE), "LogExpression")
	SDErrs = AggreGather(TempData, function(x) stats::sd(x, na.rm=TRUE), "SD")
	Exprs = merge(Exprs, SDErrs, by=c("Line", "Probe"))
	print(ggplot(data=Exprs, mapping=aes(x=Line, y=LogExpression)) +
	    ggtitle(paste("Expression of", GeneName)) +
		geom_point(size=5) +
		geom_errorbar(aes(ymin=LogExpression-SD, ymax=LogExpression+SD), width=.1) +
		theme(panel.grid.major=element_line(colour="gray", linetype="dashed"), axis.text.x=element_text(angle=90, hjust=1), plot.title=element_text(size=24, face="bold")) +
		scale_y_continuous(trans='log2'))
}

# MAIN

aCGH.Analysis = function(Options, log=NULL) {
	Report = toJSON(Options, indent=6, method="C")
	cat(Report)
	Vectorize = function(Line) { return(unlist(lapply(Line, function(x) { if (typeof(x) == "list") { paste(unlist(x), collapse='; ') } else { as.character(x) } }), recursive=FALSE)) }
	if (!is.null(log)) {
		if (file.exists(log)) {
			write(Vectorize(Options), ncolumns=length(Options), sep="\t", file=log, append=TRUE)
		} else {
			write(names(Options), ncolumns=length(Options), sep="\t", file=log, append=FALSE)
			write(Vectorize(Options), ncolumns=length(Options), sep="\t", file=log, append=TRUE)
		}
	}
	PlotsPDF = file.path(Options$OutputDir, paste(Options$TimeStamp, "_plots.pdf"))
	ComparisonXLS = file.path(Options$OutputDir, paste(Options$TimeStamp, "_comparisons.xls"))
	Design = ReadDesign(Options$TargetsList)
	if (length(Options$Excluded) > 0) Design = Design[- unlist(Options$Excluded), ]
	Data = ReadTargets(Design)
	if (Options$Plots) pdf(PlotsPDF, width=8, height=11)
	Data = NormalizeData(Data, averageDups=T)
	print(Data)
	if (Options$Plots) {
		PlotPCA(Data$E)
		PlotCorrelations(Data$E)
		for (GN in unlist(Options$ExpressionGenes)) { PlotGeneExpression(Data, Design, GN, Options$GeneNameColumn) }
		dev.null = dev.off()
	}
	if (Options$Comparisons) {
		ComparisonWorkBook = createWorkbook(type="xls")
		ComparisonDEG(Data, Design, Options$Level, unlist(Options$Contrasts), ComparisonWorkBook)
		saveWorkbook(ComparisonWorkBook, ComparisonXLS)
	}
}
