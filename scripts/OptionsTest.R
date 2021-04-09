source("aCGH-analysis.R")

LogFile = "default.log.tsv"

Options = list(
	TimeStamp = as.character(as.hexmode(as.integer((as.numeric(Sys.time()) * 1000) %% .Machine$integer.max))),
	TargetsList = "/dev/datasets/FairWind/.cloud/core/CNTN6-microarray-analyzis/datasets/Tomsk_202101/Targets.list",
	OutputDir = "/dev/datasets/FairWind/.cloud/core/CNTN6-microarray-analyzis/datasets/Tomsk_202101/Results",
	ExcludeExpressedGenesIn = list(21, 22),
	
	Plots = FALSE,
	Excluded = list(),
	ExpressionGenes = list("A_33_P3393200", "A_33_P3369153"),
	GeneNameColumn = "ProbeName",
	
	Comparisons = TRUE,
	Level = "Name",
	Contrasts = list("(s14_8dupCNTN6+s14_10dupCNTN6)/2-(iTAF1nor36+iTAF1nor25+iTAF2nor24)/3", "(s13_26delCNTN6+s13_27delCNTN6)/2-(iTAF1nor36+iTAF1nor25+iTAF2nor24)/3"),
	AnnotationFile = "/dev/datasets/FairWind/.cloud/core/CNTN6-microarray-analyzis/datasets/Tomsk_202101/Agilent_G2600D_Annotation.txt"
)

# Excluded = list(1:8)
# Excluded = list(9:29)

aCGH.Analysis(Options, log=LogFile)
