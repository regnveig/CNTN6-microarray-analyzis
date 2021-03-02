source("aCGH-analysis.R")

LogFile = "default.log.tsv"

Options = list(
	TimeStamp = as.character(as.hexmode(as.integer((as.numeric(Sys.time()) * 1000) %% .Machine$integer.max))),
	TargetsList = "/dev/datasets/FairWind/.cloud/core/CNTN6-microarray-analyzis/datasets/Tomsk_202101/Targets.list",
	OutputDir = "/dev/datasets/FairWind/.cloud/core/CNTN6-microarray-analyzis/datasets/Tomsk_202101/Results",
	Plots = TRUE,
	Comparisons = TRUE,
	Level = "Condition",
	Contrasts = list("DEL-DUP"),
	Excluded = list()
)

aCGH.Analysis(Options, log=LogFile)
