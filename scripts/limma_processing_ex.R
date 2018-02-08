library(limma)
library(xlsx)
setwd("~/Projects/CNTN6_microarray/") #Optionally - set working dir here
source("scripts/limma_processing.R") #Path to limma_processing.R file, accodring to working dir

my_design = read.table("data/2018-01-10-targets.txt",header=T)
#my_design = my_design[my_design$Condition %in% c("WT","DEL"),] #Optionally - remove GLIA, DUP, etc.

#pdf("results/DensityPlots.pdf")
#my_data = my_loadAndNormalize_data(my_design,averageDups = T, doDensityPlot = T)
#dev.off()
#save(my_data,file="results/2018-02-07_my_data.Robj")

#> colnames(my_data$E)
#[1] "1.iTAF1nor36"     "2.iTAF1nor36"     "3.iTAF1nor25"     "4.iTAF1nor25"
#[5] "5.iTAF2nor24"     "6.iTAF2nor24"     "7.iTAF3del17"     "8.iTAF3del17"
#[9] "9.iTAF3del21"     "10.iTAF3del21"    "11.iTAF3del17C98" "12.iTAF3del17C98"
#[13] "13.mGlia"         "14.mGlia"         "15.iTAF4dup22"    "16.iNAF1nor6"
#[17] "17.iTAF4dup14"    "18.iTAF4dup14"    "19.iTAF4dup22"    "20.iTAF4dup22"
#[21] "21.iTAF2nor5"

my_data = my_loadAndNormalize_data(my_design,
                                   infile = "results/2018-02-07_my_data.Robj") #load my_data
#my_data_no_glia = remove_glia(my_data)
#save(my_data_noGlia,file="results/2018-02-07_my_data_no_GLIA.Robj")
load(file="results/2018-02-07_my_data_no_GLIA.Robj") #load my_data_noGlia


wb<-createWorkbook(type="xlsx")
b=doDEG(my_data,"Condition","WT-DEL",exclude=c("13.mGlia","14.mGlia"),deGenesWb=wb)
b=doDEG(my_data,"Name","iTAF3del17C98 - (iTAF3del21+iTAF3del17)/2",exclude=c("13.mGlia","14.mGlia"),deGenesWb=wb)
saveWorkbook(wb, paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-DEGexSamplesGLIA.xlsx",sep=""))

wb<-createWorkbook(type="xlsx")
excluded_samples = c("13.mGlia","14.mGlia","21.iTAF2nor5","16.iNAF1nor6")
b=doDEG(my_data,"Condition","WT-DEL",exclude=excluded_samples,deGenesWb=wb)
b=doDEG(my_data,"Name","iTAF3del17C98 - (iTAF3del21+iTAF3del17)/2",exclude=excluded_samples,deGenesWb=wb)
b=doDEG(my_data,"Name","(iTAF1nor36+iTAF2nor24+iTAF1nor25)/3 - (iTAF3del21+iTAF3del17)/2",
        exclude=excluded_samples,deGenesWb=wb)
saveWorkbook(wb, paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-DEGexSamplesGLIA_5_6.xlsx",sep=""))

wb<-createWorkbook(type="xlsx")
#excluded_samples = c()
excluded_samples = c("13.mGlia","14.mGlia","21.iTAF2nor5","16.iNAF1nor6","15.iTAF4dup22",
                     "17.iTAF4dup14","18.iTAF4dup14","19.iTAF4dup22","20.iTAF4dup22")
b=doDEG(my_data,"Condition","WT-DEL",exclude=excluded_samples,lfc=3,deGenesWb=wb)
b=doDEG(my_data,"Name","(iTAF1nor36+iTAF2nor24+iTAF1nor25)/3 - (iTAF3del21+iTAF3del17)/2",
        exclude=excluded_samples,lfc=3,deGenesWb=wb)
b=doDEG(my_data,"Name","iTAF3del17C98 - (iTAF3del21+iTAF3del17)/2",exclude=excluded_samples,lfc=3,deGenesWb=wb)
saveWorkbook(wb, paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-DEGexSamplesGLIA_DUP_5_6_lfc_3.xlsx",sep=""))


#wb<-createWorkbook(type="xlsx")
excluded_samples = c("13.mGlia","14.mGlia","21.iTAF2nor5","16.iNAF1nor6","15.iTAF4dup22",
                     "17.iTAF4dup14","18.iTAF4dup14","19.iTAF4dup22","20.iTAF4dup22","11.iTAF3del17C98",
                     "12.iTAF3del17C98")
#If you are, as me, surprized that followin these comparisons give different results,
#read this post: https://support.bioconductor.org/p/26251/
b=doDEG(my_data,"Condition","WT-DEL",exclude=excluded_samples)
b=doDEG(my_data,"Name","(iTAF1nor36+iTAF2nor24+iTAF1nor25)/3 - (iTAF3del21+iTAF3del17)/2",
        exclude=excluded_samples)
#saveWorkbook(wb, paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-DEGexSamplesGLIA_DUP_5_6_lfc_3.xlsx",sep=""))


pdf(paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-logExpAll_PCAplotNoDupProbes.pdf",sep=""))
my_plotPCA(my_data$E)
dev.off()

pdf(paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-logExpAll_HistAndHmNoDupProbes.pdf",sep=""))
my_correlations(my_data$E,top=0.3)
dev.off()


##############All the same w/o GLIA genes##################
wb<-createWorkbook(type="xlsx")
b=doDEG(my_data_noGlia,"Condition","WT-DEL",exclude=c("13.mGlia","14.mGlia"),deGenesWb=wb)
b=doDEG(my_data_noGlia,"Name","iTAF3del17C98 - (iTAF3del21+iTAF3del17)/2",exclude=c("13.mGlia","14.mGlia"),deGenesWb=wb)
saveWorkbook(wb, paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-DEGexSamplesGLIAexGLIAGenes.xlsx",sep=""))

wb<-createWorkbook(type="xlsx")
excluded_samples = c("13.mGlia","14.mGlia","21.iTAF2nor5","16.iNAF1nor6")
b=doDEG(my_data_noGlia,"Condition","WT-DEL",exclude=excluded_samples,deGenesWb=wb)
b=doDEG(my_data_noGlia,"Name","iTAF3del17C98 - (iTAF3del21+iTAF3del17)/2",exclude=excluded_samples,deGenesWb=wb)
b=doDEG(my_data_noGlia,"Name","(iTAF1nor36+iTAF2nor24+iTAF1nor25)/3 - (iTAF3del21+iTAF3del17)/2",
        exclude=excluded_samples,deGenesWb=wb)
saveWorkbook(wb, paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-DEGexSamplesGLIA_5_6exGLIAGenes.xlsx",sep=""))

wb<-createWorkbook(type="xlsx")
excluded_samples = c("13.mGlia","14.mGlia","21.iTAF2nor5","16.iNAF1nor6","15.iTAF4dup22",
                     "17.iTAF4dup14","18.iTAF4dup14","19.iTAF4dup22","20.iTAF4dup22")
b=doDEG(my_data_noGlia,"Condition","WT-DEL",exclude=excluded_samples,lfc=3,deGenesWb=wb)
b=doDEG(my_data_noGlia,"Name","(iTAF1nor36+iTAF2nor24+iTAF1nor25)/3 - (iTAF3del21+iTAF3del17)/2",
        exclude=excluded_samples,lfc=3,deGenesWb=wb)
b=doDEG(my_data_noGlia,"Name","iTAF3del17C98 - (iTAF3del21+iTAF3del17)/2",exclude=excluded_samples,lfc=3,deGenesWb=wb)
saveWorkbook(wb, paste("results/",format(Sys.Date(),"%Y-%m-%d"),"-DEGexSamplesGLIA_DUP_5_6_lfc_3exGLIAGenes.xlsx",sep=""))
