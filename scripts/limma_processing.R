my_loadAndNormalize_data = function(data_path,my_design,infile = NULL,averageDups = T){
#Load and normalize data
  if (!is.null(infile)){
     print ("Loading data from file")
     load(file=infile)    
     return (my_data)
#    load(file=)
  }
  else{
	  my_data = read.maimages(my_design$FileName,path=data_path,source="agilent",green.only=T)
	  my_data = backgroundCorrect(my_data,method="normexp")
	  my_data = normalizeBetweenArrays(my_data,method="quantile")
	  if (averageDups){
	    my_data = avereps(my_data,my_data$genes[,"ProbeUID"])
	  }
	  return (my_data)
  }	
}

#################### By genotype comparisons ######################

by_genotype_comparisons = function(my_data, printDEgenes=F)
{
  print ("DGE comparison")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  
	my_design_matrix <<- model.matrix(~0+my_design$Name)

	#my_design_matrix
	#  my_design$NameiTAF1nor25 my_design$NameiTAF1nor36 my_design$NameiTAF2nor24
	#1                        0                        1                        0
	#2                        0                        1                        0
	#3                        1                        0                        0
	#4                        1                        0                        0
	#5                        0                        0                        1
	#6                        0                        0                        1
	#7                        0                        0                        0
	#8                        0                        0                        0
	#  my_design$NameiTAF3del17
	#1                        0
	#2                        0
	#3                        0
	#4                        0
	#5                        0
	#6                        0
	#7                        1
	#8                        1


	colnames(my_design_matrix) = levels(my_design$Name)

	#  iTAF1nor25 iTAF1nor36 iTAF2nor24 iTAF3del17
	#1          0          1          0          0
	#2          0          1          0          0
	#3          1          0          0          0
	#4          1          0          0          0
	#5          0          0          1          0
	#6          0          0          1          0
	#7          0          0          0          1
	#8          0          0          0          1
	#attr(,"assign")
	#[1] 1 1 1 1
	#attr(,"contrasts")
	#attr(,"contrasts")$`my_design$Name`
	#[1] "contr.treatment"


	PairWizeWTvsDelContrasts <<- makeContrasts(
  	        vs25 = iTAF1nor25 - iTAF3del17,	  
						vs36 = iTAF1nor36 - iTAF3del17,
						vs24 =  iTAF2nor24 - iTAF3del17,
						levels = my_design_matrix)

	fit_genotype <<- lmFit(my_data,my_design_matrix)

	fit_genotype_pairs <<- contrasts.fit(fit_genotype,PairWizeWTvsDelContrasts)
	eBayess <<- list()
	for (trend in c(F,T)) {
	  for (robust in c(F,T)){
	    name = paste("trend",as.character(trend),"robust",as.character(robust),sep=".")
      eBayess[[name]] <<- eBayes(fit_genotype_pairs,trend = trend, robust = robust)
      print (paste("Options trend =",trend,"robust=",robust))
      print(summary(decideTests(eBayess[[name]],adjust.method="BH",p.value=0.05)))
      if (printDEgenes) {
        for (coeff in colnames(eBayess[[name]]$contrasts)){
          print (paste("Comparison ",coeff))
          write.table (topTable(eBayess[[name]],number=Inf,p.value = 0.05,coef=coeff),sep="\t")
        }
      }
	  }
	}

	print ("Compare del with all of wt simultaneously")
	allNormVSDel <<- makeContrasts(iTAF1nor25+iTAF1nor36+iTAF2nor24-iTAF3del17*3,levels = my_design_matrix)
	fit_WTvsDEL <<- contrasts.fit(fit_genotype,allNormVSDel)
	eBayess_WTvsDEL <<- list()
	for (trend in c(F,T)) {
	  for (robust in c(F,T)){
	    name = paste("trend",as.character(trend),"robust",as.character(robust),sep=".")
	    eBayess_WTvsDEL[[name]] <<- eBayes(fit_WTvsDEL,trend = trend, robust = robust)
	    print (paste("Options trend =",trend,"robust=",robust))
	    print(summary(decideTests(eBayess_WTvsDEL[[name]],adjust.method="BH",p.value=0.05)))
	    if (printDEgenes) {
	      for (coeff in colnames(fit_WTvsDEL[[name]]$contrasts)){
	        print (paste("Comparison ",coeff))
	        write.table (topTable(fit_WTvsDEL[[name]],number=Inf,p.value = 0.05,coef=coeff),sep="\t")
	      }
	    }
	  }
	}
	
	PairWizeWTvsWTContrasts = makeContrasts(c36vs24 = iTAF1nor36 - iTAF2nor24,
	                                           c36vs25 = iTAF1nor36 - iTAF1nor25,
	                                           c25vs24 = iTAF1nor25 - iTAF2nor24,
	                                           levels = my_design_matrix)
	fit_genotype_pairs_bwWT = contrasts.fit(fit_genotype,PairWizeWTvsWTContrasts)
	eBayess_bwWt = list()
	for (trend in c(F,T)) {
	  for (robust in c(F,T)){
	    name = paste("trend",as.character(trend),"robust",as.character(robust),sep=".")
	    eBayess_bwWt[[name]] = eBayes(fit_genotype_pairs_bwWT,trend = trend, robust = robust)
	    print (paste("Options trend =",trend,"robust=",robust))
	    print(summary(decideTests(eBayess_bwWt[[name]],adjust.method="BH",p.value=0.05)))
#	    if (printDEgenes) {
#	      for (coeff in colnames(eBayess[[name]]$contrasts)){
#	        print (paste("Comparison ",coeff))
#	        write.table (topTable(eBayess[[name]],number=Inf,p.value = 0.05,coef=coeff),sep="\t")
#	      }
#	    }
	  }
	}
	
}

#################### WT ws DEL ######################
wt_vs_del_comparisons = function(my_data){
	condition = factor(my_design$Condition,levels=c("WT","DEL"))
	designWTvsDEL <<- model.matrix(~condition)
	fit_WTvsDEL2 <<- lmFit(my_data,designWTvsDEL)
	eBayes_WTvsDEL2_FF <<- eBayes(fit_WTvsDEL2)
	eBayes_WTvsDEL2_TF <<- eBayes(fit_WTvsDEL2, trend = T)
	eBayes_WTvsDEL2_TT <<- eBayes(fit_WTvsDEL2, trend = T, robust = T)

	summary(decideTests(eBayes_WTvsDEL2_FF,adjust.method="BH",p.value=0.05))
	#       (Intercept) conditionDEL
	#Down             0            1
	#NotSig          34        62975
	#Up           62942            0

	topTable(eBayes_WTvsDEL2_FF,adjust.method="BH",p.value=0.05)
	#Removing intercept from test coefficients
	#      Row Col Start
	#46707 285 131   828
	#                                                         Sequence ProbeUID
	#46707 GCTTCCAGGGCTTATAACTTAATTGGTACGGTGTCATTTAAGTTTCGCATTAGCCCATGA    43540
	#     ControlType     ProbeName  GeneName SystematicName
	#46707           0 A_33_P3310966 LINC01503       AK092192
	#                                                            Description
	#46707 gb|Homo sapiens cDNA FLJ34873 fis, clone NT2NE2014950. [AK092192]
	#          logFC  AveExpr        t      P.Value  adj.P.Val        B
	#46707 -3.954283 1.800412 -11.6055 6.966373e-07 0.04387143 -4.28405



	summary(decideTests(eBayes_WTvsDEL2_TF,adjust.method="BH",p.value=0.05))
	#       (Intercept) conditionDEL
	#Down             0            0
	#NotSig          32        62976
	#Up           62944            0


	summary(decideTests(eBayes_WTvsDEL2_TT,adjust.method="BH",p.value=0.05))
	#       (Intercept) conditionDEL
	#Down             0            0
	#NotSig          39        62976
	#Up           62937            0

}

#pdf("Test.pdf")
#plotSA(fit, main="Probe-level")
#dev.off()


my_plotPCA = function(fit_genotype){
	#This code is based on example https://rstudio-pubs-static.s3.amazonaws.com/98999_50e28d4bc1324523899f9b27949ba4fd.html
	data_transposed = aperm(fit_genotype$coefficients)
	pca <- prcomp(data_transposed, scale=T, center = T)
	xlab = paste("PC1",as.character(pca$sdev[1]^2/sum(pca$sdev^2)))
	ylab = paste("PC2",as.character(pca$sdev[2]^2/sum(pca$sdev^2)))
	print(xlab)
	plot(pca$x, type="n", xlab = xlab, ylab = ylab)
	     
	text(pca$x, rownames(pca$x), cex=0.5)
	print(summary(pca))
}

my_correlations = function(fit_genotype, top = 0.15){
	#This code is based on example https://rstudio-pubs-static.s3.amazonaws.com/98999_50e28d4bc1324523899f9b27949ba4fd.html
  library(gplots)
	mydata = fit_genotype$coefficients
	#remove 0 and negative values
	
	#lets use only top % if DEGs
	vars = apply(mydata,1,sd)
	thrhold = quantile(vars, 1-top) 
	mydata = mydata[vars>thrhold,]
	pearsonCorr <- as.dist(1 - cor(mydata))
	hC <- hclust(pearsonCorr)
	plot(hC)
	
	#now let's clip data range so we have better color-code
	
	
	heatmap.2(as.matrix(mydata),dendrogram = "column",trace="none",col=greenred(100),
	        main = paste("Heatmap of top",as.character(top*100),"% variative genes"),
	        symm=F,symkey=F,symbreaks=F, scale="none")
}

my_plotExpression = function(fit,geneName,save_path="results/",plotAll=T){
  library(tidyr)
  library(ggplot2)
  SE <- sqrt(fit$s2.post) * fit$stdev.unscaled  
  genePos = which(fit$genes$GeneName %in% c(geneName))

  if (length(genePos) == 1){
    plot_data = as.data.frame(as.matrix(fit$coefficients[genePos,]))
    colnames(plot_data) = c(genePos)
    plot_data$names = rownames(plot_data)
    SE = as.data.frame(as.matrix(SE[genePos,]))    
  }
  else{
    plot_data = as.data.frame(aperm(fit$coefficients[genePos,]))
    plot_data$names = rownames(plot_data)
  }

  plot_data = gather(plot_data,key=names)
  colnames(plot_data) = c("Line","Probe","LogExpression")
  SE = gather(as.data.frame(SE))
  plot_data$SE = SE$value
  pdf(paste(paste(save_path,geneName,sep="/"),"pdf",sep="."))
  print(ggplot(data=plot_data,mapping = aes(x=Line,y=LogExpression))+
    geom_point()+
    geom_errorbar(aes(ymin=LogExpression-SE, ymax=LogExpression+SE),width=.1)+
    theme(panel.grid.major = element_line(colour = "gray",linetype = "dashed"))
    )
  dev.off()
}

#https://support.bioconductor.org/p/70175/
#The effect sizes are contained in fit$coefficients.
#The standard errors can be obtained from
#SE <- sqrt(fit$s2.post) * fit$stdev.unscaled
#This gives a matrix containing standard errors for every coefficient and every gene.



