my_loadAndNormalize_data = function(my_design,infile = NULL, averageDups = T, doDensityPlot = F){
#Load and normalize data
  if (!is.null(infile)){
     print ("Loading data from file")
     load(file=infile)    
     return (my_data)
  }
  else{
    	my_data = read.maimages(my_design$FileName,source="agilent",green.only=T,
    							names=paste(1:length(my_design$Name),my_design$Name,sep="."))
    	my_data = backgroundCorrect(my_data,method="normexp")
    	if (doDensityPlot){
    		plotDensities(my_data)
    		title("Before")
    	}
    	my_data = normalizeBetweenArrays(my_data,method="quantile")
    	if (averageDups){
    		my_data = avereps(my_data,my_data$genes[,"ProbeUID"])
    	}
    	if (doDensityPlot){
    		plotDensities(my_data)
    		title("After")
    	}
    	return (my_data)
  }	
}


remove_glia = function (my_data,my_design){
  glia = (my_design$Condition == "GLIA")
  neg = max(my_data$E[my_data$genes$ControlType==-1,])
  glia_genes = apply(my_data$E[,c("14.mGlia","13.mGlia")],1,max)>neg
  
  glia_cols = !(colnames(my_data$E) %in% c("14.mGlia","13.mGlia"))
  
}

#################### By genotype comparisons ######################

do_eBayess = function(fit, deGenesWb = NULL, sheetNames = c()){
	eBayess = list()
	sheetInd = 0

	for (trend in c(F,T)) {
	  for (robust in c(F,T)){
	    name = paste("trend",as.character(trend),"robust",as.character(robust),sep=".")
      eBayess[[name]] = eBayes(fit,trend = trend, robust = robust)
      print (paste("Options trend =",trend,"robust=",robust))
      print(summary(decideTests(eBayess[[name]],adjust.method="BH",p.value=0.05)))
      if (!is.null(deGenesWb)) {
        ind = 0
        if (length(sheetNames) != length(colnames(eBayess[[name]]$contrasts))){
          sheetNames = colnames(eBayess[[name]]$contrasts)
          }
        for (coeff in colnames(eBayess[[name]]$contrasts)){
          sheetInd = sheetInd + 1
          ind = ind + 1
          sheetName = paste(as.character(sheetInd),sheetNames[ind],sep="-")
          print (sheetName)
          sheet <- createSheet(deGenesWb,sheetName)
          rows <-createRow(sheet,rowIndex=1)
          sheetTitle <-createCell(rows, colIndex=1)
          setCellValue(sheetTitle[[1,1]], coeff)
          rows <-createRow(sheet,rowIndex=2)
          sheetTitle <-createCell(rows, colIndex=1)
          setCellValue(sheetTitle[[1,1]], paste("Trend=",as.character(trend),
                                                "Robust=",as.character(robust)))
          addDataFrame(as.data.frame(
                      topTable(eBayess[[name]],number=Inf,p.value = 0.05,coef=coeff)
                        ),sheet,startRow=3,startColumn=1)
        }
      }
	  }
	}
	return(eBayess)
}

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
	#colnames(my_design_matrix) = c("(Intercept)",levels(my_design$Name))

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
	eBayess_by_genotype <<- do_eBayess(fit_genotype_pairs,printDEgenes)

	print ("Compare del with all of wt simultaneously")
	allNormVSDel <<- makeContrasts(iTAF1nor25+iTAF1nor36+iTAF2nor24-(iTAF3del17*3),
									levels = my_design_matrix)
									
	fit_WTvsDEL <<- contrasts.fit(fit_genotype,allNormVSDel)
	eBayess_WTvsDEL <<- do_eBayess(fit_WTvsDEL,printDEgenes)
	
	PairWizeWTvsWTContrasts = makeContrasts(c36vs24 = iTAF1nor36 - iTAF2nor24,
	                                           c36vs25 = iTAF1nor36 - iTAF1nor25,
	                                           c25vs24 = iTAF1nor25 - iTAF2nor24,
	                                           levels = my_design_matrix)
	fit_genotype_pairs_bwWT = contrasts.fit(fit_genotype,PairWizeWTvsWTContrasts)
	eBayess_bwWt <<- do_eBayess(fit_genotype_pairs_bwWT,printDEgenes)
}

#################### WT ws DEL ######################
wt_vs_del_comparisons = function(my_data,exclude=c(),printDEgenes=FALSE){
	if (length(exclude) > 0) {
		data_colnames = colnames(my_data)
		to_exclude = !(data_colnames %in% exclude)
		temp_data = my_data
		temp_data$targets = as.data.frame(my_data$targets[to_exclude,])
		rownames(temp_data$targets) = rownames(my_data$targets)[to_exclude]
		colnames(temp_data$targets) = colnames(my_data$targets)
		temp_data$E = temp_data$E[,to_exclude]

		temp_design = my_design[to_exclude,]
	}
	else {
		temp_data = my_data
		temp_design = my_design
	}

	condition = factor(temp_design$Condition,levels=c("WT","DEL"))
	designWTvsDEL <<- model.matrix(~condition)
	fit_WTvsDEL2 <<- lmFit(temp_data,designWTvsDEL)
	eBayess_WTvsDEL2 <<- do_eBayess(fit_WTvsDEL2,printDEgenes)
}

#pdf("Test.pdf")
#plotSA(fit, main="Probe-level")
#dev.off()

exclude_data = function(my_data, exclude){
  data_colnames = colnames(my_data)
  to_exclude = !(data_colnames %in% exclude)
  if (sum(to_exclude) < length(colnames(my_data))) {
    temp_data = my_data
    temp_data$targets = as.data.frame(my_data$targets[to_exclude,])
    rownames(temp_data$targets) = rownames(my_data$targets)[to_exclude]
    colnames(temp_data$targets) = colnames(my_data$targets)
    temp_data$E = temp_data$E[,to_exclude]
    temp_design = my_design[to_exclude,]
    for (i in names(temp_design)){
      temp_design[[i]] = factor(temp_design[[i]])
    }
  }
  else {
    temp_data <<- my_data
    temp_design = my_design
    for (i in names(temp_design)){
      temp_design[[i]] = factor(temp_design[[i]])
    }
  }
  res = list(temp_data=temp_data,temp_design=temp_design)
  return(res)
}

doDEG = function(my_data, cmpLevel, str_contrasts, exclude = c(), deGenesWb=NULL){
	excluded = exclude_data(my_data,exclude)
	temp_data = excluded$temp_data
	temp_design = excluded$temp_design
	
	print (paste("Running DEG for contrasts ",str_contrasts," on levels"))
	print (levels(temp_design[[cmpLevel]]))
	my_model_matrix = model.matrix(~0+temp_design[[cmpLevel]])
	colnames(my_model_matrix) = levels(temp_design[[cmpLevel]])
	temp_contrasts = makeContrasts(contrasts = str_contrasts, levels = my_model_matrix)
	fitted_model <<- lmFit(temp_data,my_model_matrix)
	fitted_contrasts <<- contrasts.fit(fitted_model,temp_contrasts)
	return(do_eBayess(fitted_contrasts,deGenesWb = deGenesWb))
}

my_plotPCA = function(some_data){

  data_transposed = aperm(some_data)
	pca <- prcomp(data_transposed, scale=T, center = T)
	xlab = paste("PC1",as.character(pca$sdev[1]^2/sum(pca$sdev^2)))
	ylab = paste("PC2",as.character(pca$sdev[2]^2/sum(pca$sdev^2)))
	print(xlab)
	plot(pca$x, type="n", xlab = xlab, ylab = ylab)
	     
	text(pca$x, rownames(pca$x), cex=0.5)
	print(summary(pca))
}

my_correlations = function(some_data, top = 0.15){
	#This code is based on example https://rstudio-pubs-static.s3.amazonaws.com/98999_50e28d4bc1324523899f9b27949ba4fd.html
	library(gplots)
	mydata = some_data

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

my_plotExpression = function(fit,geneName,save_path="results/"){
  library(tidyr)
  library(ggplot2)
  if (is.null(fit$s2.post)){
	b = eBayes(fit)
  }
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

my_plotExpression2 = function(some_data,geneName,save_path="results/",exclude=c()){
	library(tidyr)
	library(dplyr)
	library(ggplot2)
  excluded = exclude_data(some_data,exclude)
  some_data = excluded$temp_data
  temp_design = excluded$temp_design
	genePos = which(some_data$genes$GeneName %in% c(geneName))
	temp = 2^some_data$E[genePos,,drop=F]
	temp = t(temp)
	cnames = colnames(temp)
	temp = cbind(as.data.frame(temp),as.matrix(temp_design$Name))
	colnames(temp) = c(cnames,"Replica")
	library(dplyr)
	exprs = aggregate(temp[,-ncol(temp)],list(temp$Replica), mean)
	sderrs = aggregate(temp[,-ncol(temp)],list(temp$Replica), sd)
	exprs = gather(exprs,key="Group.1")
	colnames(exprs) = c("Line","Probe","LogExpression")
	sderrs = gather(sderrs,key="Group.1")
	colnames(sderrs) = c("Line","Probe","SD")
	exprs$SD = sderrs$SD
	
	pdf(paste(paste(save_path,geneName,sep="/"),"Log2.pdf",sep="."))
	print(ggplot(data=exprs,mapping = aes(x=Line,y=LogExpression))+
		geom_point()+
		geom_errorbar(aes(ymin=LogExpression-SD, ymax=LogExpression+SD),width=.1)+
		theme(panel.grid.major = element_line(colour = "gray",linetype = "dashed"),
		      axis.text.x = element_text(angle = 90, hjust = 1))+
		scale_y_continuous(trans='log2')
		)
	dev.off()

}

#https://support.bioconductor.org/p/70175/
#The effect sizes are contained in fit$coefficients.
#The standard errors can be obtained from
#SE <- sqrt(fit$s2.post) * fit$stdev.unscaled
#This gives a matrix containing standard errors for every coefficient and every gene.



