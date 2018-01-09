> library(limma)
> setwd("~/Projects/CNTN6_microarray/")
                          > my_design = read.table(paste(data_path,"2017-12-targets.txt",sep=""),header=T) 
                          > data_path = "/mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data/"    )
my_data = read.maimages(my_design$FileName,path=data_path,source="agilent",g 
                            
                            Read /mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data//SG11464135_257236314153_S001_GE1_1105_Oct12_1_1.txt 
                            Read /mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data//SG11464135_257236314153_S001_GE1_1105_Oct12_1_2.txt 
                            Read /mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data//SG11464135_257236314153_S001_GE1_1105_Oct12_1_3.txt 
                            Read /mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data//SG11464135_257236314153_S001_GE1_1105_Oct12_1_4.txt 
                            Read /mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data//SG11464135_257236314153_S001_GE1_1105_Oct12_2_1.txt 
                            Read /mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data//SG11464135_257236314153_S001_GE1_1105_Oct12_2_2.txt 
                            Read /mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data//SG11464135_257236314153_S001_GE1_1105_Oct12_2_3.txt 
                            Read /mnt/storage/home/vsfishman/Projects/CNTN6_microarray/data//SG11464135_257236314153_S001_GE1_1105_Oct12_2_4.txt 
                            > pdf("results/unlog_unnormed_BackCorrectExpVals_PCA.pdf")
                            > my_data = backgroundCorrect(my_data,method="normexp")
                            Array 1 corrected
                            Array 2 corrected
                            Array 3 corrected
                            Array 4 corrected
                            Array 5 corrected
                            Array 6 corrected
                            Array 7 corrected
                            Array 8 corrected
> temp=my_data
> colnames(temp$E) = my_design$Name
> pdf("results/plotDensitiesBeforeNorm.pdf")
> plotDensities(temp)
> dev.off()
null device 
1 
> colnames(temp$E) = paste(1:8,my_design$Name)
> pdf("results/plotDensitiesBeforeNorm.pdf")
> plotDensities(temp)
> dev.off()
null device 
1 
> temp.bw = normalizeBetweenArrays(my_data,method="quantile")
> colnames(temp.bw$E) = paste(1:8,my_design$Name)
> pdf("results/plotDensitiesAfterNorm.pdf")
> plotDensities(temp.bw)
> dev.off()
null device 
1 
