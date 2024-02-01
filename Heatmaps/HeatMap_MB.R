#HeatMap plotting in R
#Version 1.0 adapted by Martha Bhattacharya, last updated 5-2-2020

#first time: install gplots (NOT ggplot!)
install.packages("gplots")
install.packages("RColorBrewer")

#load libraries

library(ggplot2)
library(plyr)
library(ggrepel)
library(calibrate)
library(reshape2)
library(gridExtra)
library(stats)
library(dplyr)
library(readr)
library(gplots)
library(RColorBrewer)

################Loading and labeling rows/columns of your samples########################

#set your working directory for your files - EDIT THIS
setwd("M:/Martha/Programming/R/Volcanos and Heatmaps")

#load data frame that has first col geneID, next few are normcounts, last is adjP with cutoff  - EDIT THIS
e13fulldata <- read.csv("M:/PAPER ASSEMBLY/Itch paper/Timecourse RNAseq/E13 normcounts for heatmap AdjP01.csv")
e13fulldata <- as.data.frame(e13fulldata)

#check the column names of this data frame
colnames(e13fulldata)

#rename the column names so that the list of genes is called GeneID
#you'll need to change the numbers in brackets if you have more than n=3 x 2 samples in your analysis
newcolnames <- c("GeneID", colnames(e13fulldata[2:8]))
colnames(e13fulldata) <- newcolnames

#make the row names the gene names
newrownames <- e13fulldata$GeneID
rownames(e13fulldata) <- newrownames

######## This section is an attempt to make colored labels for each genotype, not working yet and not required for program to run#####

#make the subsets of samples and geneid
sampleid <- colnames(e13fulldata[2:7])
geneid <- as.character(e13fulldata[,1])
head(geneid)
str(geneid)

#create a vector of genotypes corresponding to the column names - only needed for 
genotype <- c("WT", "WT","WT","Mut","Mut","Mut")

#put this together with the sampleID so there is a code for what the sample is
gAnnotate <- cbind(sampleid,genotype)
gAnnotate2 <- as.character(gAnnotate)
#gAnnotate

#map samples to a color for the upper legend by making a function 
#NOT WORKING YET 
mapgeno2color<-function(annotations){
  colorsVector = ifelse(annotations["genotype"]=="WT", 
                        "red", ifelse(annotations["genotype"]=="Mut", 
                                       "blue", "black"))
  return(colorsVector)
}
mapgeno2color(gAnnotate2)
gAnnotate
gAnnotate2

#################### Getting data in the right formats/log/normalized/scaled ############################

#need to add 1 and then take the log2of all the normcounts
#the numbers here assume 3 samples of each genotype
normplus1 <- e13fulldata[,2:7] + 1
normlog <- log2(normplus1)
normlog2 <- as.matrix(normlog)
e13fulldata2 <- cbind(e13fulldata[1],normlog, e13fulldata[8])

#now I will try to Z-scale transform data so it is with a given range.
#the correct order of operations as mentioned here (http://bonsai.hgc.jp/~mdehoon/software/cluster/manual/Data.html#Data) is log, then center, then normalize. Centering and normalization are both done using the "scale" function in R.
normlog3 <- t(scale(t(normlog2)))
#this works but gives some NAs.To remove all rows with NaN:
normlog4 = normlog3[complete.cases(normlog3), ]
normlog4 <- as.matrix(normlog4)

######################### Setting colors and establishing the heatmap function######################

#make the colors scale from blue to yellow
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 1000)

#from the website: http://compbio.ucsd.edu/making-heat-maps-r/
testHeatmap3<-function(logNC) {    
  #sampleColors = mapgeno2color(annotations)
  #heat<-heatmap.2(logNC)
  heatmap.2(logNC, margins=c(5,5), 
            col = my_palette,
            #ColSideColors = sampleColors,
            #labRow = e13fulldata2$GeneID,
            key.xlab="logNormCounts",
            key=TRUE, symkey=TRUE, density.info="none", trace="none")
  
  #return(ordered_gene_list) NOT WORKING
  }

################# Making the heatmap with the heatmap function and your transformed data #######################
dev.off()

testHeatmap3(normlog4)



