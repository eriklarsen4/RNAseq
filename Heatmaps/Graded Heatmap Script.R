

##### Heatmaps ####

library(readr)
library(tidyverse)
library(ggplot)
library(reshape2)
library(dplyr)
library(plyr)
library(readr)

##### GPCR Signaling Heatmap #####

## Import the Normalized Transcript Counts Excel CSV downloaded from Galaxy; contains the replicate normalized transcript reads that were used in a Fisher's Exact Test to determine the difference between WT and Tmem184b-mutant mice
RNASeqRepResults =  read_csv("M:/Erik/Data/Omics/RNAseq/Processed Galaxy Output/Counts Files to Upload/RNASeqRepResults.csv")

## Give the dataframe appropriate column names
colnames(RNASeqRepResults) = c("GeneID", "WT1", "WT2", "WT3", "WT4", "Mut1", "Mut2", "Mut3", "Mut4")

## Subset the dataframe by the results obtained in a ConsensusPathDB search, identifying GPCR Signaling as a significantly enriched pathway with the following genes enriched:
RNASeqRepResultsGPCRSig = RNASeqRepResults[RNASeqRepResults$GeneID %in% c("Dgkz","Gnaq","Dgkh","Lpar3","Dgkg","Cd55","Gnao1","Npy2r","Trpc3","Gna14","Gng2","Trpc6","Cysltr2","Adcy3","Prkca","Prkce","Camk4","Htr4","Prkcq","Gabbr2","Htr1a","Plcb3","Htr1f","Abr","Sst","F2rl2"),]


## Create dataframes to calculate Z-scores across replicates of each gene.
  ## These will hold the original data frame values until filled in later commands

  ## "Big" will store the means and sds of each gene
RNASeqRepResultsGPCRSigBIG = RNASeqRepResultsGPCRSig
  ## "3" will store only the means by genes for all replicates
RNASeqRepResultsGPCRSig3 = RNASeqRepResultsGPCRSig
  ## "4" will store only the sds by gene for all replicates
RNASeqRepResultsGPCRSig4 = RNASeqRepResultsGPCRSig
  ## "Z" will store only the Z-scores, which will be used directly to create the heatmaps
RNASeqRepResultsGPCRSigZ = RNASeqRepResultsGPCRSig

  ## Create a new column to fill with the means and standard deviations of each gene's transcript counts per million ("TPM")
RNASeqRepResultsGPCRSigBIG$mean = 0
RNASeqRepResultsGPCRSigBIG$sd = 0

  ## Loop through the data frame and fill the "mean" and "sd" columns with their appropriate values
for (i in 1:nrow(RNASeqRepResultsGPCRSigBIG)){
  RNASeqRepResultsGPCRSigBIG$mean[i] = (sum(RNASeqRepResultsGPCRSig[i,c(2:9)])/ncol(RNASeqRepResultsGPCRSig[,c(2:9)]))
  for (j in 1:nrow(RNASeqRepResultsGPCRSigBIG)){
    RNASeqRepResultsGPCRSigBIG$sd[j] = sd(RNASeqRepResultsGPCRSig[j,c(2:9)])
  }
}


  ## Create a dataframe storing the gene-specific mean normalized TPM in all columns/replicates for Z-score calculating
RNASeqRepResultsGPCRSig3[,c(2:9)] = RNASeqRepResultsGPCRSigBIG$mean
  ## Create a dataframe storing the gene-specific normalized TPM standard deviationsin all columns/replicates for Z-score calculating
RNASeqRepResultsGPCRSig4[,c(2:9)] = RNASeqRepResultsGPCRSigBIG$sd
  ## Create the Z-score dataframe
RNASeqRepResultsGPCRSigZ[,c(2:9)] = (RNASeqRepResultsGPCRSig[,c(2:9)] - RNASeqRepResultsGPCRSig3[,c(2:9)])/RNASeqRepResultsGPCRSig4[,c(2:9)]



## Re-arrange the data so that columns and rows can be run appropriately in a heatmap; give appropriate column names to this new dataframe
RNASeqRepResultsGPCRSig2 = melt(RNASeqRepResultsGPCRSigZ, id = "GeneID")
colnames(RNASeqRepResultsGPCRSig2) = c("GeneID", "Genotype", "TPM Z-score")

## Create and store the heatmap core as a variable
## Use geom_tile to map the transcript counts by Z-score
## Fill the scale gradient, with red = downregulation, green = upregulation; Title the legend
## Add labels

## Customize by removing the filler surrounding the graph and the tickmarks
## Center the Title
HEATER = ggplot(RNASeqRepResultsGPCRSig2, aes(RNASeqRepResultsGPCRSig2$Genotype, RNASeqRepResultsGPCRSig2$GeneID)) + geom_tile(aes(fill = RNASeqRepResultsGPCRSig2$`TPM Z-score`), color = "black") + scale_fill_gradient(low = "red", high = "green", name = "TPM Z-score", limits = c(-2,2)) + labs(x = "Genotype", y = "Gene", title = "Transcription Profile of GPCR Signaling Genes Significantly Enriched in Tmem184b Mutants") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0))

# Plot the heatmap
HEATER