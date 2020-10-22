
## Script for obtaining WGCNA and enrichr-based bioinformatics data
## Developed by Erik Larsen

  ## Condense lines of code by clicking all arrows next to code line #. Click again to open, step-wise to maintain organization

##### Environment Prep #####

  ## Install packages to use the "enrichR", "PANTHER.db" packages
    ## At least for the first installation, use the following code (THIS REQUIRES R version 4.0.0+)
    ## Most of the following packages require BiocManager for installation. Examples of installing the Bioconductor package for installation of subsequent Bioconductor packages are given below.

  ## Uncomment the following if/install statements if these packages have not previously been installed, and install them.
    #if (!requireNamespace("BiocManager", quietly = TRUE))
    #  install.packages("BiocManager")

    #if (!requireNamespace("BiocManager", quietly = TRUE))
    #  install.packages("BiocManager")
    #BiocManager::install("PANTHER.db")

    #BiocManager::install("EnhancedVolcano")

  ## Load the Panther 2016 Pathway Database package and use "mouse" as reference genome/specified species for downstream bioinformatics analysis. Wait for "PANTHER.db" to show in the Global Environment before continuing to load other packages.
library(PANTHER.db)
pthOrganisms(PANTHER.db) = "MOUSE"
# ## (Don't update if the option pops up
# n)

  ## Load the enrichR package
library(enrichR)
  ## Load previously installed packages from the library (if not previously installed, install first)
#library(tidyr)
#library(ggplot2)
library(tidyverse)
library(stringr)
library(plyr)
#library(ggrepel)
#library(calibrate)
library(reshape2)
#library(gridExtra)
#library(stats)
library(dplyr)
library(readr)
#library(readxl)
#library(tibble)
#library(AnnotationDbi)
library(stats4)
library(BiocGenerics)
#library(parallel)
library(dendextend)
library(ggdendro)
library(grid)
#library(EnhancedVolcano)
library(pheatmap)
library(matrixStats)
#library(Hmisc)
#library(splines)
#library(foreach)
#library(doParallel)
library(fastcluster)
library(dynamicTreeCut)
#library(survival)
#library(clValid)
library(WGCNA)
#library(biomaRt)
library(org.Mm.eg.db)
library(GO.db)
#library(bigmemory)
#library(anRichmentMethods)

#library(RColorBrewer)

  ## Connect live to the Enrichr server/master database (website) and store the "space" as a variable. This contains a vast array of bioinformatics databases/websites.
    ## Store it as a variable for better visualization and eventual subsetting
DBs = listEnrichrDbs()

  ## View the variable in the editor to find the relevant databases and their indeces for subsetting (click on the dataframe in the Global Environment and manually peruse)

    ## e.g. "GO_Biological_Process_2018" (index # 130), and "Panther_2016" (index # 102)
  ## Concatenate them in a list for subsetting, or slice directly
diff_DBs = DBs$libraryName[c(50,127,91,92)]
DBs = DBs$libraryName[c(130,131,132,102,14,110)]


##### Data Prep; import DGEA file and TPM file #####
  ## Import the dataset you want to analyze (Tmem184b-GT/GT DRG DGE)
DESeq2_Adults = read.csv("M:/Erik/Data/Omics/RNAseq/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
DESeq2_Adults3 = subset(DESeq2_Adults, (!is.na(DESeq2_Adults[,"AdjP"])))

  ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
DESeq2_Adults9 = DESeq2_Adults3 %>% filter(!grepl(DESeq2_Adults3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

#mt.+.?$|

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
DESeq2_Adults9$labs = DESeq2_Adults9$GeneID
    ## Replace DEGs with the string, "DEGs"
DESeq2_Adults9$labs[which(DESeq2_Adults9$AdjP <= 0.05)]= "DEGs"
    ## Replace the remaining genes with "Non-DEGs"
DESeq2_Adults9$labs[which(DESeq2_Adults9$AdjP >= 0.05)]= "Non-DEGs"

DEG_list = c(DESeq2_Adults9$GeneID[which(DESeq2_Adults9$labs == "DEGs")])
#keytypes(org.Mm.eg.db)
#columns(org.Mm.eg.db)
#keytypes(GO.db)
#columns(GO.db)

  ## Find Entrez IDs (Ensembl IDs) from the DEG list
Ensembl_DEG_IDs = as.integer(
  mapIds(org.Mm.eg.db,
         keys = as.character(
           c(DEG_list[])
           ),
         keytype = "SYMBOL",
         column = "ENTREZID")
)
  ## Extract GO Terms for all genes in the DEG list
DEG_GO_IDs = c(
  mapIds(org.Mm.eg.db, 
         keys = as.character(
           c(DEG_list[])
         ),
         keytype = "SYMBOL",
         column = "GOALL",
         multiVals = "list")
)

  ## Find the genes without Ensembl/Entrez IDs
which(is.na(Ensembl_DEG_IDs[]) == TRUE)
DEG_indeces_without_Ensembl_IDs = c(
  which(is.na(Ensembl_DEG_IDs[]) == TRUE)
)
  ## Repeat for GO IDs
DEG_indeces_without_GO_IDs = as.integer(c(
  which(is.na(DEG_GO_IDs[]) == TRUE)
))
  ## There should be more genes that **don't** have GO terms than genes that are un-annotated ("without GO ID" should be longer than "without Ensembl ID"; check)
DESeq2_Adults9$GeneID[c(DEG_indeces_without_Ensembl_IDs, DEG_indeces_without_GO_IDs)]
  ## Concatenate the genes into a list and rename the variable of the DESeq2 indeces with their gene names
DEGs_without_IDs = DESeq2_Adults9$GeneID[c(DEG_indeces_without_GO_IDs)]
  ## Use the indeces of the DEG_indeces_without_GO_IDs for looping if statement in the creation of the data frame below
#DEG_indeces_without_GO_IDs

#which(DEG_indeces_without_Ensembl_IDs == DEG_indeces_without_GO_IDs)
#match(DEG_indeces_without_Ensembl_IDs, DEG_indeces_without_GO_IDs)

  ## Create a matrix to house GO information for each gene
GO_INFO = matrix(nrow = length(c(DEG_list)), ncol = 3)
rownames(GO_INFO) = c(DEG_list)
colnames(GO_INFO) = c("# of GO Terms", "GO Term ID", "GO Term")

  ## Fill the first column with zeros for genes without GO Terms
GO_INFO[c(DEG_indeces_without_GO_IDs), 1] = 0
  ## Fill the second and third columns with "None"s for genes without GO Terms
GO_INFO[c(DEG_indeces_without_GO_IDs), c(2,3)] = "None"
  ## Fill the second column with the concatenated GO Term IDs associated with each gene, skipping over the genes with 0
for (i in 1:length(DEG_GO_IDs)){
    if (i == 58) {
      next
    }
    if (i == 121) {
      next
    }
    if (i == 140) {
      next
    }
    if (i == 158) {
      next
    }
    if (i == 159) {
      next
    }
    if (i == 172) {
      next
    }
    if (i == 196) {
      next
    }
    if (i == 234) {
      next
    }
    if (i == 316) {
      next
    }
    if (i == 334) {
      next
    }
    GO_INFO[i,2] = paste(names(mapIds(GO.db, DEG_GO_IDs[[i]][], "TERM", "GOID")), collapse = ";")
}
  ## Extract the GO IDs and remove the redundant ones;
    ## Fill the first column with the numbers of GO Terms associated with each gene, skipping over the genes with 0
for (i in 1:length(DEG_GO_IDs)){
  if (i == 58) {
    next
  }
  if (i == 121) {
    next 
  }
  if (i == 140) {
    next
  }
  if (i == 158) {
    next
  }
  if (i == 159) {
    next
  }
  if (i == 172) {
    next
  }
  if (i == 196) {
    next
  }
  if (i == 234) {
    next
  }
  if (i == 316) {
    next
  }
  if (i == 334) {
    next
  }
  x = str_split(as.vector(GO_INFO[i,2]), pattern = ";", simplify = TRUE)
  x = c(as.character(x[which(x %in% x == TRUE)]))
  GO_INFO[i,2] = paste(unique(x), collapse = ";")
  GO_INFO[i,1] = paste(as.numeric(length(unique(x))))
}
  ## Fill the third column with the concatenated GO Terms associated with each gene, skipping over the genes with 0
for (i in 1:length(DEG_GO_IDs)){
    if (i == 58) {
      next
    }
    if (i == 121) {
      next
    }
    if (i == 140) {
      next
    }
    if (i == 158) {
      next
    }
    if (i == 159) {
      next
    }
    if (i == 172) {
      next
    }
    if (i == 196) {
      next
    }
    if (i == 234) {
      next
    }
    if (i == 316) {
      next
    }
    if (i == 334) {
      next
    }
  x = str_split(as.vector(GO_INFO[i,2]), pattern = ";", simplify = TRUE)
  GO_INFO[i,3] = paste(as.character(mapIds(GO.db, keys = x, keytype = "GOID", "TERM")), collapse = ";")
}

  ## Extract all of the terms into one vector to find unique GO terms
for (i in 1:length(DEG_GO_IDs)){
  if (i == 58) {
    next
  }
  if (i == 121) {
    next
  }
  if (i == 140) {
    next
  }
  if (i == 158) {
    next
  }
  if (i == 159) {
    next
  }
  if (i == 172) {
    next
  }
  if (i == 196) {
    next
  }
  if (i == 234) {
    next
  }
  if (i == 316) {
    next
  }
  if (i == 334) {
    next
  }
  Unique_GOs = str_split(as.vector(GO_INFO[,3]), pattern = ";", simplify = TRUE)
  Unique_GOs = c(as.character(Unique_GOs[which(Unique_GOs %in% Unique_GOs == TRUE)]))
  Unique_GOs = paste(unique(Unique_GOs), collapse = ";")
  Unique_GOs = c(str_split(as.vector(Unique_GOs), pattern = ";", simplify = TRUE))
}

  ## Make the matrix a data frame, and add "GeneIDs" as the first column
GO_INFO_df = as.data.frame(GO_INFO)
GO_INFO_df$GeneID = rownames(GO_INFO_df)
rownames(GO_INFO_df) = NULL
col_idx = grep("GeneID", names(GO_INFO_df))
GO_INFO_df = GO_INFO_df[ , c(col_idx, (1:ncol(GO_INFO_df))[-c(col_idx)]) ]


  ## Search for terms of individual genes
#searchidx = as.numeric(which(GO_INFO_df$GeneID == "Lgals9"))
#y = str_split(as.vector(GO_INFO[searchidx,4]), pattern = ";", simplify = TRUE)

  ## Search for terms associated with multiple genes
#str_split(as.vector(GO_INFO[2,4]), pattern = ";", simplify = TRUE)[
#  which(
#  (str_split(
#    as.vector(GO_INFO[2,4]), pattern = ";", simplify = TRUE) %in%
#     str_split(
#       as.vector(GO_INFO[4,4]), pattern = ";", simplify = TRUE)
#   ) == TRUE)
#]

#str_split(as.vector(GO_INFO[104,3]), pattern = ";", simplify = TRUE)[]
#DEG_list[104]

str_split(as.vector(GO_INFO[ which(GO_INFO_df$GeneID == "Rab3d") , 3]), pattern = ";", simplify = TRUE)
GO_INFO_df$GeneID[
  str_split(as.vector(GO_INFO[ which(GO_INFO_df$GeneID == "Rab3d") , 3]), pattern = ";", simplify = TRUE)[148] %in% str_split(as.vector(GO_INFO[,3]), pattern = ";", simplify = TRUE)  %in% str_split(as.vector(GO_INFO[,3]), pattern = ";", simplify = TRUE)]

  ## Find the GO Terms of a desired gene from the DEG list
str_split(
  as.vector(GO_INFO[ which(GO_INFO_df$GeneID == "Rab3d") ,   3]), pattern = ";", simplify = TRUE)[
    which(str_split(
      as.vector(GO_INFO[ which(GO_INFO_df$GeneID == "Rab3d") ,    3]), pattern = ";", simplify = TRUE) %in% Unique_GOs)
    ]

  ## Find GO Terms in a gene of interest also in any other genes
GO_INFO_df$GeneID[which(c(str_split(
  as.vector(GO_INFO[ which(GO_INFO_df$GeneID == "Trpa1") ,   3]), pattern = ";", simplify = TRUE)[
    which(str_split(
      as.vector(GO_INFO[ which(GO_INFO_df$GeneID == "Trpa1") ,    3]), pattern = ";", simplify = TRUE) %in% Unique_GOs)
  ]
  ) %in% str_split(as.vector(GO_INFO[,3]), pattern = ";", simplify = TRUE)
  )
]




    ## Import the TPM file. These are expression estimates for each gene, for each sample/replicate, where each gene's value is normalized to its sample's effect size
Adult_Normalized_Counts = read_csv("M:/Erik/Data/Omics/RNAseq/Processed Galaxy Output/Counts Files to Upload/RNASeqRepResults.csv", col_names = TRUE)
  ## Rename columns; should know this ahead of time
colnames(Adult_Normalized_Counts) = c("GeneID", "WT1", "WT2", "WT3", "WT4", "Mut1", "Mut2", "Mut3", "Mut4")


  ## Subset the genes most affected by Tmem mutation (4000+ is computationally dense)
DESeq2_Short = DESeq2_Adults9[c(1:5000),]
  ## Subset the TPM file by the 4000 most-affected genes
Short_Profile = Adult_Normalized_Counts[Adult_Normalized_Counts$GeneID %in% c(DESeq2_Short$GeneID),]


##### Generate heatmaps from hierarchically clustered data (Z-scored or not; based on Euclidean distance) #####
  ##### For Z #####
  ## Create a function that determines the Zscores of the expression profile (row Z-scores)
Find_Row_Z = function(Expression_Profile){
  ## Create dataframes to calculate Z-scores across replicates of each gene.
  ## These will hold the original data frame values until filled in later commands
  ## "Big" will store the means and sds of each gene
  NetBIG = Expression_Profile
  ## "3" will store only the means by genes for all replicates
  Net3 = Expression_Profile
  ## "4" will store only the sds by gene for all replicates
  Net4 = Expression_Profile
  ## "Z" will store only the Z-scores, which will be used directly to create the heatmaps
  NetZ = Expression_Profile
  
  ## Create a new column to fill with the means and standard deviations of each gene's transcript counts per million ("TPM")
  NetBIG$mean = 0
  NetBIG$sd = 0
  
  ## Loop through the dataframe and fill the "mean" and "sd" columns with their appropriate values
  for (i in 1:nrow(NetBIG)){
    NetBIG$mean[i] = (sum(Expression_Profile[i,c(2:ncol(Expression_Profile))])/ncol(Expression_Profile[,c(2:ncol(Expression_Profile))]))
  }
  for (j in 1:nrow(NetBIG)){
    NetBIG$sd[j] = sd(Expression_Profile[j,c(2:ncol(Expression_Profile))], na.rm = TRUE)
  }
  ## Create a dataframe storing the gene-specific mean normalized TPM in all columns/replicates for Z-score calculating
  Net3[,c(2:ncol(Expression_Profile))] = NetBIG$mean
  ## Create a dataframe storing the gene-specific normalized TPM standard deviationsin all columns/replicates for Z-score calculating
  Net4[,c(2:ncol(Expression_Profile))] = NetBIG$sd
  ## Create the Z-score dataframe
  NetZ[,c(2:ncol(Expression_Profile))] = (Expression_Profile[,c(2:ncol(Expression_Profile))] - Net3[,c(2:ncol(Expression_Profile))])/Net4[,c(2:ncol(Expression_Profile))]
  ## Remove the genes that were detected to have 0 TPMs across all samples
  ## Create a function that evaluates a vector/dataframe's (x's) numerical values. Returns equal length vector with T/F bools.
  ## Subset the dataframe that filters those genes.
  row_has_na = apply(NetZ, 1, function(x){any(is.na(x))})
  NetZ = NetZ[!row_has_na,]
  Z <<- NetZ
}
  ## Determine Row Z-score (gene-wise Z)
Find_Row_Z(Expression_Profile = Short_Profile)

  ## Perform hierarchical clustering on the Z-score transcriptional profiles; Use Euclidean distance (Correlation distance is equivalent to Euclidean Z)
Euclid_dist_order_Z = hclust(dist(Z[,2:9], method = "euclidean"))$order
  ## Find the gene names of Z-scored, hierarchically-clustered-and-ordered-by-Euclidean-distance TPM data
Euclid_dist_ord_Z_Genes = c(Z$GeneID[Euclid_dist_order_Z])
  ## Re-arrange the Z-clustered profiles into a format compatible with heatmap construction
Corr_Z = Z %>%
  mutate(GeneID = factor(GeneID, levels = Euclid_dist_ord_Z_Genes)) %>%
  arrange(GeneID)
  ## Find where Tmem is
which(Euclid_dist_ord_Z_Genes == "Tmem184b") # 2904

  ## Create a heatmap around Tmem
pheatmap(mat = Corr_Z[2884:2924,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, labels_row = Euclid_dist_ord_Z_Genes[2884:2924])
  ## Create a heatmap of all 4000 genes
pheatmap(mat = Corr_Z[,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, show_rownames = F)

  ##### For non-centered data #####
  ## Perform hierarchical clustering on the un-centered transcriptional profiles; use Euclicdean distance
Euclid_dist_order = hclust(dist(Short_Profile[,2:9], method = "euclidean"))$order
  ## Find gene names
Euclid_dist_ord_Genes = c(Short_Profile$GeneID[Euclid_dist_order])
  ## Re-arrange the clustered profiles into a format compatible with heatmap construction
Corr = Short_Profile %>%
  mutate(GeneID = factor(GeneID, levels = Euclid_dist_ord_Genes)) %>%
  arrange(GeneID)
  ## Find where Tmem is
which(Euclid_dist_ord_Genes == "Tmem184b") # 3213

  ## Create a heatmap around Tmem
pheatmap(mat = Corr[3193:3233,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, labels_row = Euclid_dist_ord_Genes[3193:3233])



##### Data prep for manual network creation; find the beta value/power to raise the adjacency #####
  ## Remove genes that generate NANs ##
  ## Create a function that takes the Expression Data, which downstream analysis to perform (WT, Mut, or All), and removes genes with too many 0s that will create NANs (dividing by 0)
Remove_NAN_Generators = function(Expression_Profile, Samples_to_Compare){
  
  ## Classify which analysis to determine dataframe subset
  
  #MutAdultItchExp = RNASeqRepResultsAdultAll[,grep(names(RNASeqRepResultsAdultAll), pattern = "GeneID|Mut")]
  #WTAdultItchExp = RNASeqRepResultsAdultAll %>% select(starts_with(GeneID) & ends_with(WT4))
  
  if (Samples_to_Compare == "WT" ){
    Expression_Profile = Expression_Profile[,grep(names(Expression_Profile), pattern = "GeneID|WT")]
  } else if (Samples_to_Compare == "Mut"){
    Expression_Profile = Expression_Profile[,grep(names(Expression_Profile), pattern = "GeneID|Mut")]
  } else if (Samples_to_Compare == "All"){
    Expression_Profile = Expression_Profile
  } else {
    Expression_Profile = Expression_Profile
  }
  ## Run a for-loop to evaluate whether any genes in the dataframe contain all 0s. Store it in the above initialized vector
  ## For an entire profile
  if (ncol(Expression_Profile) > 5){
    for (i in 1:nrow(Expression_Profile)){
      expcheck[i] = if_else(Expression_Profile[i,2] == 0 & 
                              Expression_Profile[i,3] == 0 & 
                              Expression_Profile[i,4] == 0 & 
                              Expression_Profile[i,5] == 0 &
                              Expression_Profile[i,6] == 0 &
                              Expression_Profile[i,7] == 0 &
                              Expression_Profile[i,8] == 0 &
                              Expression_Profile[i,9] == 0, "Remove", "Keep")
    }
  } 
  ## For Genotype-specific profile
  else {
    for (i in 1:nrow(Expression_Profile)){
      expcheck[i] = if_else(Expression_Profile[i,2] == 0 & 
                              Expression_Profile[i,3] == 0 & 
                              Expression_Profile[i,4] == 0 & 
                              Expression_Profile[i,5] == 0, "Remove", "Keep")
    }
  }
  ## Find the indeces of those genes to remove (they create NaNs and prevent correlation calculation and hierarchical clustering)
  which(expcheck == "Remove")
  Expression_Profile$GeneID[7816]
  
  #Removables = Expression_Profile[c(1873,2665,2730,2841,3103,3141,3989,4478,4515,4745,5027,5288,5501,5745,5761,5974,5988,6428,6776,7630,7816,7925,8143,9102,9147,9670,9825,10205,10530,10935),]  
  
  ## Filter out "Rest", a developmental TF with 0s across all WT replicates
  which(Expression_Profile$GeneID == "Rest")
  
  if (Samples_to_Compare == "WT") {
    Expression_Profile_cleaned = Expression_Profile[-7816,]
  }else if (Samples_to_Compare == "Mut") {
    Expression_Profile_cleaned = Expression_Profile[-c(4956,8485),]
  }else if (Samples_to_Compare == "All") {
    Expression_Profile_cleaned = Expression_Profile
  }else{
    Expression_Profile_cleaned = Expression_Profile
  }
  
  Expression_Profile_cleaned <<- Expression_Profile_cleaned
}

  ## Create a vector as long as the expression dataframe for determining additional genes to remove upstream of correlation calculations
expcheck = as.array(vector(length = nrow(Short_Profile)))
  ## Remove NAN-generating genes
Remove_NAN_Generators(Expression_Profile = Short_Profile, Samples_to_Compare = "All")

  ## Take a range of potential powers to determine the appropriate power by which the correlation/adjacency matrix will be raised. This can be thought of as increasing the "connectivity" of more related genes
powers = c(seq(4,10,by=1), seq(12,20, by=2))
  ## Use the WGCNA function, "pickSoftThreshold", to compute model fits of various powers, and to compute the "connectivity" metrics of the data based on these power values
powerTable = list(data = pickSoftThreshold(data = t(Expression_Profile_cleaned[,2:9]), powerVector = powers, verbose = 2))
collectGarbage()
colors = c("black")
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")

  ## Re-arrange the data output of the soft thresholding function into a dataframe for graphics and analysis
Topology_df = data.frame(powerTable$data$fitIndices[1],
                         powerTable$data$fitIndices[2], 
                         powerTable$data$fitIndices[3], 
                         powerTable$data$fitIndices[4], 
                         powerTable$data$fitIndices[5], 
                         powerTable$data$fitIndices[6], 
                         powerTable$data$fitIndices[7], stringsAsFactors = F)
powerTable$data$powerEstimate
powerlabels = c(Topology_df$Power)

sizeGrWindow(8,6)
par(mfcol = c(2,2))
par(mar = c(4.2,4.2,2.2,0.5))
cex1 = 0.7

  ## Graph the Scale-Free model fit as a function of soft threshold (power)
plot(Topology_df$Power, -Topology_df$slope*Topology_df$SFT.R.sq,
     xlab = "Soft Threshold (power)",
     ylab = colNames[1], type = "n", 
     xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
     ylim = c(0,5), main = colNames[1])
addGrid()
text(Topology_df$Power, -Topology_df$slope*Topology_df$SFT.R.sq, labels = powers)

  ## Graph the mean connectivity as a function of soft threshold (power)
plot(Topology_df$mean.k., -Topology_df$slope*Topology_df$mean.k.,
     xlab = "Soft Threshold (power)",
     ylab = colNames[2], type = "n",
     xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
     ylim = c(min(-Topology_df$slope*Topology_df$mean.k.), max(-Topology_df$slope*Topology_df$mean.k.)), main = colNames[2])
addGrid()
text(Topology_df$Power, -(Topology_df$slope)*Topology_df$mean.k., labels = powers)

  ## Graph the median connectivity as a function of soft threshold (power)
plot(Topology_df$median.k., -Topology_df$slope*Topology_df$median.k.,
     xlab = "Soft Threshold (power)", ylab = colNames[3], type = "n",
     xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
     ylim = c(min(-Topology_df$slope*Topology_df$median.k.), max(-Topology_df$slope*Topology_df$median.k.)), main = colNames[3])
addGrid()
text(Topology_df$Power, -(Topology_df$slope)*Topology_df$median.k., labels = powers)

  ## Graph the maximum connectivity as a function of soft threshold (power)
plot(Topology_df$max.k., -Topology_df$slope*Topology_df$max.k.,
     xlab = "Soft Threshold (power)", ylab = colNames[4], type = "n",
     xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
     ylim = c(min(-Topology_df$slope*Topology_df$max.k.), max(-Topology_df$slope*Topology_df$max.k.)), main = colNames[4])
addGrid()
text(Topology_df$Power, -(Topology_df$slope)*Topology_df$max.k., labels = powers)
  #### Balance the maximum connectivity loss with the scale free topology gain. --> Use power of 4-7

  ##### Create the adjacency/TOM #####
  ## Create a weighted adjacency (no absolute value of the correlations raised to a power)
ADJ = cor(t(Expression_Profile_cleaned[,2:9]), method = "pearson")^6
rownames(ADJ) = Expression_Profile_cleaned$GeneID
colnames(ADJ) = Expression_Profile_cleaned$GeneID

  ## Create the topological overlap matrix using the WGCNA function "TOMsimilarity": converts an adjacency into a TOM
TOM = TOMsimilarity(adjMat = ADJ, TOMType = "signed")
rownames(TOM) = rownames(ADJ)
colnames(TOM) = colnames(ADJ)
  ## Create the topological overlap dissimilarity matrix
disTOM = 1-TOMsimilarity(adjMat = ADJ, TOMType = "signed")
rownames(disTOM) = rownames(ADJ)
colnames(disTOM) = rownames(ADJ)

  ##### Visualize #####
  ## Create a hierarchically clustered dendrogram based on the DISsimilarity (of the TOM); distinguishes clusters better than similarity
Tree = hclust(as.dist(disTOM), method = "average")
  ## Determine clusters to plot based on a "dynamic dendrogram cut"
unmergedLabs = cutreeDynamic(dendro = Tree, distM = disTOM,
                             deepSplit = 1,
                             cutHeight = 0.995,
                             pamRespectsDendro = F)

  ## Convert the clusters from index labels to colors
unmergedCols = labels2colors(unmergedLabs)
  ## Merge clusters that are similar enough based on the subjective/arbitrary "cutHeight" parameter provided within the WGCNA function, "mergeCloseModules"
merge = mergeCloseModules(t(Expression_Profile_cleaned[,2:9]), unmergedLabs, cutHeight = 0.10, verbose = 2)
  ## Convert color indeces into a string
moduleLabs = merge$colors
  ## Convert the string of labels into colors
moduleCols = labels2colors(moduleLabs)
  ## Find the new, merged eigenvectors (modules)
MEs = merge$newMEs

#t(Expression_Profile_cleaned[,2:9])

  ## Find Tmem index
which(Expression_Profile_cleaned$GeneID == "Tmem184b")
  ## Plot the dendrogram and associated clusters/modules
plotDendroAndColors(Tree, cbind(unmergedCols, moduleCols), 
                    c("Unmerged Modules", "Merged Modules"), 
                    dendroLabels = F, 
                    hang = 0.001,
                    addGuide = T,
                    guideHang = 0.001,
                    main = "Tmem GT DRG Dendrogram")

which(Expression_Profile_cleaned$GeneID == "Il31ra")
  ## Determine Tmem's modules
moduleCols[4395]
unmergedCols[4395]
#moduleCols[1537]

#Expression_Profile_cleaned$GeneID[c(which(moduleCols == "sienna3"))]
  ## Find the indeces from the clustered data of genes within Tmem's clusters
TmemUnmergedModule = c(which(unmergedCols == "yellow"))
TmemModule = c(which(moduleCols == "brown"))
  ## Find the genes of those indeces
TmemUnmergedModule = c(Expression_Profile_cleaned$GeneID[TmemUnmergedModule])
TmemModule = c(Expression_Profile_cleaned$GeneID[TmemModule])


  ## Create a dendrogram compatible with TOM-TOM plot creation
colorDynamicTOM = labels2colors(cutreeDynamic(Tree, method ="tree"))
  ## Make the diagonal of the TOM NA to re-scale the plot
diag(disTOM) = 0
  ## Create the TOM-TOM plot
TOMplot((disTOM^3), Tree, as.character(moduleCols), main = "Dissim. TOM Plot Top5000")

#colorDynamicADJ = labels2colors(cutreeDynamic(manTree, method = "tree"))
#diag(ADJ) = 1
#TOMplot((ADJ), manTree, as.character(colorDynamicADJ), main = "TOM Plot based on ADJ Sim.")

##### Functional Analysis #####
  ## Create a function that will return Enrichr-based functional analysis on a given gene set.
  ## First, create dataframes to store returned analyses as variables
GO_Processes = data.frame()
GO_Cell_Comps = data.frame()
GO_Mol_Funcs = data.frame()
PATHWAYS = data.frame()
PPIs = data.frame()
DISEASES = data.frame()
ENCODE_TFs = data.frame()
Enrichr_TFs = data.frame()
GPCR_Down_Perturbs = data.frame()
GPCR_Up_Perturbs = data.frame()
GO_PROCESS_GENES = data.frame()
GO_MOL_FUNC_GENES = data.frame()
GO_CELL_COMP_GENES = data.frame()
PATHWAY_GENES = data.frame()
PPI_GENES = data.frame()
DISEASE_GENES = data.frame()
ENCODE_TF_GENES = data.frame()
ENRICHR_TF_GENES = data.frame()
GPCR_DOWN_PERTURB_GENES = data.frame()
GPCR_UP_PERTURB_GENES = data.frame()

Enrichr_Functional_Output = data.frame()

  ## Create the function
Functional_Analysis = function(Genes){
  ## first: Gene ontology biological processes; this identifies cell/biological processes with which the genes in the dataset are associated
  GO_Processes = as.data.frame(
    enrichr(
      c(Genes), DBs[1])
  )
  ## second: Gene ontology cellular components; this identifies the cellular organelles in which the genes in the dataset reside
  GO_Cell_Comps = as.data.frame(
    enrichr(
      c(Genes), DBs[2])
  )
  ## third: Gene ontology molecular functions; this identifies molecular functions of the protein products related to the genes in the dataset
  GO_Mol_Funcs = as.data.frame(
    enrichr(
      c(Genes), DBs[3])
  )
  ## fourth: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
  PATHWAYS = as.data.frame(
    enrichr(
      c(Genes), DBs[4])
  )
  ## fifth: PPIs; this identifies major "hub" genes' protein products- what they interact with, etc.
  PPIs = as.data.frame(
    enrichr(
      c(Genes), DBs[5])
  )
  ## last: diseases; this identifies phenotypes or diseases associated with the genes provided
  DISEASES = as.data.frame(
    enrichr(
      c(Genes), DBs[6])
  )
  
  ## Additional, probably less robust bioinformatics analyses, involving ChIP seq/ChIP, and GPCR/kinase regulation based on provided gene lists
  ENCODE_TFs = as.data.frame(
    enrichr(
      c(Genes), diff_DBs[1])
  )
  
  Enrichr_TFs = as.data.frame(
    enrichr(
      c(Genes), diff_DBs[2])
  )
  
  GPCR_Down_Perturbs = as.data.frame(
    enrichr(
      c(Genes), diff_DBs[3])
  )
  
  GPCR_Up_Perturbs = as.data.frame(
    enrichr(
      c(Genes), diff_DBs[4])
  )
  
  ## Remove unnecessary columns
  GO_Processes = GO_Processes[,-c(5,6,7)]
  GO_Cell_Comps = GO_Cell_Comps[,-c(5,6,7)]
  GO_Mol_Funcs = GO_Mol_Funcs[,-c(5,6,7)]
  PATHWAYS = PATHWAYS[,-c(5,6,7)]
  PPIs = PPIs[,-c(5,6,7)]
  DISEASES = DISEASES[,-c(5,6,7)]
  ENCODE_TFs = ENCODE_TFs[,-c(5,6,7)]
  Enrichr_TFs = Enrichr_TFs[,-c(5,6,7)]
  GPCR_Down_Perturbs = GPCR_Down_Perturbs[,-c(5,6,7)]
  GPCR_Up_Perturbs = GPCR_Up_Perturbs[,-c(5,6,7)]
  
  
  ## Clean up the dataframes
  GO_Processes = GO_Processes %>%
    separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
  GO_Processes$Overlap = as.numeric(GO_Processes$Overlap)
  GO_Processes$`Process Size` = as.numeric(GO_Processes$`Process Size`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  GO_Processes = add_column(GO_Processes, EnrichrZscore = GO_Processes$GO_Biological_Process_2018.Combined.Score/log(GO_Processes$GO_Biological_Process_2018.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  GO_Processes = add_column(GO_Processes, Weighted_Overlap_Ratio = GO_Processes$Overlap*(GO_Processes$Overlap/GO_Processes$`Process Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  GO_Processes = add_column(GO_Processes, Modified_Combined_Score = (abs(GO_Processes$EnrichrZscore))*(-log(GO_Processes$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)
  ## Remove GO IDs
  GO_Processes$GO_Biological_Process_2018.Term = GO_Processes$GO_Biological_Process_2018.Term %>% gsub(x = GO_Processes$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")
  
  
  GO_Cell_Comps = GO_Cell_Comps %>%
    separate(GO_Cellular_Component_2018.Overlap, c("Overlap", "Number of Components"), "/")
  ## Convert the values into integers for computation
  GO_Cell_Comps$Overlap = as.numeric(GO_Cell_Comps$Overlap)
  GO_Cell_Comps$`Number of Components` = as.numeric(GO_Cell_Comps$`Number of Components`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  GO_Cell_Comps = add_column(GO_Cell_Comps, EnrichrZscore = GO_Cell_Comps$GO_Cellular_Component_2018.Combined.Score/log(GO_Cell_Comps$GO_Cellular_Component_2018.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  GO_Cell_Comps = add_column(GO_Cell_Comps, Weighted_Overlap_Ratio = GO_Cell_Comps$Overlap*(GO_Cell_Comps$Overlap/GO_Cell_Comps$`Number of Components`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  GO_Cell_Comps = add_column(GO_Cell_Comps, Modified_Combined_Score = (abs(GO_Cell_Comps$EnrichrZscore))*(-log(GO_Cell_Comps$GO_Cellular_Component_2018.Adjusted.P.value)), .before = 6)
  ## Remove GO IDs
  GO_Cell_Comps$GO_Cellular_Component_2018.Term = GO_Cell_Comps$GO_Cellular_Component_2018.Term %>% gsub(x = GO_Cell_Comps$GO_Cellular_Component_2018.Term, pattern = " \\(.*\\)$", replacement = "")
  
  
  GO_Mol_Funcs = GO_Mol_Funcs %>%
    separate(GO_Molecular_Function_2018.Overlap, c("Overlap", "Function Size"), "/")
  ## Convert the values into integers for computation
  GO_Mol_Funcs$Overlap = as.numeric(GO_Mol_Funcs$Overlap)
  GO_Mol_Funcs$`Function Size` = as.numeric(GO_Mol_Funcs$`Function Size`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  GO_Mol_Funcs = add_column(GO_Mol_Funcs, EnrichrZscore = GO_Mol_Funcs$GO_Molecular_Function_2018.Combined.Score/log(GO_Mol_Funcs$GO_Molecular_Function_2018.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  GO_Mol_Funcs = add_column(GO_Mol_Funcs, Weighted_Overlap_Ratio = GO_Mol_Funcs$Overlap*(GO_Mol_Funcs$Overlap/GO_Mol_Funcs$`Function Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  GO_Mol_Funcs = add_column(GO_Mol_Funcs, Modified_Combined_Score = (abs(GO_Mol_Funcs$EnrichrZscore))*(-log(GO_Mol_Funcs$GO_Molecular_Function_2018.Adjusted.P.value)), .before = 6)
  ## Remove GO IDs
  GO_Mol_Funcs$GO_Molecular_Function_2018.Term = GO_Mol_Funcs$GO_Molecular_Functiont_2018.Term %>% gsub(x = GO_Mol_Funcs$GO_Molecular_Function_2018.Term, pattern = " \\(.*\\)$", replacement = "")
  
  
  ## Repeat the above code for Panther Pathways
  PATHWAYS = PATHWAYS %>%
    separate(Panther_2016.Overlap, c("Overlap", "Pathway Size"), "/")
  ## Convert the values into integers for computation
  PATHWAYS$Overlap = as.numeric(PATHWAYS$Overlap)
  PATHWAYS$`Pathway Size` = as.numeric(PATHWAYS$`Pathway Size`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  PATHWAYS = add_column(PATHWAYS, EnrichrZscore = PATHWAYS$Panther_2016.Combined.Score/log(PATHWAYS$Panther_2016.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  PATHWAYS = add_column(PATHWAYS, Weighted_Overlap_Ratio = PATHWAYS$Overlap*(PATHWAYS$Overlap/PATHWAYS$`Pathway Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  PATHWAYS = add_column(PATHWAYS, Modified_Combined_Score = (abs(PATHWAYS$EnrichrZscore))*(-log(PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)
  ## Remove the Panther Pathways species name
  PATHWAYS$Panther_2016.Term = PATHWAYS$Panther_2016.Term %>% gsub(x = PATHWAYS$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")
  
  
  ## Repeat the above code for PPIs
  PPIs = PPIs %>%
    separate(PPI_Hub_Proteins.Overlap, c("Overlap", "Hub Size"), "/")
  ## Convert the values into integers for computation
  PPIs$Overlap = as.numeric(PPIs$Overlap)
  PPIs$`Hub Size` = as.numeric(PPIs$`Hub Size`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  PPIs = add_column(PPIs, EnrichrZscore = PPIs$PPI_Hub_Proteins.Combined.Score/log(PPIs$PPI_Hub_Proteins.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway/ ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  PPIs = add_column(PPIs, Weighted_Overlap_Ratio = PPIs$Overlap*(PPIs$Overlap/PPIs$`Hub Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  PPIs = add_column(PPIs, Modified_Combined_Score = (abs(PPIs$EnrichrZscore))*(-log(PPIs$PPI_Hub_Proteins.Adjusted.P.value)), .before = 6)
  # ## Change the terms to lower case for future analysis
  # for (i in 1:nrow(PPIs)){
  #   PPIs[i,1] = paste(substr(PPIs[i,1], 1, 1),
  #                    tolower(substr(PPIs[i,1], 2, 7)), sep = "")
  #}
  
  ## Repeat for diseases
  DISEASES = DISEASES %>%
    separate(Jensen_DISEASES.Overlap, c("Overlap", "Total Disease Genes"), "/")
  ## Convert the values into integers for computation
  DISEASES$Overlap = as.numeric(DISEASES$Overlap)
  DISEASES$`Total Disease Genes` = as.numeric(DISEASES$`Total Disease Genes`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  DISEASES = add_column(DISEASES, EnrichrZscore = DISEASES$Jensen_DISEASES.Combined.Score/log(DISEASES$Jensen_DISEASES.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  DISEASES = add_column(DISEASES, Weighted_Overlap_Ratio = DISEASES$Overlap*(DISEASES$Overlap/DISEASES$`Total Disease Genes`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  DISEASES = add_column(DISEASES, Modified_Combined_Score = (abs(DISEASES$EnrichrZscore))*(-log(DISEASES$Jensen_DISEASES.Adjusted.P.value)), .before = 6)
  
  
  
  ## Repeat for TFs identified by ENCODE (Encyclopedia of DNA Elements)
  ENCODE_TFs = ENCODE_TFs %>%
    separate(ENCODE_TF_ChIP.seq_2015.Overlap, c("Overlap", "Size"), "/")
  ## Convert the values into integers for computation
  ENCODE_TFs$Overlap = as.numeric(ENCODE_TFs$Overlap)
  ENCODE_TFs$`Size` = as.numeric(ENCODE_TFs$`Size`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  ENCODE_TFs = add_column(ENCODE_TFs, EnrichrZscore = ENCODE_TFs$ENCODE_TF_ChIP.seq_2015.Combined.Score/log(ENCODE_TFs$ENCODE_TF_ChIP.seq_2015.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway/ ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  ENCODE_TFs = add_column(ENCODE_TFs, Weighted_Overlap_Ratio = ENCODE_TFs$Overlap*(ENCODE_TFs$Overlap/ENCODE_TFs$`Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  ENCODE_TFs = add_column(ENCODE_TFs, Modified_Combined_Score = (abs(ENCODE_TFs$EnrichrZscore))*(-log(ENCODE_TFs$ENCODE_TF_ChIP.seq_2015.Adjusted.P.value)), .before = 6)
  
  
  ## Repeat for TFs that emerged in other Enrichr submissions containing similar genes
  Enrichr_TFs = Enrichr_TFs %>%
    separate(Enrichr_Submissions_TF.Gene_Coocurrence.Overlap, c("Overlap", "Genes from Submissions Regulated by a TF"), "/")
  ## Convert the values into integers for computation
  Enrichr_TFs$Overlap = as.numeric(Enrichr_TFs$Overlap)
  Enrichr_TFs$`Genes from Submissions Regulated by a TF` = as.numeric(Enrichr_TFs$`Genes from Submissions Regulated by a TF`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  Enrichr_TFs = add_column(Enrichr_TFs, EnrichrZscore = Enrichr_TFs$Enrichr_Submissions_TF.Gene_Coocurrence.Combined.Score/log(Enrichr_TFs$Enrichr_Submissions_TF.Gene_Coocurrence.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway/ ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  Enrichr_TFs = add_column(Enrichr_TFs, Weighted_Overlap_Ratio = Enrichr_TFs$Overlap*(Enrichr_TFs$Overlap/Enrichr_TFs$`Genes from Submissions Regulated by a TF`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  Enrichr_TFs = add_column(Enrichr_TFs, Modified_Combined_Score = (abs(Enrichr_TFs$EnrichrZscore))*(-log(Enrichr_TFs$Enrichr_Submissions_TF.Gene_Coocurrence.Adjusted.P.value)), .before = 6)
  
  
  ## Repeat the above code for TFs that result in other submissions containing similar genes
  GPCR_Down_Perturbs = GPCR_Down_Perturbs %>%
    separate(L1000_Kinase_and_GPCR_Perturbations_down.Overlap, c("Overlap", "Genes Specifically Involved in given GPCR/Kinase Activity"), "/")
  ## Convert the values into integers for computation
  GPCR_Down_Perturbs$Overlap = as.numeric(GPCR_Down_Perturbs$Overlap)
  GPCR_Down_Perturbs$`Genes Specifically Involved in given GPCR/Kinase Activity` = as.numeric(GPCR_Down_Perturbs$`Genes Specifically Involved in given GPCR/Kinase Activity`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  GPCR_Down_Perturbs = add_column(GPCR_Down_Perturbs, EnrichrZscore = GPCR_Down_Perturbs$L1000_Kinase_and_GPCR_Perturbations_down.Combined.Score/log(GPCR_Down_Perturbs$L1000_Kinase_and_GPCR_Perturbations_down.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway/ ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  GPCR_Down_Perturbs = add_column(GPCR_Down_Perturbs, Weighted_Overlap_Ratio = GPCR_Down_Perturbs$Overlap*(GPCR_Down_Perturbs$Overlap/GPCR_Down_Perturbs$`Genes Specifically Involved in given GPCR/Kinase Activity`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  GPCR_Down_Perturbs = add_column(GPCR_Down_Perturbs, Modified_Combined_Score = (abs(GPCR_Down_Perturbs$EnrichrZscore))*(-log(GPCR_Down_Perturbs$L1000_Kinase_and_GPCR_Perturbations_down.Adjusted.P.value)), .before = 6)
  
  
  ## Repeat the above code for TFs that result in other submissions containing similar genes
  GPCR_Up_Perturbs = GPCR_Up_Perturbs %>%
    separate(L1000_Kinase_and_GPCR_Perturbations_up.Overlap, c("Overlap", "Genes Specifically Involved in given GPCR/Kinase Activity"), "/")
  ## Convert the values into integers for computation
  GPCR_Up_Perturbs$Overlap = as.numeric(GPCR_Up_Perturbs$Overlap)
  GPCR_Up_Perturbs$`Genes Specifically Involved in given GPCR/Kinase Activity` = as.numeric(GPCR_Up_Perturbs$`Genes Specifically Involved in given GPCR/Kinase Activity`)
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  GPCR_Up_Perturbs = add_column(GPCR_Up_Perturbs, EnrichrZscore = GPCR_Up_Perturbs$L1000_Kinase_and_GPCR_Perturbations_up.Combined.Score/log(GPCR_Up_Perturbs$L1000_Kinase_and_GPCR_Perturbations_up.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
  ## The fraction of differentially expressed genes belonging to a term/pathway/ ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  GPCR_Up_Perturbs = add_column(GPCR_Up_Perturbs, Weighted_Overlap_Ratio = GPCR_Up_Perturbs$Overlap*(GPCR_Up_Perturbs$Overlap/GPCR_Up_Perturbs$`Genes Specifically Involved in given GPCR/Kinase Activity`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  GPCR_Up_Perturbs = add_column(GPCR_Up_Perturbs, Modified_Combined_Score = (abs(GPCR_Up_Perturbs$EnrichrZscore))*(-log(GPCR_Up_Perturbs$L1000_Kinase_and_GPCR_Perturbations_up.Adjusted.P.value)), .before = 6)
  
  ## Extract a concatenated string containing genes involved in each GO Biological Process
  GO_PROCESS_GENES = str_split(GO_Processes$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(GO_PROCESS_GENES)){
    GO_PROCESS_GENES[i,] = paste(substr(GO_PROCESS_GENES[i,], 1, 1),
                                 tolower(substr(GO_PROCESS_GENES[i,], 2, 7)), sep = "") 
  }
  ## Extract a concatenated string containing genes involved in each GO Cellular Component
  GO_MOL_FUNC_GENES = str_split(GO_Mol_Funcs$GO_Molecular_Function_2018.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(GO_MOL_FUNC_GENES)){
    GO_MOL_FUNC_GENES[i,] = paste(substr(GO_MOL_FUNC_GENES[i,], 1, 1),
                                  tolower(substr(GO_MOL_FUNC_GENES[i,], 2, 7)), sep = "") 
  }
  ## Extract a concatenated string containing genes involved in each GO Molecular Function
  GO_CELL_COMP_GENES = str_split(GO_Cell_Comps$GO_Cellular_Component_2018.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(GO_CELL_COMP_GENES)){
    GO_CELL_COMP_GENES[i,] = paste(substr(GO_CELL_COMP_GENES[i,], 1, 1),
                                   tolower(substr(GO_CELL_COMP_GENES[i,], 2, 7)), sep = "") 
  }
  ## Extract a concatenated string containing genes involved in each Panther Pathway
  PATHWAY_GENES = str_split(PATHWAYS$Panther_2016.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(PATHWAY_GENES)){
    PATHWAY_GENES[i,] = paste(substr(PATHWAY_GENES[i,], 1, 1),
                              tolower(substr(PATHWAY_GENES[i,], 2, 7)), sep = "") 
  }
  ## Extract a concatenated string containing genes involved in each protein hub
  PPI_GENES = str_split(PPIs$PPI_Hub_Proteins.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(PPI_GENES)){
    PPI_GENES[i,] = paste(substr(PPI_GENES[i,], 1, 1),
                          tolower(substr(PPI_GENES[i,], 2, 7)), sep = "") 
  }
  ## Extract a concatenated string containing genes involved in each disease
  DISEASE_GENES = str_split(DISEASES$Jensen_DISEASES.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(DISEASE_GENES)){
    DISEASE_GENES[i,] = paste(substr(DISEASE_GENES[i,], 1, 1),
                              tolower(substr(DISEASE_GENES[i,], 2, 7)), sep = "") 
  }
  
  ## Extract a concatenated string containing genes by which a TF is identified to possibly regulate
  ENCODE_TF_GENES = str_split(ENCODE_TFs$ENCODE_TF_ChIP.seq_2015.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(ENCODE_TF_GENES)){
    ENCODE_TF_GENES[i,] = paste(substr(ENCODE_TF_GENES[i,], 1, 1),
                                tolower(substr(ENCODE_TF_GENES[i,], 2, 7)), sep = "") 
  }
  
  ## Extract a concatenated string containing genes identified in other submissions relating similar genes to specific TFs
  ENRICHR_TF_GENES = str_split(Enrichr_TFs$Enrichr_Submissions_TF.Gene_Coocurrence.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(ENRICHR_TF_GENES)){
    ENRICHR_TF_GENES[i,] = paste(substr(ENRICHR_TF_GENES[i,], 1, 1),
                                 tolower(substr(ENRICHR_TF_GENES[i,], 2, 7)), sep = "") 
  }
  
  ## Extract a concatenated string containing genes involved in each affected GPCR
  GPCR_DOWN_PERTURB_GENES = str_split(GPCR_Down_Perturbs$L1000_Kinase_and_GPCR_Perturbations_down.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(GPCR_DOWN_PERTURB_GENES)){
    GPCR_DOWN_PERTURB_GENES[i,] = paste(substr(GPCR_DOWN_PERTURB_GENES[i,], 1, 1),
                                        tolower(substr(GPCR_DOWN_PERTURB_GENES[i,], 2, 7)), sep = "") 
  }
  
  ## Extract a concatenated string containing genes involved in each affected GPCR
  GPCR_UP_PERTURB_GENES = str_split(GPCR_Up_Perturbs$L1000_Kinase_and_GPCR_Perturbations_up.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(GPCR_UP_PERTURB_GENES)){
    GPCR_UP_PERTURB_GENES[i,] = paste(substr(GPCR_UP_PERTURB_GENES[i,], 1, 1),
                                      tolower(substr(GPCR_UP_PERTURB_GENES[i,], 2, 7)), sep = "") 
  }
  
  ## Make sure the Enrichr reference results are correctly ordered
  GO_Processes = GO_Processes[order(GO_Processes$Modified_Combined_Score, decreasing = TRUE), ]
  ## Re-set the index to make sure indeces are correct
  row.names(GO_Processes) = NULL
  
  GO_Mol_Funcs = GO_Mol_Funcs[order(GO_Mol_Funcs$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(GO_Mol_Funcs) = NULL
  
  GO_Cell_Comps = GO_Cell_Comps[order(GO_Cell_Comps$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(GO_Cell_Comps) = NULL
  
  PATHWAYS = PATHWAYS[order(PATHWAYS$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(PATHWAYS) = NULL
  
  PPIs = PPIs[order(PPIs$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(PPIs) = NULL
  
  DISEASES = DISEASES[order(DISEASES$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(DISEASES) = NULL
  
  ENCODE_TFs = ENCODE_TFs[order(ENCODE_TFs$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(ENCODE_TFs) = NULL
  
  Enrichr_TFs = Enrichr_TFs[order(Enrichr_TFs$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(Enrichr_TFs) = NULL
  
  GPCR_Down_Perturbs = GPCR_Down_Perturbs[order(GPCR_Down_Perturbs$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(GPCR_Down_Perturbs) = NULL
  
  GPCR_Up_Perturbs = GPCR_Up_Perturbs[order(GPCR_Up_Perturbs$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(GPCR_Up_Perturbs) = NULL
  
  
  ## Export these data frames to the global environment
  GO_Processes <<- GO_Processes
  GO_PROCESS_GENES <<- GO_PROCESS_GENES
  GO_Cell_Comps <<- GO_Cell_Comps
  GO_CELL_COMP_GENES <<- GO_CELL_COMP_GENES
  GO_Mol_Funcs <<- GO_Mol_Funcs
  GO_MOL_FUNC_GENES <<- GO_MOL_FUNC_GENES
  PATHWAYS <<- PATHWAYS
  PATHWAY_GENES <<- PATHWAY_GENES
  PPIs <<- PPIs
  PPI_GENES <<- PPI_GENES
  DISEASES <<- DISEASES
  DISEASE_GENES <<- DISEASE_GENES
  ENCODE_TFs <<- ENCODE_TFs
  ENCODE_TF_GENES <<- ENCODE_TF_GENES
  Enrichr_TFs <<- Enrichr_TFs
  ENRICHR_TF_GENES <<- ENRICHR_TF_GENES
  GPCR_Down_Perturbs <<- GPCR_Down_Perturbs
  GPCR_DOWN_PERTURB_GENES <<- GPCR_DOWN_PERTURB_GENES
  GPCR_Up_Perturbs <<- GPCR_Up_Perturbs
  GPCR_UP_PERTURB_GENES <<- GPCR_UP_PERTURB_GENES
  
  
  
  
  return(c("GO BIOLOGICAL PROCESSES", GO_Processes[1:10,1], "GO CELL COMPONENTS", GO_Cell_Comps[1:10,1], "GO MOLECULAR FUNCTIONS", GO_Mol_Funcs[1:10,1], "PANTHER PATHWAYS", PATHWAYS[1:10,1], "PPIs", PPIs[1:10,1], "JENSEN DISEASE PHENOTYPES", DISEASES[1:10,1]))
}
  ## Perform the analyses
Functional_Analysis(Genes = TmemUnmergedModule)
#Functional_Analysis(Genes = TmemModule)

  ##### Export the Enrichr-based functional analyses of the unmerged TMEM module #####
    ## For the Unmerged module
Functional_Analysis(Genes = TmemUnmergedModule)
#Enrichr_Functional_Output = cbind(GO_Processes[1:20,c(1,4:8,10)],
#                                  GO_Cell_Comps[1:20,c(1,4:8,10)],
#                                  GO_Mol_Funcs[1:20,c(1,4:8,10)],
#                                  PATHWAYS[1:20,c(1,4:8,10)],
#                                  PPIs[1:20,c(1,4:8,10)])
#colnames(Enrichr_Functional_Output) = c("GO Biological Processes",
#                                        "GO Cellular Components",
#                                        "GO Molecular Functions",
#                                        "Panther Pathways",
#                                        "Protein-Protein Interaction Hub Genes")

#write.table(Enrichr_Functional_Output, file = "M:\\Erik\\Data\\Omics\\RNAseq\\Enrichr Analysis on TMEM Unmerged Module.csv", sep = ",", row.names = F, col.names = T)

  ##### Export the Enrichr-based functional analyses of the merged TMEM module #####
Functional_Analysis(Genes = TmemModule)
#Enrichr_Functional_Output = cbind(GO_Processes[1:20,c(1,4:8,10)],
#                                  GO_Cell_Comps[1:20,c(1,4:8,10)],
#                                  GO_Mol_Funcs[1:20,c(1,4:8,10)],
#                                  PATHWAYS[1:20,c(1,4:8,10)],
#                                  PPIs[1:20,c(1,4:8,10)])
#colnames(Enrichr_Functional_Output) = c("GO Biological Processes",
#                                        "GO Cellular Components",
#                                        "GO Molecular Functions",
#                                        "Panther Pathways",
#                                        "Protein-Protein Interaction Hub Genes")

#write.table(Enrichr_Functional_Output, file = "M:\\Erik\\Data\\Omics\\RNAseq\\Enrichr Analysis on TMEM Merged Module.csv", sep = ",", row.names = F, col.names = T)

##### Evaluate clustering by TOM plots #####
  ##### Extract the genes and TOM similarity values within Tmem Modules for Cytoscape analysis #####
    ## Find the "Tmem unmerged module" and its TOM similarity values
      ## The genes
UnmergedTmemTOMnames = c(names(TOM[4395,c(TmemUnmergedModule)]))
      ## The values
UnmergedTmemTOMnums = c(as.numeric(TOM[4395,c(TmemUnmergedModule)]))
    ## Concatenate these vectors into a dataframe for export, visualization in Cytoscape
TmemTOM_Unmerged_Cyto = as.data.frame(cbind(UnmergedTmemTOMnames, UnmergedTmemTOMnums))
colnames(TmemTOM_Unmerged_Cyto)[1] = "GeneID"
colnames(TmemTOM_Unmerged_Cyto)[2] = "TOM Similarity"

    ## Repeat for the merged module
MergedTmemTOMnames = c(names(disTOM[4395,c(TmemModule)]))
MergedTmemTOMnums = c(as.numeric(disTOM[4395,c(TmemModule)]))
TmemTOM_Merged_Cyto = as.data.frame(cbind(MergedTmemTOMnames, MergedTmemTOMnums))
colnames(TmemTOM_Merged_Cyto)[1] = "GeneID"
colnames(TmemTOM_Merged_Cyto)[2] = "TOM Similarity"

    ## Find the clustered genes within the TPM file (expression data) to subset that TPM file (expression data) for more computationally manageable plots
      ## Create a variable to store the indeces of these genes and loop through the list
exprTmemTOMnamesIndeces = vector()
for (i in 1:length(TmemTOMnames)){
  exprTmemTOMnamesIndeces[i] = which(Expression_Profile_cleaned$GeneID == TmemTOMnames[i])
}

    ## Create a "not in" operator
"%ni%" = Negate("%in%")

    ## Find which of the DEGs are within the TmemTOM (unmerged module)
which(DESeq2_Adults3$GeneID[1:405] %in% Expression_Profile_cleaned$GeneID[c(exprTmemTOMnamesIndeces)] == TRUE)
    ## Find the DEGs that are within the TmemTOM unmerged module
which(
  DESeq2_Adults9$GeneID[
    which(DESeq2_Adults9$labs == "DEGs")
    ] %in% Expression_Profile_cleaned$GeneID[c(exprTmemTOMnamesIndeces)] == TRUE)

    ## Find which of the genes within the TmemTOM unmerged module are not DEGs
which(
  Expression_Profile_cleaned$GeneID[
    c(exprTmemTOMnamesIndeces)
    ] %ni% DESeq2_Adults9$GeneID[
      which(DESeq2_Adults9$labs == "DEGs")
      ]
)
    ## Put the genes within the TmemTOM unmerged module that are not DE into a list
Non_DEGs_in_TmemUnmergedModule = c(
  which(
    Expression_Profile_cleaned$GeneID[
      c(exprTmemTOMnamesIndeces)
      ] %ni% DESeq2_Adults9$GeneID[
        which(DESeq2_Adults9$labs == "DEGs")
        ]
  )
)
    ## See which they are
Non_DEGs_in_TmemUnmergedModule = 
  Expression_Profile_cleaned$GeneID[
    c(exprTmemTOMnamesIndeces[
      c(Non_DEGs_in_TmemUnmergedModule)
      ]
      )
    ]
Non_DEGs_in_TmemUnmergedModule

  ##### Evaluate clustering by a smaller TOM plot #####
    ## Create a new TOM based on a subset that comprises the genes with the 405-lowest Q-values, including genes within the TMEM cluster
ShortTOM = disTOM[DESeq2_Adults9$GeneID[1:405], DESeq2_Adults9$GeneID[1:405]]
    ## Concatenate the list of 405 genes
ShortList = c(DESeq2_Adults9$GeneID[1:405])

    ## Turn the subsetted list of genes into a concatenated vector of their indeces within the TPM file (expression dataset)
exprShortTOMnamesIndeces = vector()
for (i in 1:nrow(ShortTOM)){
  exprShortTOMnamesIndeces[i] = which(Expression_Profile_cleaned$GeneID == ShortList[i])
}

    ## Create a new dendrogram within the shortened TPM file (expression dataset) to investigate clustering
shortTree = hclust(as.dist(ShortTOM), method = "average")
    ## Determine clusters to plot based on a "dynamic dendrogram (tree) cut"
unmergedShortLabs = cutreeDynamic(dendro = shortTree, distM = ShortTOM,
                             deepSplit = 3,
                             cutHeight = 0.995,
                             pamRespectsDendro = F)

    ## Convert the clusters from index labels to colors
unmergedShortCols = labels2colors(unmergedShortLabs)
    ## Merge clusters that are similar enough based on the subjective/arbitrary "cutHeight" parameter provided within the WGCNA function, "mergeCloseModules"
mergeShort = mergeCloseModules(t(Expression_Profile_cleaned[exprShortTOMnamesIndeces,2:9]), unmergedShortLabs, cutHeight = 0.01, verbose = 2)
    ## Convert color indeces into a string
moduleShortLabs = mergeShort$colors
    ## Convert the string of labels into colors
moduleShortCols = labels2colors(moduleShortLabs)
    ## Find the new, merged eigenvectors (modules)
ShortMEs = mergeShort$newMEs

    ## Find Tmem index
which(exprShortTOMnamesIndeces == 4395)
    ## Plot the dendrogram and associated clusters/modules
plotDendroAndColors(shortTree, cbind(unmergedShortCols, moduleShortCols), 
                    c("Unmerged", "Merged"), 
                    dendroLabels = F, 
                    hang = 0.001,
                    addGuide = T,
                    guideHang = 0.001,
                    main = "Tmem GT DRG DEG Dendrogram")

#which(Expression_Profile_cleaned$GeneID == "Il31ra")
#which(exprShortTOMnamesIndeces == 1914)

    ## Determine Tmem's modules
moduleShortCols[1]
unmergedShortCols[1]
    ## Find the indeces from the clustered data of genes within Tmem's clusters
shortTmemUnmergedModule = c(which(unmergedShortCols == "turquoise"))
shortTmemModule = c(which(moduleShortCols == "turquoise"))
    ## Find the genes of those indeces
shortTmemUnmergedModule = c(Expression_Profile_cleaned$GeneID[shortTmemUnmergedModule])
shortTmemModule = c(Expression_Profile_cleaned$GeneID[shortTmemModule])

    ## Make the diagonal of the TOM NA to re-scale the plot
diag(ShortTOM) = 0
    ## Create the disTOM plot of the 405 most-differentially expressed genes
TOMplot((ShortTOM^6), shortTree, as.character(unmergedShortCols), main = "disTOM plot Top 405 DEGs")

# ## Create a dendrogram compatible with TOM-TOM plot creation
#colorDynamicShortTOM = labels2colors(cutreeDynamic(shortTree, method ="tree"))
#TOMplot((ShortTOM^6), shortTree, as.character(colorDynamicShortTOM), main = "disTOM plot Top 405 DEGs")

# ## Create an adjacency-based TOMplot
#colorDynamicADJ = labels2colors(cutreeDynamic(shortTree, method = "tree"))
#diag(ADJ) = 1
#TOMplot((ADJ), manTree, as.character(colorDynamicADJ), main = "TOM Plot based on ADJ Sim.")

  ## Create a list of the genes that cluster with Tmem
ShortList[which(unmergedShortCols == "turquoise")]

  ## Perform functional analysis on that cluster
#Functional_Analysis(Genes = ShortList[which(unmergedShortCols == "turquoise")])

##### Export unmerged module for Cytoscape analysis #####
  ## The Unmerged Tmem Module of the 5000-gene net
write.table(TmemTOM_Unmerged_Cyto, file = "M:\\Erik\\Data\\Omics\\RNAseq\\TMEM Unmerged Module Top5000.csv", sep = ",", row.names = F, col.names = T)
##### Export the merged module for Cytoscape analysis #####
  ## The 0.1 CutHeight Merged Tmem Module of the 5000-gene net
write.table(TmemTOM_Merged_Cyto, file = "M:\\Erik\\Data\\Omics\\RNAseq\\TMEM 0.1 CutHeight Merged Module Top5000.csv", sep = ",", row.names = F, col.names = T)

  ##### Export the Enrichr-based functional analysis of the unmerged Tmem module from the shortened net #####
  ## The Unmerged Tmem Module of a shortened net (including the genes within the Tmem module and the top 405 DEGs)
Functional_Analysis(Genes = ShortList[which(unmergedShortCols == "turquoise")])
#Enrichr_Functional_Output = cbind(GO_Processes[1:20,c(1,4:8,10)],
#                                  GO_Cell_Comps[1:20,c(1,4:8,10)],
#                                  GO_Mol_Funcs[1:20,c(1,4:8,10)],
#                                  PATHWAYS[1:20,c(1,4:8,10)],
#                                  PPIs[1:20,c(1,4:8,10)])
#colnames(Enrichr_Functional_Output) = c("GO Biological Processes",
#                                        "GO Cellular Components",
#                                        "GO Molecular Functions",
#                                        "Panther Pathways",
#                                        "Protein-Protein Interaction Hub Genes") 

#write.table(Enrichr_Functional_Output, file = "M:\\Erik\\Data\\Omics\\RNAseq\\TMEM Unmerged Module Top405.csv", sep = ",", row.names = F, col.names = T)



  ##### Export the Enrichr-based functional analysis of the merged Tmem module from the shortened net #####
  ## The 0.1 CutHeight Merged Tmem Module of the shortened net (including the genes within the Tmem module and the top 405 DEGs)
Functional_Analysis(Genes = ShortList[which(moduleShortCols == "turquoise")])
#Enrichr_Functional_Output = cbind(GO_Processes[1:20,c(1,4:8,10)],
#                                  GO_Cell_Comps[1:20,c(1,4:8,10)],
#                                  GO_Mol_Funcs[1:20,c(1,4:8,10)],
#                                  PATHWAYS[1:20,c(1,4:8,10)],
#                                  PPIs[1:20,c(1,4:8,10)])
#colnames(Enrichr_Functional_Output) = c("GO Biological Processes",
#                                        "GO Cellular Components",
#                                        "GO Molecular Functions",
#                                        "Panther Pathways",
#                                        "Protein-Protein Interaction Hub Genes")

#write.table(Enrichr_Functional_Output, file = "M:\\Erik\\Data\\Omics\\RNAseq\\TMEM 0.1 CutHeight Merged Module Top405.csv", sep = ",", row.names = F, col.names = T)



exportNetworkToCytoscape()