




  ## Script for obtaining WGCNA and enrichr-based bioinformatics data
  ## Developed by Erik Larsen

## Condense lines of code by clicking all arrows next to code line #. Click again to open, step-wise to maintain organization

##### Environment Prep #####

  ## Install packages to use the "enrichr", "panther" packages
  ## At least for the first installation, use the following code; THIS REQUIRES R version 4.0.0+ !!!
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

  ## Load the enrichr package
library(enrichR)
  ## Load previously installed packages from the library (if not previously installed, install!)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(plyr)
library(ggrepel)
#library(calibrate)
library(reshape2)
#library(gridExtra)
library(stats)
library(dplyr)
library(readr)
library(readxl)
library(tibble)
library(AnnotationDbi)
library(stats4)
library(BiocGenerics)
library(parallel)
library(dendextend)
library(ggdendro)
library(grid)
#library(EnhancedVolcano)
library(pheatmap)
library(matrixStats)
library(Hmisc)
library(splines)
library(foreach)
library(doParallel)
library(fastcluster)
library(dynamicTreeCut)
library(survival)
library(clValid)
library(WGCNA)
#library(RColorBrewer)
# ("darkgoldenrod4", "white", and "navy" with n = 50 is preferable)

  ## Connect live to the Enrichr server/master database (website) and store the "space" as a variable. This contains a vast array of bioinformatics databases/websites
  ## Store it as a variable for better visualization and eventual subsetting
DBs = listEnrichrDbs()

  ## View the variable in the editor to find the relevant databases and their indeces for subsetting (click on the dataframe in the global environment and manually peruse)

  ## "GO_Biological_Process_2018" (index # 130), and "Panther_2016" (index # 102)
  ## Concatenate them in a list for subsetting, or slice directly
DBs = DBs$libraryName[c(130,102,14,110)]



##### Perform Gene Ontology and Pathway Analysis using DGEA of Tmem184b mutant DRG #####
  ## Import the dataset you want to analyze (Tmem184b-GT/GT DRG DGE)
DESeq2_Adults = read.csv("M:/Erik/Data/Omics/RNAseq/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
DESeq2_Adults3 = subset(DESeq2_Adults, (!is.na(DESeq2_Adults[,"AdjP"])))
#  ## Include Nppb?
#DESeq2_Adults9 = DESeq2_Adults[c(1:12955,25308),]
#  ## Put a placeholder Q value in for now
#DESeq2_Adults9[12956,7] = 0.9999999
  ## Re-order the index
#rownames(DESeq2_Adults9) = NULL

  ## Perform relevant analysis through enrichr; in this case, FDR <= 0.01
  ## first: GO Bio Processes; this analyzes genes from your dataset in enrichr, which draws terms from geneontology.org

  ## Define the input for high dimensional analysis
Itch_Genes = c(DESeq2_Adults3$GeneID[which(DESeq2_Adults3$AdjP <= 0.05)])
  ## Define the global environment dataframes to fill
GOs = data.frame()
PATHWAYS = data.frame()
PPIs = data.frame()
GO_GENES = data.frame()
PATHWAY_GENES = data.frame()
PPI_GENES = data.frame()
DISEASES = data.frame()
DISEASE_GENES = data.frame()


  ## Create a function that will return the WT Exp profile of 4 C57Bl6 mouse DRG, normalized by gene across replicates and by sample across mice (Z-score) to enable a Euclidean clustering algorithm to organize expression by related clusters and permit Weighted Gene Co-expression Analysis (WGCNA).
Functional_Analysis = function(Genes){
  GOs = as.data.frame(
    enrichr(
      c(Genes), DBs[1])
  )
    ## second: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
  PATHWAYS = as.data.frame(
    enrichr(
      c(Genes), DBs[2])
  )
    ## third: PPIs; this identifies major "hub" genes' protein products- what they interact with, etc.
  PPIs = as.data.frame(
    enrichr(
      c(Genes), DBs[3])
  )
    ## last: diseases; this identifies phenotypes or diseases associated with the genes provided
  DISEASES = as.data.frame(
    enrichr(
      c(Genes), DBs[4])
  )
  
    ## Remove unnecessary columns
  GOs = GOs[,-c(5,6,7)]
  PATHWAYS = PATHWAYS[,-c(5,6,7)]
  PPIs = PPIs[,-c(5,6,7)]
  DISEASES = DISEASES[,-c(5,6,7)]
  
    ## Clean up both dataframes and export both dataframes.
    ## The pathway analysis dataframe needs additional cleaning in Python;
    ## Following that manipulation, data will be re-imported in other scripts for visualization using ggplot
  GOs = GOs %>%
    separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
    ## Convert the values into integers for computation
  GOs$Overlap = as.numeric(GOs$Overlap)
  GOs$`Process Size` = as.numeric(GOs$`Process Size`)
  
    ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  GOs = add_column(GOs, EnrichrZscore = GOs$GO_Biological_Process_2018.Combined.Score/log(GOs$GO_Biological_Process_2018.P.value), .before = 6)
    ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
      ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  GOs = add_column(GOs, Weighted_Overlap_Ratio = GOs$Overlap*(GOs$Overlap/GOs$`Process Size`), .before = 4)
    ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  GOs = add_column(GOs, Modified_Combined_Score = (abs(GOs$EnrichrZscore))*(-log(GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)
    ## Remove GO IDs
  GOs$GO_Biological_Process_2018.Term = GOs$GO_Biological_Process_2018.Term %>% gsub(x = GOs$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")
  

    ## Repeat the above code for Panther Pathways
  PATHWAYS = PATHWAYS %>%
    separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
    ## Convert the values into integers for computation
  PATHWAYS$Overlap = as.numeric(PATHWAYS$Overlap)
  PATHWAYS$`Process Size` = as.numeric(PATHWAYS$`Process Size`)
  
    ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
  PATHWAYS = add_column(PATHWAYS, EnrichrZscore = PATHWAYS$Panther_2016.Combined.Score/log(PATHWAYS$Panther_2016.P.value), .before = 6)
    ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
  PATHWAYS = add_column(PATHWAYS, Weighted_Overlap_Ratio = PATHWAYS$Overlap*(PATHWAYS$Overlap/PATHWAYS$`Process Size`), .before = 4)
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
  
  
    ## Extract a concatenated string containing genes involved in each GO Biological Process
  GO_GENES = str_split(GOs$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(GO_GENES)){
    GO_GENES[i,] = paste(substr(GO_GENES[i,], 1, 1),
                               tolower(substr(GO_GENES[i,], 2, 7)), sep = "") 
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
  
  
  
    ## Make sure the Enrichr reference results are correctly ordered
  GOs = GOs[order(GOs$Modified_Combined_Score, decreasing = TRUE), ]
    ## Re-set the index to make sure indeces are correct
  row.names(GOs) = NULL
  
  PATHWAYS = PATHWAYS[order(PATHWAYS$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(PATHWAYS) = NULL
  
  PPIs = PPIs[order(PPIs$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(PPIs) = NULL
  
  DISEASES = DISEASES[order(DISEASES$Modified_Combined_Score, decreasing = TRUE), ]
  row.names(DISEASES) = NULL
  
  GOs <<- GOs
  GO_GENES <<- GO_GENES
  PATHWAYS <<- PATHWAYS
  PATHWAY_GENES <<- PATHWAY_GENES
  PPIs <<- PPIs
  PPI_GENES <<- PPI_GENES
  DISEASES <<- DISEASES
  DISEASE_GENES <<- DISEASE_GENES
  return(c(GOs[1:10,1], PATHWAYS[1:10,1], PPIs[1:10,1], DISEASES[1:10,1]))
}
Functional_Analysis(Genes = Itch_Genes)


##### Analyze the mutant Tmem184b gene network by creating heatmaps (Euclidean Hierarchical clustering) #####
##### Data prep #####
  ## Import the normalized counts file. These are expression estimates for each gene, for each sample/replicate, where each gene's value is normalized to its sample's effect size
Adult_Normalized_Counts = read_csv("M:/Erik/Data/Omics/RNAseq/Processed Galaxy Output/Counts Files to Upload/RNASeqRepResults.csv", col_names = TRUE)
  ## Rename columns; should know this ahead of time
colnames(Adult_Normalized_Counts) = c("GeneID", "WT1", "WT2", "WT3", "WT4", "Mut1", "Mut2", "Mut3", "Mut4")

  ## Subset the Normalized Counts file by genes from the DGEA results dataframe containing only genes with Adjusted P-values
RNASeqRepResultsAdultAll = Adult_Normalized_Counts[Adult_Normalized_Counts$GeneID %in% c(DESeq2_Adults3$GeneID),]
  ## Filter rRNAs or any other unidentified mRNAs
RNASeqRepResultsAdultAll = RNASeqRepResultsAdultAll %>% filter(!grepl(RNASeqRepResultsAdultAll$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

  ## For DEGs
DESeq2_AdultsItchDEGs = DESeq2_Adults3[c(1:405),]
    ## Filter rRNAs or any other unidentified mRNAs
DESeq2_AdultsItchDEGs = DESeq2_AdultsItchDEGs %>% filter(!grepl(DESeq2_AdultsItchDEGs$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

  ## Set up potential dataframes for downstream analyses (WT, Mut)
  ## Duplicate the entire transcriptional profile
AdultItchExp = RNASeqRepResultsAdultAll
  ## Prep the Mutant subset
MutAdultItchExp = RNASeqRepResultsAdultAll[RNASeqRepResultsAdultAll$GeneID %in% c(DESeq2_AdultsItchDEGs$GeneID),]
  ## Slice it
MutAdultItchExp = MutAdultItchExp[,c(1,6:9)]
  ## Slice for WT
WTAdultItchExp = AdultItchExp[,c(1:5)]

##### Find Mutant Z-scores and graph the profile #####
  ## Create a function that determines the Zscores of the expression profile
Find_Z = function(Expression_Profile){
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
  ## Run the function
Find_Z(Expression_Profile = MutAdultItchExp)


  ## Create a list of gene names ordered by euclidean distance by which to re-order the Z-scored dataframe
Mut_Euclid_dist_order_Z = hclust(dist(Z[,c(2:5)], method = "euclidean"))$order
  ## The names (not the numbers)
Mut_Euclid_dist_ord_Z_Genes = c(Z$GeneID[Mut_Euclid_dist_order_Z])


  ## Transform to order by clusters (Euclidean distances)
Mut_Itch_Corr_Z = Z %>%
  mutate(GeneID =  factor(GeneID, levels = Mut_Euclid_dist_ord_Z_Genes)) %>%
  arrange(GeneID)

  ## Find relevant reference genes to verify clustering and provide reference indeces for plotting
which(Mut_Euclid_dist_ord_Z_Genes == "Il31ra") # 227
which(Mut_Euclid_dist_ord_Z_Genes == "Cysltr2") # 186
which(Mut_Euclid_dist_ord_Z_Genes == "Npy2r") # 323
which(Mut_Euclid_dist_ord_Z_Genes == "Sst") # 263
which(Mut_Euclid_dist_ord_Z_Genes == "Htr1a") # 338
which(Mut_Euclid_dist_ord_Z_Genes == "P2rx3") # 101
which(Mut_Euclid_dist_ord_Z_Genes == "Lpar3") # 79
which(Mut_Euclid_dist_ord_Z_Genes == "Lpar5") # 346
which(Mut_Euclid_dist_ord_Z_Genes == "Scn11a") # 76
which(Mut_Euclid_dist_ord_Z_Genes == "Scn10a") # 58
which(Mut_Euclid_dist_ord_Z_Genes == "Mrgprd") # 103
which(Mut_Euclid_dist_ord_Z_Genes == "Trpc6") # 113
which(Mut_Euclid_dist_ord_Z_Genes == "Trpc3") # 85
which(Mut_Euclid_dist_ord_Z_Genes == "F2rl2") # 29
which(Mut_Euclid_dist_ord_Z_Genes == "Htr1f") # 367
which(Mut_Euclid_dist_ord_Z_Genes == "Osmr") # 145
which(Mut_Euclid_dist_ord_Z_Genes == "Fxyd2") # 173
which(Mut_Euclid_dist_ord_Z_Genes == "Htr4") # 138
which(Mut_Euclid_dist_ord_Z_Genes == "Mrgprx1") # 325
which(Mut_Euclid_dist_ord_Z_Genes == "Ptgdr") # 120
which(Mut_Euclid_dist_ord_Z_Genes == "Trpa1") # 353
which(Mut_Euclid_dist_ord_Z_Genes == "Trpm6") # 84
which(Mut_Euclid_dist_ord_Z_Genes == "Hrh1") # 224
which(Mut_Euclid_dist_ord_Z_Genes == "Mrgpra3") # 57
#which(Mut_Euclid_dist_ord_Genes == "Nppb") # 134
which(Mut_Euclid_dist_ord_Z_Genes == "Tmem184b") # 233
which(Mut_Euclid_dist_ord_Z_Genes == "Noxo1") # 86


  ## Create a list of gene names ordered by euclidean distance by which to re-order the dataframe
Mut_Euclid_dist_order = hclust(dist(MutAdultItchExp[-219,c(2:5)], method = "euclidean"))$order
  ## The names (not the numbers)
Mut_Euclid_dist_ord_Genes = c(MutAdultItchExp$GeneID[Mut_Euclid_dist_order])


  ## Transform to order by clusters (Euclidean distances)
Mut_Itch_Corr = MutAdultItchExp %>%
  mutate(GeneID =  factor(GeneID, levels = Mut_Euclid_dist_ord_Genes)) %>%
  arrange(GeneID)

  ## Find relevant reference genes to verify clustering and provide reference indeces for plotting
which(Mut_Euclid_dist_ord_Genes == "Il31ra") # 147
which(Mut_Euclid_dist_ord_Genes == "Cysltr2") # 145
which(Mut_Euclid_dist_ord_Genes == "Npy2r") # 183
which(Mut_Euclid_dist_ord_Genes == "Sst") # 138
which(Mut_Euclid_dist_ord_Genes == "Htr1a") # 148
which(Mut_Euclid_dist_ord_Genes == "P2rx3") # 222
which(Mut_Euclid_dist_ord_Genes == "Lpar3") # 261
which(Mut_Euclid_dist_ord_Genes == "Lpar5") # 172
which(Mut_Euclid_dist_ord_Genes == "Scn11a") # 366
which(Mut_Euclid_dist_ord_Genes == "Scn10a") # 7
which(Mut_Euclid_dist_ord_Genes == "Mrgprd") # 26
which(Mut_Euclid_dist_ord_Genes == "Trpc6") # 208
which(Mut_Euclid_dist_ord_Genes == "Trpc3") # 231
which(Mut_Euclid_dist_ord_Genes == "F2rl2") # 234
which(Mut_Euclid_dist_ord_Genes == "Htr1f") # 154
which(Mut_Euclid_dist_ord_Genes == "Osmr") # 143
which(Mut_Euclid_dist_ord_Genes == "Fxyd2") # 362
which(Mut_Euclid_dist_ord_Genes == "Htr4") # 224
which(Mut_Euclid_dist_ord_Genes == "Mrgprx1") # 185
which(Mut_Euclid_dist_ord_Genes == "Ptgdr") # 115
which(Mut_Euclid_dist_ord_Genes == "Trpa1") # 228
which(Mut_Euclid_dist_ord_Genes == "Trpm6") # 142
which(Mut_Euclid_dist_ord_Genes == "Hrh1") # 146
which(Mut_Euclid_dist_ord_Genes == "Mrgpra3") # 219
#which(Mut_Euclid_dist_ord_Genes == "Nppb") # 
which(Mut_Euclid_dist_ord_Genes == "Tmem184b") # 164
which(Mut_Euclid_dist_ord_Genes == "Noxo1") # 161


## Confirm by indexing; copy and paste to customize the clustered heatmap
Mut_Euclid_dist_ord_Z_Genes[c(227,186,323,263,338,101,79,346,76,58,103,113,85,29,367,145,173,138,325,120,353,84,224,57,233,86)]
## Confirm the location/order of the transformed matrix/dataframe
## Should be "Il31ra"
Mut_Itch_Corr_Z[c(227,186,323,263,338,101,79,346,76,58,103,113,85,29,367,145,173,138,325,120,353,84,224,57,233,86),]



  ## Confirm by indexing; copy and paste to customize the clustered heatmap
Mut_Euclid_dist_Genes[c(147,145,183,138,148,222,261,172,366,7,26,208,231,234,154,143,362,224,185,115,228,142,146,219,164,161)]
  ## Confirm the location/order of the transformed matrix/dataframe
  ## Should be "Il31ra"
Mut_Itch_Corr[c(147,145,183,138,148,222,261,172,366,7,26,208,231,234,154,143,362,224,185,115,228,142,146,219,164,161),]


  ## Visualize the Z-scored, Euclidean-clustered and ordered DEG TPMs across mutants with "pheatmap"
pheatmap(mat = Mut_Itch_Corr_Z[1:375,2:5], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, show_rownames = F)
  ## Compare to un-Z-scored, Euclidean-clustered and ordered DEG TPMs across mutants
pheatmap(mat = Mut_Itch_Corr[-161,2:5], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, show_rownames = F)
  ## --> Obscured by genes at the bottom of the dataframe that comprise the maximal transcriptional values

  ## Visualize the Tmem184b cluster of Z-scored data
pheatmap(mat = Mut_Itch_Corr[198:248,2:5], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, labels_row = Mut_Euclid_dist_ord_Genes[198:248])


##### Analyze the normal Tmem184b gene network by Pearson's Correlation on the Z-scored and Euclidean clustered data #####
  ## Z-score the WT profiles
Find_Z(Expression_Profile = AdultItchExp)
  ## Conduct a Pearson correlation on the Z-scored WT samples
#Itch_Pear = cor(x = Z[,2:9], method = "pearson", use = "p")
  ## Create a new dataframe that clusters the Z-scored WT expression data based on Euclidean distance (basically Pearson correlation)
Itch_Clust = hclust(dist(Z[,c(2:9)], method = "euclidean"))$order
  ## The names (not the numbers)
Clust_Genes = c(Z$GeneID[Itch_Clust])

  ## Transform to order by clusters (Euclidean distances)
Itch_Corr = Z %>%
  mutate(GeneID =  factor(GeneID, levels = Clust_Genes)) %>%
  arrange(GeneID)

  ## Find relevant reference genes to verify clustering and provide reference indeces for plotting
which(Clust_Genes == "Il31ra") # 1614
which(Clust_Genes == "Cysltr2") # 1613
which(Clust_Genes == "Npy2r") # 1947
which(Clust_Genes == "Sst") # 1623
which(Clust_Genes == "Htr1a") # 1616
which(Clust_Genes == "P2rx3") # 1307
which(Clust_Genes == "Lpar3") # 1871
which(Clust_Genes == "Lpar5") # 1689
which(Clust_Genes == "Scn11a") # 1883
which(Clust_Genes == "Scn10a") # 1893
which(Clust_Genes == "Mrgprd") # 1882
which(Clust_Genes == "Trpc6") # 2427
which(Clust_Genes == "Trpc3") # 1543
which(Clust_Genes == "F2rl2") # 1901
which(Clust_Genes == "Htr1f") # 1696
which(Clust_Genes == "Osmr") # 1596
which(Clust_Genes == "Fxyd2") # 2023
which(Clust_Genes == "Htr4") # 1785
which(Clust_Genes == "Mrgprx1") # 1702
which(Clust_Genes == "Ptgdr") # 2431
which(Clust_Genes == "Trpa1") # 1698
which(Clust_Genes == "Trpm6") # 1612
which(Clust_Genes == "Hrh1") # 1546
which(Clust_Genes == "Mrgpra3") # 2061
#which(WT_Clust_Genes == "Nppb") # 
which(Clust_Genes == "Tmem184b") # 1617
which(Clust_Genes == "Noxo1") # 4306


  ## Confirm by indexing; copy and paste to customize the clustered heatmap
Clust_Genes[c(1614,1613,1947,1623,1616,1307,1871,1689,1883,1893,1882,2427,1543,1901,1696,1596,2023,1785,1702,2431,1698,1612,1546,2061,1617,4306)]
  ## Confirm the location/order of the transformed matrix/dataframe
    ## Should be "Il31ra"
Itch_Corr[c(1614,1613,1947,1623,1616,1307,1871,1689,1883,1893,1882,2427,1543,1901,1696,1596,2023,1785,1702,2431,1698,1612,1546,2061,1617,4306),]


  ## Visualize the Z-scored, Euclidean-clustered and ordered TPMs across WTs with "pheatmap"
pheatmap(mat = Itch_Corr[,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, show_rownames = F)
  ## Visualize around the Tmem184b cluster (ballparked)
pheatmap(mat = Itch_Corr[1600:1650,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, labels_row = Clust_Genes[1600:1650])

#Z[5457:5558,]


  ## Perform gene ontology and pathway analysis of the Tmem184b cluster (ballparked)
Functional_Analysis(Genes = c(Clust_Genes[1600:1700]))



##### Hack the WGCNA Tutorial. Determine some parameters.. #####
  ## Hack the WGCNA Tutorial. Determine the soft thresholding power for "inclusivity" within eigengene modules
  ## Create a data structure to house the data for downstream manipulation and analysis
multiExpr = vector(mode = "list", length = 1)
  # Create a vector as long as the expression dataframe for determining additional genes to remove upstream of correlation calculations
expcheck = as.array(vector(length = nrow(WTAdultItchExp)))

Remove_NAN_Generators = function(Expression_Profile, Samples_to_Compare){
  
  ## Classify which analysis to determine dataframe subset
  
  #MutAdultItchExp = RNASeqRepResultsAdultAll[,grep(names(RNASeqRepResultsAdultAll), pattern = "GeneID|Mut")]
  #WTAdultItchExp = RNASeqRepResultsAdultAll %>% select(starts_with(GeneID) & ends_with(WT4))
  
  if (Samples_to_Compare == "WT" ){
    Expression_Profile = Expression_Profile[,grep(names(Expression_Profile), pattern = "GeneID|WT")]
  } else if (Samples_to_Compare == "Mut"){
    Expression_Profile = Expression_Profile[,grep(names(Expression_Profile), pattern = "GeneID|Mut")]
  } else if (Samples_to_Compare == "All"){
    EXpression_Profile = Expression_Profile
  } else {
    Expression_Profile = Expression_Profile
  }
  ## Run a for-loop to evaluate whether any genes in the dataframe contain all 0s. Store it in the above initialized vector
    ## For an entire profile
  if (ncol(Expression_Profile) > 5){
    for (i in 1:nrow(Expression_Profile)){
      expcheck[i] = if_else(Expression_Profile[i,2] == 0 & Expression_Profile[i,3] == 0 & Expression_Profile[i,4] == 0 & Expression_Profile[i,5] == 0, Expression_Profile[i,6] == 0 & Expression_Profile[i,7] == 0 & Expression_Profile[i,8] == 0 & Expression_Profile[i,9] == 0, "Remove", "Keep")
    }
  } 
   ## For Genotype-specific profile
    else {
    for (i in 1:nrow(Expression_Profile)){
      expcheck[i] = if_else(Expression_Profile[i,2] == 0 & Expression_Profile[i,3] == 0 & Expression_Profile[i,4] == 0 & Expression_Profile[i,5] == 0, "Remove", "Keep")
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
    Expression_Profile_cleaned = Expression_Profile
  }else if (Samples_to_Compare == "All") {
    Expression_Profile_cleaned = Expression_Profile
  }else{
    Expression_Profile_cleaned = Expression_Profile
  }
  
  Expression_Profile_cleaned <<- Expression_Profile_cleaned
}
#  ## Run the function
Remove_NAN_Generators(Expression_Profile = AdultItchExp, Samples_to_Compare = "WT")


    ## Fill the data structure with the transposed Adult DRG WT Expression profile
multiExpr[[1]] = list(data = as.data.frame(t(Expression_Profile_cleaned[,2:5])))
    ## Rename the columns by their GeneIDs
names(multiExpr[[1]]$data) = Expression_Profile_cleaned$GeneID
    ## Rename the rows by their Sample IDs
rownames(multiExpr[[1]]$data) = names(Expression_Profile_cleaned[2:5])
  
  ## Take a range of potential powers to determine the appropriate power by which the correlation matrix will be raised. This can be   thought of as increasing the "connectivity" of more related genes
  
powers = c(seq(4,10,by=1), seq(12,20, by=2))
    ## Use the WGCNA function, "pickSoftThreshold", to compute model fits of various powers, and to compute the "connectivity" metrics of the data based on these power values
powerTable = list(data = pickSoftThreshold(data = multiExpr[[1]]$data, powerVector = powers, verbose = 2))
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
  plot(Topology_df$Power, Topology_df$slope*Topology_df$SFT.R.sq,
       xlab = "Soft Threshold (power)",
       ylab = colNames[1], type = "n", 
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(0,1), main = colNames[1])
  addGrid()
  text(Topology_df$Power, Topology_df$slope*Topology_df$SFT.R.sq, labels = powers)
  
    ## Graph the mean connectivity as a function of soft threshold (power)
  plot(Topology_df$mean.k., Topology_df$slope*Topology_df$mean.k.,
       xlab = "Soft Threshold (power)",
       ylab = colNames[2], type = "n",
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(min(Topology_df$slope*Topology_df$mean.k.), max(Topology_df$slope*Topology_df$mean.k.)), main = colNames[2])
  addGrid()
  text(Topology_df$Power, (Topology_df$slope)*Topology_df$mean.k., labels = powers)
  
    ## Graph the median connectivity as a function of soft threshold (power)
  plot(Topology_df$median.k., Topology_df$slope*Topology_df$median.k.,
       xlab = "Soft Threshold (power)", ylab = colNames[3], type = "n",
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(min(Topology_df$slope*Topology_df$median.k.), max(Topology_df$slope*Topology_df$median.k.)), main = colNames[3])
  addGrid()
  text(Topology_df$Power, (Topology_df$slope)*Topology_df$median.k., labels = powers)
  
    ## Graph the maximum connectivity as a function of soft threshold (power)
  plot(Topology_df$max.k., -Topology_df$slope*Topology_df$max.k.,
       xlab = "Soft Threshold (power)", ylab = colNames[4], type = "n",
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(min(Topology_df$slope*Topology_df$max.k.), max(Topology_df$slope*Topology_df$max.k.)), main = colNames[4])
  addGrid()
  text(Topology_df$Power, (Topology_df$slope)*Topology_df$max.k., labels = powers)
  ### -->> beta value of 18  seems optimal in the over-fit model ###
  
  WTnet = blockwiseConsensusModules(multiExpr = multiExpr,
                                  checkMissingData = T,
                                  maxBlockSize = 11226,
                                  blockSizePenaltyPower = Inf,
                                  power = 6,
                                  networkType = "signed",
                                  deepSplit = 1,
                                  detectCutHeight = 0.995,
                                  minModuleSize = 30,
                                  pamRespectsDendro = T,
                                  mergeCutHeight = 0.10,
                                  numericLabels = F,
                                  minCoreKME = 0.7,
                                  minCoreKMESize = 30,
                                  minKMEtoStay = 0.6,
                                  saveTOMs = F,
                                  verbose = 4)
  

  WTModuleEigengenes = c(names(WTnet$multiMEs[[1]]$data))
  ## Create a variable for the merged eigenvectors for plotting
  WTModuleLabels = WTnet$colors
  ## Create a variable of the dendrogram for plotting
  WTGeneDendrogram = WTnet$dendrograms[[1]]
  ## Create a variable (vector) of colors corresponding to the eigenvectors
  WTModuleColors = labels2colors(WTModuleLabels)
  
  ## Find Tmem's index to find it in the net
  which(RNASeqRepResultsAdultAll$GeneID == "Tmem184b")
  ## Use that index to find Tmem in the network
  WTnet$colors[9726]
  ## Find Tmem in the pre-merged network
  WTnet$unmergedColors[9726]
  ## Find the Tmem module
  WTTmem_module = c(which(WTnet$colors == "midnightblue"))
  ## Remove the indeces/numbers for downstream analysis
  WTTmem_module = c(names(WTTmem_module))
  ## Check out the module
  WTTmem_module
  ## Perform downstream bioinformatics (GO, pathway, PPI analysis)
  Functional_Analysis(Genes = WTTmem_module)
  
  
  ## Plot the dendrogram
  plotDendroAndColors(dendro = WTGeneDendrogram, colors = WTModuleColors, "Module Colors", dendroLabels = F, hang = 0.01, addGuide = T, guideHang = 0.01, main = "Gene Dendrogram and Module Colors", marAll = c(4,4,6,4))
  
  
  

###### Build a hierarchically clustered dendrogram for the entire network using the WGCNA package "blockwiseConsensusModules ######

    ## Initialize the dataset
  multiExpr = vector(mode = "list", length = 1)
  multiExpr[[1]] = list(data = as.data.frame(t(AdultItchExp[,2:9])))
    ## Rename the columns by their GeneIDs
  names(multiExpr[[1]]$data) = AdultItchExp$GeneID
    ## Rename the rows by their Sample IDs
  rownames(multiExpr[[1]]$data) = names(AdultItchExp[,2:9])
  
    ## Determine the power by which to multiply the adjacency matrix
  powers = c(seq(4,10,by=1), seq(12,20, by=2))
  powerTable = list(data = pickSoftThreshold(data = multiExpr[[1]]$data, powerVector = powers, verbose = 2))
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
       ylim = c(0,4), main = colNames[1])
  addGrid()
  text(Topology_df$Power, -Topology_df$slope*Topology_df$SFT.R.sq, labels = powers)
  
  ## Graph the mean connectivity as a function of soft threshold (power)
  plot(Topology_df$mean.k., -Topology_df$slope*Topology_df$mean.k.,
       xlab = "Soft Threshold (power)",
       ylab = colNames[2], type = "n",
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(min(-Topology_df$slope*Topology_df$mean.k.), max(-Topology_df$slope*Topology_df$mean.k.)), main = colNames[2])
  addGrid()
  text(Topology_df$Power, -Topology_df$slope*Topology_df$mean.k., labels = powers)
  
  ## Graph the median connectivity as a function of soft threshold (power)
  plot(Topology_df$median.k., Topology_df$slope*Topology_df$median.k.,
       xlab = "Soft Threshold (power)", ylab = colNames[3], type = "n",
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(min(-Topology_df$slope*Topology_df$median.k.), max(-Topology_df$slope*Topology_df$median.k.)), main = colNames[3])
  addGrid()
  text(Topology_df$Power, -Topology_df$slope*Topology_df$median.k., labels = powers)
  
  ## Graph the maximum connectivity as a function of soft threshold (power)
  plot(Topology_df$max.k., -Topology_df$slope*Topology_df$max.k.,
       xlab = "Soft Threshold (power)", ylab = colNames[4], type = "n",
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(min(-Topology_df$slope*Topology_df$max.k.), max(-Topology_df$slope*Topology_df$max.k.)), main = colNames[4])
  addGrid()
  text(Topology_df$Power, -Topology_df$slope*Topology_df$max.k., labels = powers)
  
    ## Confirm the optimized power value
  powerTable$data[1]
  ### -->> beta value of 6 (index #3; covariances raised to power of 6) seems optimal ###
  
net = blockwiseConsensusModules(multiExpr = multiExpr,
                                checkMissingData = T,
                                maxBlockSize = 11226,
                                blockSizePenaltyPower = Inf,
                                corType = "pearson",
                                power = 6,
                                networkType = "signed",
                                TOMtype = "signed",
                                deepSplit = 2,
                                detectCutHeight = 0.995,
                                minModuleSize = 30,
                                pamRespectsDendro = T,
                                mergeCutHeight = 0.2,
                                numericLabels = F,
                                minCoreKME = 0.8,
                                minCoreKMESize = 10,
                                minKMEtoStay = 0.7,
                                saveTOMs = F,
                                verbose = 4)


  ## Create a variable for the eigengenes for plotting
ModuleEigengenes = c(names(net$multiMEs[[1]]$data))
  ## Create a variable for the merged eigenvectors for plotting
ModuleLabels = net$colors
  ## Create a variable of the dendrogram for plotting
GeneDendrogram = net$dendrograms[[1]]
  ## Create a variable (vector) of colors corresponding to the eigenvectors
ModuleColors = labels2colors(ModuleLabels)

  ## Find Tmem's index to find it in the net
which(RNASeqRepResultsAdultAll$GeneID == "Tmem184b")
  ## Use that index to find Tmem in the network
net$colors[9727]
  ## Find Tmem in the pre-merged network
net$unmergedColors[9727]
  ## Find the Tmem module
Tmem_module = c(which(net$colors == "turquoise"))
  ## Remove the indeces/numbers for downstream analysis
Tmem_module = c(names(Tmem_module))
  ## Check out the module
Tmem_module
  ## Perform downstream bioinformatics (GO, pathway, PPI analysis)
Functional_Analysis(Genes = Tmem_module)


  ## Plot the dendrogram
plotDendroAndColors(dendro = GeneDendrogram, colors = ModuleColors, "Module Colors", dendroLabels = F, hang = 0.01, addGuide = T, guideHang = 0.01, main = "Tmem184b GT DRG Gene Dendrogram", marAll = c(4,4,6,4))





##### Clustering validation ... In progress #####
  ## Evaluate the hierarchical clustering
a = clValid(obj = as.data.frame(t(multiExpr[[1]]$data)), nClust = 131, clMethods = "hierarchical", validation = "internal", method = "average", maxitems = nrow(t(multiExpr[[1]]$data)))

  ## Dunn's index should be high, although given the noise in the dataset, low values are not surprising; could be overfit, or relationships might be too "blurry". The network is weighted and signed..
a@measures[]

#dunn(distance = NULL, clusters = net$colors)

  ## Create a network from the ground, up
  ## Create the adjacency matrix, which houses the co-relational (similarity/distance) measures of genes between genes.
    ## Raise to the power that maximizes the scale free topology model fit and the mean or median connectivity (6 in this case)
adjacencyM = cor(multiExpr[[1]]$data, use = "p")^6
  ## Convert the adjacency to a distance/DISsimilarity matrix
distTOMs = TOMsimilarity(adjacencyM)
  ## Cluster the data for a dendrogram
Tree = hclust(as.dist(distTOMs), method = "average")
  ## Duplicate the dissimilarity matrix for investigation of values within a slice
newTOMs = distTOMs
  ## Now find the minimally distanced genes related to e.g. Tmem
which(newTOMs[1,2] == min(newTOMs[9727,-9727]))
which(newTOMs[,] == max(newTOMs[9727,-9727]))


which(newTOMs[,] == 0.2237042)
  ## Divide results by the length or width of the matrix to obtain values
109194946/nrow(newTOMs)-1
122025121/nrow(newTOMs)-1



which(newTOMs[,] == )

13851385/nrow(newTOMs)
109185310/nrow(newTOMs)
  ## Round up and search, by index, the dataframe by which the distances were originally calculated
AdultItchExp[10869,1]
AdultItchExp[9727,1]
AdultItchExp[1235,1]

      ## The same can be done for any row/column in the entire matrix
  ## Slice the Tmem adjacency
Tmem_adjacencies = newTOMs[9727,]
  ## Arbitrarily pick genes most closely related to Tmem
which(Tmem_adjacencies < 0.01)
  ## Concatenate them into a list of vector indeces even more closely related to Tmem (pare down the list above)
Tmem_adjacencies2 = which(Tmem_adjacencies < 0.007)
  ## Plug the indeces into the original dataframe to find which genes they are, and concatenate those genes into a list
Tmemmodulegenes = c(AdultItchExp[Tmem_adjacencies2,1])
  ## Peruse the list
Tmemmodulegenes$GeneID[]
Tmemmodulegenes$GeneID[1001:2000]
Tmemmodulegenes$GeneID[2001:2283]

  ## Prep a plot of this WGCNA method
    ## Create a variable of labels; use the dissimilarity distance matrix (1 - TOM) to make the dendrogram with the same base parameters as the blockwiseConsensusModule function
unmergedLabels = cutreeDynamic(dendro = Tree, distM = distTOMs, deepSplit = 4, minClusterSize = 25, pamRespectsDendro = T, method = "hybrid", respectSmallClusters = T, cutHeight = 0.25)
    ## Create a variable of the colors accompanying those labels
unmergedColors = labels2colors(unmergedLabels)
    ## Peek at the labels (eigenvectors/clusters/modules)
table(unmergedLabels)
    ## Find Tmem's index
which(unmergedLabels == 9727)
    ## Find Tmeme's color module
Tree$Colors[]

sizeGrWindow(8,6)
    ## Plot
plotDendroAndColors(Tree, unmergedColors, "Dynamic Tree Cut", dendroLabels = F, hang = 0.01, addGuide = T, guideHang = 0.03)




consMEsC = multiSetMEs(multiExpr, universalColors = ModuleColors)
MET = consensusOrderMEs(consMEsC)

sizeGrWindow(8,10)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels = , marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)