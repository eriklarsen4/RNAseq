




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
library(clValid)
library(WGCNA)
library(bigmemory)
library(anRichmentMethods)

#library(RColorBrewer)
# ("darkgoldenrod4", "white", and "navy" with n = 50 is preferable)

  ## Connect live to the Enrichr server/master database (website) and store the "space" as a variable. This contains a vast array of bioinformatics databases/websites.
  ## Store it as a variable for better visualization and eventual subsetting
DBs = listEnrichrDbs()

  ## View the variable in the editor to find the relevant databases and their indeces for subsetting (click on the dataframe in the Global Environment and manually peruse)

  ## e.g. "GO_Biological_Process_2018" (index # 130), and "Panther_2016" (index # 102)
  ## Concatenate them in a list for subsetting, or slice directly
DBs = DBs$libraryName[c(130,131,132,102,14,110)]



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
GO_Processes = data.frame()
GO_Cell_Comps = data.frame()
GO_Mol_Funcs = data.frame()
PATHWAYS = data.frame()
PPIs = data.frame()
DISEASES = data.frame()
GO_PROCESS_GENES = data.frame()
GO_MOL_FUNC_GENES = data.frame()
GO_CELL_COMP_GENES = data.frame()
PATHWAY_GENES = data.frame()
PPI_GENES = data.frame()
DISEASE_GENES = data.frame()


  ## Create a function that will return functional analysis on a given gene set.
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
  
    ## Remove unnecessary columns
  GO_Processes = GO_Processes[,-c(5,6,7)]
  GO_Cell_Comps = GO_Cell_Comps[,-c(5,6,7)]
  GO_Mol_Funcs = GO_Mol_Funcs[,-c(5,6,7)]
  PATHWAYS = PATHWAYS[,-c(5,6,7)]
  PPIs = PPIs[,-c(5,6,7)]
  DISEASES = DISEASES[,-c(5,6,7)]
  
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
  return(c("GO BIOLOGICAL PROCESSES", GO_Processes[1:10,1], "GO CELL COMPONENTS", GO_Cell_Comps[1:10,1], "GO MOLECULAR FUNCTIONS", GO_Mol_Funcs[1:10,1], "PANTHER PATHWAYS", PATHWAYS[1:10,1], "PPIs", PPIs[1:10,1], "JENSEN DISEASE PHENOTYPES", DISEASES[1:10,1]))
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
DESeq2_AdultsItchDEGs = DESeq2_Adults3[c(1:579),]
    ## Filter rRNAs or any other unidentified mRNAs
DESeq2_AdultsItchDEGs = DESeq2_AdultsItchDEGs %>% filter(!grepl(DESeq2_AdultsItchDEGs$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

  ## Set up potential dataframes for downstream analyses (WT, Mut)
  ## Duplicate the entire transcriptional profile
AdultItchExp = RNASeqRepResultsAdultAll
  ## Prep the DEG subset
MutAdultItchExp = RNASeqRepResultsAdultAll[RNASeqRepResultsAdultAll$GeneID %in% c(DESeq2_AdultsItchDEGs$GeneID), ]
  ## Slice it 
MutAdultItchExp = MutAdultItchExp[,c(1,6:9)]
  ## Slice for WT
WTAdultItchExp = RNASeqRepResultsAdultAll[RNASeqRepResultsAdultAll$GeneID %in% c(DESeq2_AdultsItchDEGs$GeneID), ]
WTAdultItchExp = WTAdultItchExp[,c(1:5)]

##### Find Z-scores of the entire dataset and graph the profile centered around Tmem184b #####
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
Find_Z(Expression_Profile = AdultItchExp)


  ## Create a list of gene names ordered by euclidean distance by which to re-order the Z-scored dataframe
Euclid_dist_order_Z = hclust(dist(Z[,c(2:9)], method = "euclidean"))$order
  ## The names (not the numbers)
Euclid_dist_ord_Z_Genes = c(Z$GeneID[Euclid_dist_order_Z])


  ## Transform to order by clusters (Euclidean distances)
Itch_Corr_Z = Z %>%
  mutate(GeneID =  factor(GeneID, levels = Euclid_dist_ord_Z_Genes)) %>%
  arrange(GeneID)

  ## Find relevant reference genes to verify clustering and provide reference indeces for plotting
which(Euclid_dist_ord_Z_Genes == "Il31ra") # 1614
which(Euclid_dist_ord_Z_Genes == "Cysltr2") # 1613
which(Euclid_dist_ord_Z_Genes == "Npy2r") # 1947
which(Euclid_dist_ord_Z_Genes == "Sst") # 1623
which(Euclid_dist_ord_Z_Genes == "Htr1a") # 1616
which(Euclid_dist_ord_Z_Genes == "P2rx3") # 1307
which(Euclid_dist_ord_Z_Genes == "Lpar3") # 1871
which(Euclid_dist_ord_Z_Genes == "Lpar5") # 1689
which(Euclid_dist_ord_Z_Genes == "Scn11a") # 1883
which(Euclid_dist_ord_Z_Genes == "Scn10a") # 1893
which(Euclid_dist_ord_Z_Genes == "Mrgprd") # 1882
which(Euclid_dist_ord_Z_Genes == "Trpc6") # 2427
which(Euclid_dist_ord_Z_Genes == "Trpc3") # 1543
which(Euclid_dist_ord_Z_Genes == "F2rl2") # 1901
which(Euclid_dist_ord_Z_Genes == "Htr1f") # 1696
which(Euclid_dist_ord_Z_Genes == "Osmr") # 1596
which(Euclid_dist_ord_Z_Genes == "Fxyd2") # 2023
which(Euclid_dist_ord_Z_Genes == "Htr4") # 1785
which(Euclid_dist_ord_Z_Genes == "Mrgprx1") # 1702
which(Euclid_dist_ord_Z_Genes == "Ptgdr") # 2431
which(Euclid_dist_ord_Z_Genes == "Trpa1") # 1698
which(Euclid_dist_ord_Z_Genes == "Trpm6") # 1612
which(Euclid_dist_ord_Z_Genes == "Hrh1") # 1546
which(Euclid_dist_ord_Z_Genes == "Mrgpra3") # 2061
#which(Euclid_dist_ord_Z_Genes == "Nppb") # 134
which(Euclid_dist_ord_Z_Genes == "Tmem184b") # 1617
which(Euclid_dist_ord_Z_Genes == "Noxo1") # 4306


  ## Create a list of gene names ordered by euclidean distance by which to re-order the dataframe
Euclid_dist_order = hclust(dist(AdultItchExp[,c(2:9)], method = "euclidean"))$order
  ## The names (not the numbers)
Euclid_dist_ord_Genes = c(AdultItchExp$GeneID[Euclid_dist_order])

  ## Transform to order by clusters (Euclidean distances)
Itch_Corr = AdultItchExp %>%
  mutate(GeneID =  factor(GeneID, levels = Euclid_dist_ord_Genes)) %>%
  arrange(GeneID)

  ## Find relevant reference genes to verify clustering and provide reference indeces for plotting
which(Euclid_dist_ord_Genes == "Il31ra") # 1980
which(Euclid_dist_ord_Genes == "Cysltr2") # 5143
which(Euclid_dist_ord_Genes == "Npy2r") # 1611
which(Euclid_dist_ord_Genes == "Sst") # 5137
which(Euclid_dist_ord_Genes == "Htr1a") # 5138
which(Euclid_dist_ord_Genes == "P2rx3") # 1612
which(Euclid_dist_ord_Genes == "Lpar3") # 9989
which(Euclid_dist_ord_Genes == "Lpar5") # 4839
which(Euclid_dist_ord_Genes == "Scn11a") # 11226
which(Euclid_dist_ord_Genes == "Scn10a") # 54
which(Euclid_dist_ord_Genes == "Mrgprd") # 755
which(Euclid_dist_ord_Genes == "Trpc6") # 4841
which(Euclid_dist_ord_Genes == "Trpc3") # 10805
which(Euclid_dist_ord_Genes == "F2rl2") # 10857
which(Euclid_dist_ord_Genes == "Htr1f") # 5140
which(Euclid_dist_ord_Genes == "Osmr") # 5134
which(Euclid_dist_ord_Genes == "Fxyd2") # 72
which(Euclid_dist_ord_Genes == "Htr4") # 4840
which(Euclid_dist_ord_Genes == "Mrgprx1") # 4741
which(Euclid_dist_ord_Genes == "Ptgdr") # 1929
which(Euclid_dist_ord_Genes == "Trpa1") # 10394
which(Euclid_dist_ord_Genes == "Trpm6") # 6650
which(Euclid_dist_ord_Genes == "Hrh1") # 6283
which(Euclid_dist_ord_Genes == "Mrgpra3") # 5064
#which(Euclid_dist_ord_Genes == "Nppb") # 
which(Euclid_dist_ord_Genes == "Tmem184b") # 9884
which(Euclid_dist_ord_Genes == "Noxo1") # 6717


## Confirm by indexing; copy and paste to customize the clustered heatmap
Euclid_dist_ord_Z_Genes[c(1614,1613,1947,1623,1616,1307,1871,1689,1883,1893,1882,2427,1543,1901,1696,1596,2023,1785,1702,2431,1698,1612,1546,2061,1617,4306)]
## Confirm the location/order of the transformed matrix/dataframe
## Should be "Il31ra"
Itch_Corr_Z[c(1614,1613,1947,1623,1616,1307,1871,1689,1883,1893,1882,2427,1543,1901,1696,1596,2023,1785,1702,2431,1698,1612,1546,2061,1617,4306),]



  ## Confirm by indexing; copy and paste to customize the clustered heatmap
Euclid_dist_Genes[c(1980,5143,1611,5137,5138,1612,9989,4839,11226,54,755,4841,10805,10857,5140,5134,72,4840,4741,1929,10394,6650,6283,5064,9884,6717)]
  ## Confirm the location/order of the transformed matrix/dataframe
  ## Should be "Il31ra"
Itch_Corr[c(1980,5143,1611,5137,5138,1612,9989,4839,11226,54,755,4841,10805,10857,5140,5134,72,4840,4741,1929,10394,6650,6283,5064,9884,6717),]


  ## Visualize the Z-scored, Euclidean-clustered and ordered TPMs with "pheatmap"
pheatmap(mat = Itch_Corr_Z[,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, show_rownames = F)
  
#pheatmap(mat = Itch_Corr_Z[1607:1627,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, labels_row = Euclid_dist_ord_Z_Genes[1607:1627])

  ## Visualize the Z-scored, Euclidean-clustered and ordered TPMs around the TMEM cluster with "pheatmap"
pheatmap(mat = Itch_Corr_Z[1592:1642,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, labels_row = Euclid_dist_ord_Z_Genes[1592:1642])

  ## Visualize the non-Z-scored, Euclidean-clustered and ordered TPMs
pheatmap(mat = Itch_Corr[,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, show_rownames = F)
  ## --> Obscured by genes at the bottom of the dataframe that comprise the maximal transcriptional values

  ## Visualize the Tmem184b cluster of un-Z-scored data
pheatmap(mat = Itch_Corr[9884:9934,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, labels_row = Euclid_dist_ord_Genes[9884:9934])

  ## Peek at Tmem to see how the pheatmap values compare to the dataframe values
Itch_Corr_Z[1617,]



  ## Perform gene ontology and pathway analysis of the Tmem184b cluster (ballparked)
Functional_Analysis(Genes = c(Euclid_dist_ord_Z_Genes[1592:1642]))
#Functional_Analysis(Genes = c(Euclid_dist_ord_Z_Genes[1607:1627]))
Functional_Analysis(Genes = c(Euclid_dist_ord_Genes[9884:9934]))



##### Hack the WGCNA Tutorial. Prep data for net creation across WTs with default and optimized parameters; repeat for Mutants #####
  ## Hack the WGCNA Tutorial. Determine the soft thresholding power for "inclusivity" within eigengene modules
  ## Create a data structure to house the data for downstream manipulation and analysis
multiExpr = vector(mode = "list", length = 1)
  # Create a vector as long as the expression dataframe for determining additional genes to remove upstream of correlation calculations
expcheck = as.array(vector(length = nrow(AdultItchExp)))

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
    ## Run the function
Remove_NAN_Generators(Expression_Profile = AdultItchExp, Samples_to_Compare = "WT")


    ## Fill the data structure with the transposed Adult DRG WT Expression profile
multiExpr[[1]] = list(data = as.data.frame(t(Expression_Profile_cleaned[,2:5])))
    ## Rename the columns by their GeneIDs
names(multiExpr[[1]]$data) = Expression_Profile_cleaned$GeneID
    ## Rename the rows by their Sample IDs
rownames(multiExpr[[1]]$data) = names(Expression_Profile_cleaned[2:5])
  
  ## Take a range of potential powers to determine the appropriate power by which the correlation matrix will be raised. This can be thought of as increasing the "connectivity" of more related genes
  
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
  ### -->> Model isn't right due to lack of data.. when all samples are used, softThreshold of ~6 is used. So, use power of 6 anyway ###
  
  
  
  
  
  ##### WT nets #####
    ## "Automatically" create a network using default settings of the "blockwiseModules" WGCNA function (for datasets > 5000 genes)
  WTnet = blockwiseModules(datExpr = multiExpr[[1]]$data,
                         checkMissingData = T,
                         maxBlockSize = 11300,
                         blockSizePenaltyPower = Inf,
                         corType = "pearson",
                         power = 6,
                         networkType = "signed",
                         TOMType = "signed",
                         pamRespectsDendro = F,
                         numericLabels = F,
                         verbose = 4)
  
  
    ## Generate a variable to store the hierarchically clustered dendrogram
  WTGeneDendrogram = WTnet$dendrograms[[1]]
    ## Generate a variable to store the module labels (colors in this case)
  WTModuleEigengenes = c(names(WTnet$MEs))
    ## Generate a variable that stores all genes' associated module colors
  WTModuleLabels = WTnet$colors
    ## Generate a variable that stores the conversion of those module labels to colors
  WTModuleColors = labels2colors(WTModuleLabels)
    ## Find the Tmem module
  WTnet$colors[9726]
    ## Find the Tmem module before merging of similar modules
  WTnet$unmergedColors[9726]
    ## Find all the genes with the same module label
  WTTmem_module = c(which(WTnet$colors == "turquoise"))
    ## Find the names of all those genes
  WTTmem_module = c(names(WTTmem_module))
    ## Print out the module
  WTTmem_module
    ## Run downstream bioinformatics analysis through enrichR on that module
  Functional_Analysis(Genes = WTTmem_module)
  
  ## Find if there were any genes that confounded the clustering
  #which(WTnet$goodGenes == "FALSE")
  ## Remove those genes and re-determine modules/colors if necessary
  #WTModuleLabels = WTnet$colors[-7816]
  ## Re-calculate the module colors if necessary
  #WTModuleColors = labels2colors(WTModuleLabels)
  
  ## Plot the dendrogram and associated Module colors using the variables previously stored
  plotDendroAndColors(dendro = WTGeneDendrogram,
                      colors = WTModuleColors,
                      "Module Colors",
                      dendroLabels = F,
                      hang = 0.01,
                      addGuide = T,
                      guideHang = 0.01,
                      main = "WT DRG Gene Default Dendrogram")
  
  which(WTnet$goodGenes == "FALSE")
  which(Expression_Profile_cleaned[,2] == 0)
  which(Expression_Profile_cleaned[,3] == 0)
  which(Expression_Profile_cleaned[,4] == 0)
  which(Expression_Profile_cleaned[,5] == 0)
  Expression_Profile_cleaned[10474,]
  
    ## Compare with more strict clustering parameters
  WTOptinet = blockwiseModules(datExpr = multiExpr[[1]]$data,
                             checkMissingData = T,
                             maxBlockSize = 11300,
                             blockSizePenaltyPower = Inf,
                             corType = "pearson",
                             power = 6,
                             networkType = "signed",
                             TOMType = "signed",
                             deepSplit = 2,
                             detectCutHeight = 0.995,
                             mergeCutHeight = 0.10,
                             minCoreKME = 0.6,
                             minCoreKMESize = 15,
                             minKMEtoStay = 0.5,
                             pamRespectsDendro = F,
                             numericLabels = F,
                             verbose = 4)
  
    ## Generate a variable to store the hierarchically clustered dendrogram
  WTOptiGeneDendrogram = WTOptinet$dendrograms[[1]]
    ## Generate a variable to store the module labels (colors in this case)
  WTOptiModuleEigengenes = c(names(WTOptinet$MEs))
    ## Generate a variable that stores all genes' associated module colors
  WTOptiModuleLabels = WTOptinet$colors
    ## Generate a variable that stores the conversion of those module labels to colors
  WTOptiModuleColors = labels2colors(WTOptiModuleLabels)
    ## Find the Tmem module
  WTOptinet$colors[9726]
    ## Find the Tmem module before merging of similar modules
  WTOptinet$unmergedColors[9726]
    ## Find all the genes with the same module label
  WTOptiTmem_module = c(which(WTOptinet$colors == "royalblue"))
    ## Find the names of all those genes
  WTOptiTmem_module = c(names(WTOptiTmem_module))
    ## Print out the module
  WTOptiTmem_module
    ## Run downstream bioinformatics analysis through enrichR on that module
  Functional_Analysis(Genes = WTOptiTmem_module)
  
  
    ## Find if there were any genes that confounded the clustering
  #which(WTOptinet$goodGenes == "FALSE")
    ## Remove those genes and re-determine modules/colors if necessary
  #WTOptiModuleLabels = WTOptinet$colors[-7816]
    ## Re-calculate the module colors if necessary
  #WTOptiModuleColors = labels2colors(WTOptiModuleLabels)
  
    ## Plot the dendrogram and associated Module colors using the variables previously stored
  plotDendroAndColors(dendro = WTOptiGeneDendrogram,
                      colors = WTOptiModuleColors,
                      "Module Colors",
                      dendroLabels = F,
                      hang = 0.01,
                      addGuide = T,
                      guideHang = 0.01,
                      main = "WT DRG Gene (Optimized?) Dendrogram")
  
 
  
  
##### Repeat for Mutants #####
  multiExpr = vector(mode = "list", length = 1)
    ## Create a vector as long as the expression dataframe for determining additional genes to remove upstream of correlation calculations
  expcheck = as.array(vector(length = nrow(AdultItchExp)))
    ## Run the function to remove genes that have too low variation or no expression
  Remove_NAN_Generators(Expression_Profile = AdultItchExp, Samples_to_Compare = "Mut")
  
    ## Fill the data structure with the transposed Adult DRG WT Expression profile
  multiExpr[[1]] = list(data = as.data.frame(t(Expression_Profile_cleaned[,2:5])))
    ## Rename the columns by their GeneIDs
  names(multiExpr[[1]]$data) = Expression_Profile_cleaned$GeneID
    ## Rename the rows by their Sample IDs
  rownames(multiExpr[[1]]$data) = names(Expression_Profile_cleaned[2:5])
  
    ## Take a range of potential powers to determine the appropriate power by which the correlation matrix will be raised. This can be   thought of as increasing the intramodular "connectivity" of more related genes
  
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
  plot(Topology_df$Power, (-Topology_df$slope)*Topology_df$SFT.R.sq,
       xlab = "Soft Threshold (power)",
       ylab = colNames[1], type = "n", 
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(0,1), main = colNames[1])
  addGrid()
  text(Topology_df$Power, (-Topology_df$slope)*Topology_df$SFT.R.sq, labels = powers)
  
    ## Graph the mean connectivity as a function of soft threshold (power)
  plot(Topology_df$mean.k., (-Topology_df$slope)*Topology_df$mean.k.,
       xlab = "Soft Threshold (power)",
       ylab = colNames[2], type = "n",
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(min(-Topology_df$slope*Topology_df$mean.k.), max(-Topology_df$slope*Topology_df$mean.k.)), main = colNames[2])
  addGrid()
  text(Topology_df$Power, (-Topology_df$slope)*Topology_df$mean.k., labels = powers)
  
    ## Graph the median connectivity as a function of soft threshold (power)
  plot(Topology_df$median.k., (-Topology_df$slope)*Topology_df$median.k.,
       xlab = "Soft Threshold (power)", ylab = colNames[3], type = "n",
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(min(-Topology_df$slope*Topology_df$median.k.), max(-Topology_df$slope*Topology_df$median.k.)), main = colNames[3])
  addGrid()
  text(Topology_df$Power, (-Topology_df$slope)*Topology_df$median.k., labels = powers)
  
    ## Graph the maximum connectivity as a function of soft threshold (power)
  plot(Topology_df$max.k., (-Topology_df$slope)*Topology_df$max.k.,
       xlab = "Soft Threshold (power)", ylab = colNames[4], type = "n",
       xlim = c(min(Topology_df$Power), max(Topology_df$Power)),
       ylim = c(min(-Topology_df$slope*Topology_df$max.k.), max(-Topology_df$slope*Topology_df$max.k.)), main = colNames[4])
  addGrid()
  text(Topology_df$Power, (-Topology_df$slope)*Topology_df$max.k., labels = powers)
  ### -->> Model doesn't seem right.. use power of 6 anyway ###
  
  
  ##### Mutant nets #####
    ## "Automatically" create a network using default settings of the "blockwiseModules" WGCNA function (for datasets > 5000 genes)
  Mutnet = blockwiseModules(datExpr = multiExpr[[1]]$data,
                         checkMissingData = T,
                         maxBlockSize = 11300,
                         blockSizePenaltyPower = Inf,
                         corType = "pearson",
                         power = 6,
                         networkType = "signed",
                         TOMType = "signed",
                         pamRespectsDendro = F,
                         numericLabels = F,
                         verbose = 4)
  
    ## Generate a variable to store the hierarchically clustered dendrogram
  MutGeneDendrogram = Mutnet$dendrograms[[1]]
    ## Generate a variable to store the module labels (colors in this case)
  MutModuleEigengenes = c(names(Mutnet$MEs))
    ## Generate a variable that stores all genes' associated module colors
  MutModuleLabels = Mutnet$colors
    ## Generate a variable that stores the conversion of those module labels to colors
  MutModuleColors = labels2colors(MutModuleLabels)
    ## Find the Tmem module
  Mutnet$colors[9725]
    ## Find the Tmem module before merging of similar modules
  Mutnet$unmergedColors[9725]
    ## Find all the genes with the same module label
  MutTmem_module = c(which(Mutnet$colors == "black"))
    ## Find the names of all those genes
  MutTmem_module = c(names(MutTmem_module))
    ## Print out the module
  MutTmem_module
    ## Run downstream bioinformatics analysis through enrichR on that module
  Functional_Analysis(Genes = MutTmem_module)
  
    ## Find if there were any genes that confounded the clustering
  which(Mutnet$goodGenes == "FALSE")
  ## Remove those genes and re-determine modules/colors if necessary
  #MutModuleLabels = Mutnet$colors[-7816]
  ## Re-calculate the module colors if necessary
  #MutModuleColors = labels2colors(MutModuleLabels)
  
    ## Plot the dendrogram and associated Module colors using the variables previously stored
  plotDendroAndColors(dendro = MutGeneDendrogram,
                      colors = MutModuleColors,
                      "Module Colors",
                      dendroLabels = F,
                      hang = 0.01,
                      addGuide = T,
                      guideHang = 0.01,
                      main = "Mut DRG Gene Default Dendrogram")
  
  
  
    ## Compare with more strict clustering parameters
  MutOptinet = blockwiseModules(datExpr = multiExpr[[1]]$data,
                             checkMissingData = T,
                             maxBlockSize = 11300,
                             blockSizePenaltyPower = Inf,
                             corType = "pearson",
                             power = 6,
                             networkType = "signed",
                             TOMType = "signed",
                             deepSplit = 2,
                             detectCutHeight = 0.995,
                             mergeCutHeight = 0.10,
                             minCoreKME = 0.6,
                             minCoreKMESize = 15,
                             minKMEtoStay = 0.5,
                             pamRespectsDendro = F,
                             numericLabels = F,
                             verbose = 4)
  
  
    ## Generate a variable to store the hierarchically clustered dendrogram
  MutOptiGeneDendrogram = MutOptinet$dendrograms[[1]]
    ## Generate a variable to store the module labels (colors in this case)
  MutOptiModuleEigengenes = c(names(MutOptinet$MEs))
    ## Generate a variable that stores all genes' associated module colors
  MutOptiModuleLabels = MutOptinet$colors
    ## Generate a variable that stores the conversion of those module labels to colors
  MutOptiModuleColors = labels2colors(MutOptiModuleLabels)
    ## Find the Tmem module
  MutOptinet$colors[9725]
    ## Find the Tmem module before merging of similar modules
  MutOptinet$unmergedColors[9725]
    ## Find all the genes with the same module label
  MutOptiTmem_module = c(which(MutOptinet$colors == "darkolivegreen"))
    ## Find the names of all those genes
  MutOptiTmem_module = c(names(MutOptiTmem_module))
    ## Print out the module
  MutOptiTmem_module
    ## Run downstream bioinformatics analysis through enrichR on that module
  Functional_Analysis(Genes = MutOptiTmem_module)
  
    ## Find if there were any genes that confounded the clustering
  which(MutOptinet$goodGenes == "FALSE")
  ## Remove those genes and re-determine modules/colors if necessary
  #MutOptiModuleLabels = MutOptinet$colors[-7816]
  ## Re-calculate the module colors if necessary
  #MutOptiModuleColors = labels2colors(MutOptiModuleLabels)
  
    ## Plot the dendrogram and associated Module colors using the variables previously stored
  plotDendroAndColors(dendro = MutOptiGeneDendrogram,
                      colors = MutOptiModuleColors,
                      "Module Colors",
                      dendroLabels = F,
                      hang = 0.01,
                      addGuide = T,
                      guideHang = 0.01,
                      main = "Mut DRG Gene (Optimized?) Dendrogram")
  

###### Build a hierarchically clustered dendrogram for all samples (the entire network, technically as one genotype) using the WGCNA package "blockwiseModules ######

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
  
net = blockwiseModules(datExpr = multiExpr[[1]]$data,
                                checkMissingData = T,
                                maxBlockSize = 11300,
                                blockSizePenaltyPower = Inf,
                                corType = "pearson",
                                power = 6,
                                networkType = "signed",
                                TOMtype = "signed",
                                deepSplit = 2,
                                detectCutHeight = 0.995,
                                pamRespectsDendro = F,
                                mergeCutHeight = 0.1,
                                numericLabels = F,
                                minCoreKME = 0.6,
                                minCoreKMESize = 15,
                                minKMEtoStay = 0.5,
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
plotDendroAndColors(dendro = GeneDendrogram,
                    colors = ModuleColors,
                    "Module Colors",
                    dendroLabels = F,
                    hang = 0.01,
                    addGuide = T,
                    guideHang = 0.01,
                    main = "Tmem184b GT DRG Gene Dendrogram",
                    marAll = c(4,4,6,4))



  ##### Build a hierarchically clustered dendrogram for all samples using the "blockwiseConsensusModules" WGCNA function #####

multiExpr = vector(mode = "list", length = 2)
  ## Create a vector as long as the expression dataframe for determining additional genes to remove upstream of correlation calculations
expcheck = as.array(vector(length = nrow(AdultItchExp)))

  ## Run the NAN function
Remove_NAN_Generators(Expression_Profile = AdultItchExp, Samples_to_Compare = "All")


  ## Fill the data structure with the transposed Adult DRG WT Expression profile
multiExpr[[1]] = list(data = as.data.frame(t(Expression_Profile_cleaned[,2:5])))
multiExpr[[2]] = list(data = as.data.frame(t(Expression_Profile_cleaned[,6:9])))

  ## Rename the columns by their GeneIDs
names(multiExpr[[1]]$data) = Expression_Profile_cleaned$GeneID
names(multiExpr[[2]]$data) = Expression_Profile_cleaned$GeneID
  ## Rename the rows by their Sample IDs
rownames(multiExpr[[1]]$data) = names(Expression_Profile_cleaned[2:5])
rownames(multiExpr[[2]]$data) = names(Expression_Profile_cleaned[6:9])

  ## How many genotypes are there?
nSets = 2

  ## Take a range of potential powers to determine the appropriate power by which the correlation matrix will be raised. This can be thought of as enhancing boundaries of relatedness of genes.

powers = c(seq(4,10,by=1), seq(12,20, by=2))
  ## Use the WGCNA function, "pickSoftThreshold", to compute model fits of various powers, and to compute the "connectivity" metrics of the data based on these power values
powerTables = vector(mode = "list", length = 2)
for (set in 1:nSets){
  powerTables[[set]] = list(data = pickSoftThreshold(data = multiExpr[[set]]$data, powerVector = powers, verbose = 2)[[2]])
}

collectGarbage()
colors = c("black", "blue")
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")

  ## Re-arrange the data output of the soft thresholding function into a dataframe for graphics and analysis
Topology_df = data.frame(powerTables[[1]]$data[1],
                         powerTables[[1]]$data[2], 
                         powerTables[[1]]$data[3], 
                         powerTables[[1]]$data[4], 
                         powerTables[[1]]$data[5], 
                         powerTables[[1]]$data[6], 
                         powerTables[[1]]$data[7], 
                         powerTables[[2]]$data[1],
                         powerTables[[2]]$data[2], 
                         powerTables[[2]]$data[3], 
                         powerTables[[2]]$data[4], 
                         powerTables[[2]]$data[5], 
                         powerTables[[2]]$data[6], 
                         powerTables[[2]]$data[7],stringsAsFactors = F)
powerTables[[1]]$powerEstimate
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
  ### -->> Model doesn't seem right.. power of 6 ###

    ## Repeat for mutants
plot(Topology_df$Power.1, Topology_df$slope.1*Topology_df$SFT.R.sq.1,
     xlab = "Soft Threshold (power)",
     ylab = colNames[1], type = "n", 
     xlim = c(min(Topology_df$Power.1), max(Topology_df$Power.1)),
     ylim = c(0,1), main = colNames[1])
addGrid()
text(Topology_df$Power.1, Topology_df$slope.1*Topology_df$SFT.R.sq.1, labels = powers)

  ## Graph the mean connectivity as a function of soft threshold (power)
plot(Topology_df$mean.k..1, Topology_df$slope.1*Topology_df$mean.k..1,
     xlab = "Soft Threshold (power)",
     ylab = colNames[2], type = "n",
     xlim = c(min(Topology_df$Power.1), max(Topology_df$Power.1)),
     ylim = c(min(Topology_df$slope*Topology_df$mean.k..1), max(Topology_df$slope*Topology_df$mean.k..1)), main = colNames[2])
addGrid()
text(Topology_df$Power.1, (Topology_df$slope.1)*Topology_df$mean.k..1, labels = powers)

  ## Graph the median connectivity as a function of soft threshold (power)
plot(Topology_df$median.k..1, Topology_df$slope.1*Topology_df$median.k..1,
     xlab = "Soft Threshold (power)", ylab = colNames[3], type = "n",
     xlim = c(min(Topology_df$Power.1), max(Topology_df$Power.1)),
     ylim = c(min(Topology_df$slope.1*Topology_df$median.k..1), max(Topology_df$slope.1*Topology_df$median.k..1)), main = colNames[3])
addGrid()
text(Topology_df$Power.1, (Topology_df$slope.1)*Topology_df$median.k..1, labels = powers)

  ## Graph the maximum connectivity as a function of soft threshold (power)
plot(Topology_df$max.k..1, -Topology_df$slope.1*Topology_df$max.k..1,
     xlab = "Soft Threshold (power)", ylab = colNames[4], type = "n",
     xlim = c(min(Topology_df$Power.1), max(Topology_df$Power.1)),
     ylim = c(min(Topology_df$slope.1*Topology_df$max.k..1), max(Topology_df$slope.1*Topology_df$max.k..1)), main = colNames[4])
addGrid()
text(Topology_df$Power.1, (Topology_df$slope.1)*Topology_df$max.k..1, labels = powers)


net = blockwiseConsensusModules(multiExpr = multiExpr,
                                checkMissingData = T,
                                maxBlockSize = 12000,
                                blockSizePenaltyPower = Inf,
                                corType = "pearson",
                                power = 6,
                                networkType = "signed",
                                TOMType = "signed",
                                deepSplit = 2,
                                detectCutHeight = 0.995,
                                pamRespectsDendro = F,
                                mergeCutHeight = 0.10,
                                numericLabels = F,
                                minCoreKME = 0.6,
                                minCoreKMESize = 15,
                                minKMEtoStay = 0.5,
                                verbose = 4)


WTModuleEigengenes1 = c(names(net$multiMEs[[1]]$data))
MutModuleEigengenes1 = c(names(net$multiMEs[[2]]$data))
## Create a variable for the merged eigenvectors for plotting
ModuleLabels1 = net$colors
## Create a variable of the dendrogram for plotting
GeneDendrogram1 = net$dendrograms[[1]]
## Create a variable (vector) of colors corresponding to the eigenvectors
ModuleColors1 = labels2colors(ModuleLabels)

## Find Tmem's index to find it in the net
which(RNASeqRepResultsAdultAll$GeneID == "Tmem184b")
## Use that index to find Tmem in the network
net$colors[9727]
## Find Tmem in the pre-merged network
net$unmergedColors[9727]
## Find the Tmem module
Tmem_module1 = c(which(net$colors == "slateblue1"))
## Remove the indeces/numbers for downstream analysis
Tmem_module1 = c(names(Tmem_module1))
## Check out the module
Tmem_module1
## Perform downstream bioinformatics (GO, pathway, PPI analysis)
Functional_Analysis(Genes = Tmem_module1)


## Plot the dendrogram
plotDendroAndColors(dendro = GeneDendrogram1,
                    colors = ModuleColors1,
                    "Module Colors",
                    dendroLabels = F,
                    hang = 0.001,
                    addGuide = T,
                    guideHang = 0.001,
                    main = "Consensus Gene Dendrogram and Module Colors",
                    marAll = c(4,4,6,4))

## Find which genes are "bad" (have "too many missing samples or zero variance")
which(net$goodGenes == "FALSE")
AdultItchExp[c(4956,7816,8485),]

## Omit them from the colors to color the dendrogram
ModuleLabels1 = net$colors[-c(4956,7816,8485)]
ModuleColors1 = labels2colors(ModuleLabels1)

## Re-plot the dendrogram
plotDendroAndColors(dendro = GeneDendrogram1,
                    colors = ModuleColors1,
                    "Module Colors",
                    dendroLabels = F,
                    hang = 0.001,
                    addGuide = T,
                    guideHang = 0.001,
                    main = "Consensus Tmem184b GT DRG Dendrogram",
                    marAll = c(4,4,6,4))

  ##### Clustering validation ... In progress #####
  ## Evaluate the hierarchical clustering
a = clValid(obj = as.data.frame(t(multiExpr[[1]]$data)), nClust = 25, clMethods = "hierarchical", validation = "internal", method = "average", maxitems = nrow(t(multiExpr[[1]]$data)))

  ## Dunn's index should be high, although given the noise in the dataset, low values are not surprising; could be overfit, or relationships might be too "blurry". The network is weighted and signed..
a@measures[]

#dunn(distance = NULL, clusters = net$colors)

  ## Create a network from the ground, up
  ## Create the adjacency matrix, which houses the co-relational (similarity/distance) measures of genes between genes.
    ## Raise to the power that maximizes the scale free topology model fit and the mean or median connectivity (6 in this case)
adjacencyM = abs(cor(multiExpr[[1]]$data, use = "p"))^6
  ## Make the adjacency a dissimilarity matrix
dissADJ = 1 - adjacencyM
  ## Convert the adjacency to a distance/DISsimilarity matrix
collectGarbage()
dissTOM = TOMdist(adjacencyM)
  ## Cluster the data for a dendrogram
hierdissTOM = hclust(as.dist(dissTOM), method = "average")
hierADJ = hclust(as.dist(adjacencyM), method = "average", d = dissADJ)
hierDissADJ = hclust(as.dist(dissADJ), method = "average")

  ## Differe TOMsimilarity function
TOM = TOMsimilarityFromExpr(datExpr = multiExpr, 
                            corType = "pearson",
                            networkType = "signed",
                            power = 6,
                            TOMType = "signed")


  ## Create a grouping using the TOM dissimilarity distance matrix to derive hierarchically clustered genes; dendrogram is cut "dynamically"
branchnumber = cutreeDynamic(hierdissTOM, method = "tree")
colorDynamicTOM = labels2colors(branchnumber)

  ## Create an an alternate grouping "statically" cut (accurate at the nodes, not sensitive to weaker on fringes)
colorStaticTOM = as.character(cutreeStaticColor(hierdissTOM, cutHeight = 0.7, minSize = 20))
#colorDynamicTOM = labels2colors(cutreeDynamic(hierdissTOM), method = "tree")

  ## Create a third module grouping based on a hybrid approach (cross between static and hybrid)
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierdissTOM, distM = dissTOM, cutHeight = 0.7, deepSplit = 2, pamRespectsDendro = F))

  ## Plot the groupings
plotDendroAndColors(dendro = hierdissTOM,
                    colors = data.frame(colorStaticTOM, colorDynamicTOM, colorDynamicHybridTOM),
                    dendroLabels = FALSE,
                    marAll = c(1, 8, 3, 1),
                    main = "WT DRG Gene Dendrogram and Module Colors")

  ## Find intramodular connectivity
Alldegrees1 = intramodularConnectivity(dissTOM, colorDynamicHybridTOM)

multiExpr[[1]]$data[9726]

names(c(which(WTnet$colors == "violet")))
names(c(which(WTnet$unmergedColors == "skyblue3")))

intModules = c("violet")
for (module in intModules){
  moduleGenes = names(c(which(WTnet$colors == "violet")))
}

Mod_Enrichment_Analysis = GOenrichmentAnalysis(labels = WTModuleColors,entrezCodes = names(WTModuleLabels),organism = "mouse",ontologies = c("BP", "CC", "MF"),nBestP = 10,nBiggest = 10)


