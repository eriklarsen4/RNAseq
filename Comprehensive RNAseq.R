


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
diff_DBs = DBs$libraryName[c(50,127,91,92,109)]
DBs = DBs$libraryName[c(130,131,132,102,14,110)]



##### Data Prep; import DGEA adult DRG file  #####
  ## Import the dataset you want to analyze (adult Tmem184b-GT/GT DRG DGE)
aDRG = read.csv("M:/Erik/Data/Omics/RNAseq/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
aDRG3 = subset(aDRG, (!is.na(aDRG[,"AdjP"])))

  ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
aDRG9 = aDRG3 %>% filter(!grepl(aDRG3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

#mt.+.?$|

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
aDRG9$labs = aDRG9$GeneID
  ## Replace DEGs with the string, "DEGs"
aDRG9$labs[which(aDRG9$AdjP <= 0.05)]= "DEGs"
  ## Replace the remaining genes with "Non-DEGs"
aDRG9$labs[which(aDRG9$AdjP >= 0.05)]= "Non-DEGs"

aDRG_DEG_list = c(aDRG9$GeneID[which(aDRG9$labs == "DEGs")])



##### Data Prep; import DGEA e13 DRG file  #####
  ## Import the dataset you want to analyze (e13 Tmem184b-GT/GT DRG DGE)
e13_DRG = read.csv("M:/Erik/Data/Omics/TimeCourse/Processed Galaxy Output/Test Results to Upload/DESeq2_result_file_on_GT_e13_PARA.csv")
## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
e13_DRG3 = subset(e13_DRG, (!is.na(e13_DRG[,"AdjP"])))

## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
e13_DRG9 = e13_DRG3 %>% filter(!grepl(e13_DRG3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

#mt.+.?$|

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
e13_DRG9$labs = e13_DRG9$GeneID
    ## Replace DEGs with the string, "DEGs"
e13_DRG9$labs[which(e13_DRG9$AdjP <= 0.01)]= "DEGs"
    ## Replace the remaining genes with "Non-DEGs"
e13_DRG9$labs[which(e13_DRG9$AdjP >= 0.01)]= "Non-DEGs"

e13_DRG_DEG_list = c(e13_DRG9$GeneID[which(e13_DRG9$labs == "DEGs")])


##### Data Prep; import DGEA P0 DRG file  #####
  ## Import the dataset you want to analyze (P0 Tmem184b-GT/GT DRG DGE)
P0_DRG = read.csv("M:/Erik/Data/Omics/TimeCourse/Processed Galaxy Output/Test Results to Upload/DESeq2_result_file_on_GT_p0_PARA.csv")
## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
P0_DRG3 = subset(P0_DRG, (!is.na(P0_DRG[,"AdjP"])))

## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
P0_DRG9 = P0_DRG3 %>% filter(!grepl(P0_DRG3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

#mt.+.?$|

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
P0_DRG9$labs = P0_DRG9$GeneID
    ## Replace DEGs with the string, "DEGs"
P0_DRG9$labs[which(P0_DRG9$AdjP <= 0.05)]= "DEGs"
    ## Replace the remaining genes with "Non-DEGs"
P0_DRG9$labs[which(P0_DRG9$AdjP >= 0.05)]= "Non-DEGs"

P0_DRG_DEG_list = c(P0_DRG9$GeneID[which(P0_DRG9$labs == "DEGs")])


##### Data Prep; import DGEA P10 DRG file  #####
  ## Import the dataset you want to analyze (P10 Tmem184b-GT/GT DRG DGE)
P10_DRG = read.csv("M:/Erik/Data/Omics/TimeCourse/Processed Galaxy Output/Test Results to Upload/DESeq2_result_file_on_GT_p10_PARA.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
P10_DRG3 = subset(P10_DRG, (!is.na(P10_DRG[,"AdjP"])))

  ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
P10_DRG9 = P10_DRG3 %>% filter(!grepl(P10_DRG3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

#mt.+.?$|

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
P10_DRG9$labs = P10_DRG9$GeneID
    ## Replace DEGs with the string, "DEGs"
P10_DRG9$labs[which(P10_DRG9$AdjP <= 0.05)]= "DEGs"
    ## Replace the remaining genes with "Non-DEGs"
P10_DRG9$labs[which(P10_DRG9$AdjP >= 0.05)]= "Non-DEGs"

P10_DRG_DEG_list = c(P10_DRG9$GeneID[which(P10_DRG9$labs == "DEGs")])


##### Data Prep; import DGEA Young Mutant v WT Hippocampus file  #####
  ## Import the dataset you want to analyze (Young Tmem184b-GT/GT Hpc DGE)
YH = read.csv("M:/Erik/Data/Omics/Hippocampus/Genotype Comparisons/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results for Young Hippocampus without ymu 4141.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
YH3 = subset(YH, (!is.na(YH[,"AdjP"])))

  ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
YH9 = YH3 %>% filter(!grepl(YH3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

#mt.+.?$|

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
YH9$labs = YH9$GeneID
    ## Replace DEGs with the string, "DEGs"
YH9$labs[which(YH9$AdjP <= 0.05)]= "DEGs"
    ## Replace the remaining genes with "Non-DEGs"
YH9$labs[which(YH9$AdjP >= 0.05)]= "Non-DEGs"

YH_DEG_list = c(YH9$GeneID[which(YH9$labs == "DEGs")])


##### Data Prep; import DGEA Old Mutant v WT Hippocampus file  #####
  ## Import the dataset you want to analyze (Old Tmem184b-GT/GT Hpc DGE)
OH = read.csv("M:/Erik/Data/Omics/Hippocampus/Genotype Comparisons/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results for Old Genotype without 1849 and 1880.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
OH3 = subset(OH, (!is.na(OH[,"AdjP"])))

  ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
OH9 = OH3 %>% filter(!grepl(OH3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

#mt.+.?$|

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
OH9$labs = OH9$GeneID
    ## Replace DEGs with the string, "DEGs"
OH9$labs[which(OH9$AdjP <= 0.05)]= "DEGs"
    ## Replace the remaining genes with "Non-DEGs"
OH9$labs[which(OH9$AdjP >= 0.05)]= "Non-DEGs"

OH_DEG_list = c(OH9$GeneID[which(OH9$labs == "DEGs")])

##### Downstream analyses #####
  ## Check which genes are DE in multiple datsasets (iterate through)
aDRG9$GeneID[which(aDRG_DEG_list %in% e13_DRG_DEG_list)]
#P0_DRG_DEG_list
#P10_DRG_DEG_list
#YH_DEG_list
#OH_DEG_list

##### Individual GO Analysis #####
  ## Find Entrez IDs (Ensembl IDs) from a DEG list; remember to re-define as a key when mapping
Ensembl_DEG_IDs = as.integer(
  mapIds(org.Mm.eg.db,
         keys = as.character(
           c(DEG_list[])
         ),
         keytype = "SYMBOL",
         column = "ENTREZID")
)
  ## Extract GO Terms for all genes in the DEG list; remember to re-define as a key when mapping
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
DEG_indeces_without_Ensembl_IDs = c(
  which(is.na(Ensembl_DEG_IDs[]) == TRUE)
)
  ## Repeat for GO IDs
DEG_indeces_without_GO_IDs = as.integer(c(
  which(is.na(DEG_GO_IDs[]) == TRUE)
))
  ## There should be more genes that **don't** have GO terms than genes that are un-annotated ("without GO ID" should be longer than "without Ensembl ID"; check)
aDRG9$GeneID[c(DEG_indeces_without_Ensembl_IDs, DEG_indeces_without_GO_IDs)]
  ## Concatenate the genes into a list and rename the variable of the DESeq2 indeces with their gene names
DEGs_without_IDs = aDRG9$GeneID[c(DEG_indeces_without_GO_IDs)]
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
  ## Run the data
Functional_Analysis(Genes = )