
library(stringr)
library(tidyverse)
library(reshape2)
library(plyr)
library(dplyr)


  ##### Data pre-processing #####
  ## Import data
MS_DF = read.csv("M:/Erik/Data/Omics/MS/20201012_Samples_View_Report_4441.csv")
  ## Re-name the column names by the first row
colnames(MS_DF) = MS_DF[1,]
  ## Delete that row
MS_DF = MS_DF[-1,]
  ## Re-order the row index
rownames(MS_DF) = NULL

  ## Replace missing integers with 0s, and question marks with blanks
for (i in 4:8){
  for (j in 1:nrow(MS_DF)){
    if (MS_DF[j,i] == ""){
      MS_DF[j,i] = 0
    } else if (MS_DF[j,i] == "?"){
      MS_DF[j,i] = ""
    }
  }
}

  ## Co-erce the data to be integers, not characters
for (i in 5:8){
  for (j in 1:nrow(MS_DF)){
    MS_DF[j,i] = as.numeric(MS_DF[j,i])
  }
}

  ## Remove the last row
MS_DF = MS_DF[-(nrow(MS_DF)),]
  ## Remove the '#' column
MS_DF = MS_DF[ , -1]
  
  ## Sort the data
MS_DF = MS_DF[order(MS_DF$myc_GFP_1, decreasing = FALSE), ]

  ## Add a column that will comprise just of human gene names
MS_DF = add_column(MS_DF, HumanGeneID = "", .before = 3)
  ## Copy the strings
MS_DF$HumanGeneID = MS_DF$HumanGeneID %>% gsub(x = MS_DF$`Identified Proteins (6608)`, replacement = "")


  ## Add a column that will comprise just of mouse gene names
MS_DF = add_column(MS_DF, MouseGeneID = "", .before = 4)
  ## Copy the strings
MS_DF$MouseGeneID = MS_DF$MouseGeneID %>% gsub(x = MS_DF$`Identified Proteins (6608)`, replacement = "")

  ## Extract the human protein's human gene name
MS_DF$HumanGeneID = str_extract_all(MS_DF$HumanGeneID, "GN=[:alnum:]{1,}")
  ## Remove the identifier
MS_DF$HumanGeneID = MS_DF$HumanGeneID %>% gsub(x = MS_DF$HumanGeneID, pattern = "GN=", replacement = "")

  ## Extract the protein's gene name
MS_DF$MouseGeneID = str_extract_all(MS_DF$MouseGeneID, "GN=[:alnum:]{1,}")
  ## Remove the identifier
MS_DF$MouseGeneID = MS_DF$MouseGeneID %>% gsub(x = MS_DF$MouseGeneID, pattern = "GN=", replacement = "")
  ## Put the gene name in MM9 EntrezID form (note, however, these were human proteins; consider leaving)
MS_DF$MouseGeneID = paste(substr(MS_DF$MouseGeneID, 1, 1),
                          tolower(substr(MS_DF$MouseGeneID, 2, 10)), sep = "")

  ## Remove extra strings from the "Identified Proteins" column
MS_DF$`Identified Proteins (6608)` = MS_DF$`Identified Proteins (6608)` %>% gsub(x = MS_DF$`Identified Proteins (6608)`, pattern = " OS.+.?$", replacement = "")

  ## Fix p53 in mouse
MS_DF$MouseGeneID[which(MS_DF$`Accession Number` == "P53_HUMAN")] = "Trp53"

  ## Find the Ratio of TMEM1:TMEM2
TMEM_RATIO = as.numeric(MS_DF[ which(MS_DF$MouseGeneID == "Tmem184b"), 8 ]) / as.numeric(MS_DF[ which(MS_DF$MouseGeneID == "Tmem184b"), 9 ])
  ## Create a column that stores the given protein's ratio of spectral counts in TMEM samples (scale the counts from TMEM 1: TMEM 2)
MS_DF = add_column(MS_DF, TMEM_RATIO = as.numeric(MS_DF$TMEM184B_myc_1) / as.numeric(MS_DF$TMEM184B_myc_2))
  ## Create a column to include whether the protein's gene was DE in aDRG TMEM mutants
MS_DF$DEG_in_aDRG = ""
  ## Remove the decoy proteins
MS_DF_filter = MS_DF %>% filter(!grepl(MS_DF$`Accession Number`, pattern = "DECOY"))
  ## Remove proteins under 10 counts in TMEM sample
MS_DF_filter1 = MS_DF_filter[ -c(which(as.numeric(MS_DF_filter$TMEM184B_myc_1) < 5)) , ]
MS_DF_filter2 = MS_DF_filter1[ -c(which(as.numeric(MS_DF_filter1$TMEM184B_myc_2) < 5)) , ]
MS_DF_filter3 = MS_DF_filter2[ -c(which(as.numeric(MS_DF_filter2$TMEM184B_myc_1) < 10)) , ]
#MS_DF_filter4 = MS_DF_filter3[ -c(which(as.numeric(MS_DF_filter3$TMEM184B_myc_2) < 10)) , ]
  ## Remove proteins that have more counts in the 2nd TMEM sample
MS_DF_filter5 = subset(MS_DF_filter3, MS_DF_filter3$TMEM_RAT > 1.0)
  ## Remove Keratin and Tubulin, common contaminants or false positives
MS_DF_filter6 = MS_DF_filter5 %>% filter(!grepl(MS_DF_filter5$`Identified Proteins (6608)` , pattern = "Keratin.+?$|Tubulin.+?$"))

MS_DF_filter6$Pauls_Likely_Candidates = ""
MS_DF_filter6$Eriks_Likely_Candidates = ""


#c(MS_DF_filter$`Identified Proteins (6608)`[which(as.numeric(MS_DF_filter$TMEM184B_myc_1) > 300 )])


DESeq2_Adults = read.csv("M:/Erik/Data/Omics/RNAseq/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
DESeq2_Adults3 = subset(DESeq2_Adults, (!is.na(DESeq2_Adults[,"AdjP"])))

  ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
DESeq2_Adults9 = DESeq2_Adults3[,] %>% filter(!grepl(DESeq2_Adults3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

#mt.+.?$|

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
DESeq2_Adults9$labs = DESeq2_Adults9$GeneID
  ## Replace DEGs with the string, "DEGs"
DESeq2_Adults9$labs[which(DESeq2_Adults9$AdjP <= 0.05)]= "DEGs"
  ## Replace the remaining genes with "Non-DEGs"
DESeq2_Adults9$labs[which(DESeq2_Adults9$AdjP >= 0.05)]= "Non-DEGs"

DEG_list = c(DESeq2_Adults9$GeneID[which(DESeq2_Adults9$labs == "DEGs")])

MS_DF_filter6$DEG_in_aDRG[which(MS_DF_filter6$MouseGeneID %in% c(DEG_list))] = "Yes"

"%notin%" = Negate("%in%")

MS_DF_filter6$DEG_in_aDRG[which(MS_DF_filter6$DEG_in_aDRG != "Yes")] = "No"

  ###### Fill "likely candidates #####

#which(MS_DF_filter6$`Accession Number` == "ABI2_HUMAN")
  ## On Paul's, not on mine: NINL, GATD1, MTG8R, ARI3B, SPNDC, MTCH1, PDL17, MYOD1, SUCA
MS_DF_filter6$Pauls_Likely_Candidates[c(15,107,65,29,43,17,9,30,74,22,35,76,31,32,77,26,134,51,11,98,84,85,86,116,38,123,124,145,89,149,155,166,204,152,167,262,274,528,542,557,544,533,546,579,617,572,618,564,565,567,541,539,586,560,663,645,822,777,817,798,824,793,797,806,890,854,847,837,856,809,839,859,929,967,1013,971,1001,1003,991,995,996,1020,1052,1040,1008,1043,1087,993,1130,1115,1140,1145,1116,1149,1172,1177,1178,1232,1237,1239,1235,1327,1347,1348,1336,1345,1358,1346,1355,1339,1402,1408,1452,1459,1464,297,313,340,377)] = MS_DF_filter6$HumanGeneID[c(15,107,65,29,43,17,9,30,74,22,35,76,31,32,77,26,134,51,11,98,84,85,86,116,38,123,124,145,89,149,155,166,204,152,167,262,274,528,542,557,544,533,546,579,617,572,618,564,565,567,541,539,586,560,663,645,822,777,817,798,824,793,797,806,890,854,847,837,856,809,839,859,929,967,1013,971,1001,1003,991,995,996,1020,1052,1040,1008,1043,1087,993,1130,1115,1140,1145,1116,1149,1172,1177,1178,1232,1237,1239,1235,1327,1347,1348,1336,1345,1358,1346,1355,1339,1402,1408,1452,1459,1464,297,313,340,377)]


MS_DF_filter6$Eriks_Likely_Candidates[c(6,7,9,10,11,13,15:20,22,24,26,27,29:31,36,38,40,42,43,47,48,51,56,63,69,71,74,76,77,80,83:91,96,97,98,100,102,105,107,108,109,113:117,123,124,126,127,132,134,136,138,141,143,144,145,149,152,154:156,158,159,166,167,172:177,183,184,186,190,193:196,198,203,204,209,212,219,222:227,230,231,236:240,244,247,259,297,308,312,313,316,340,343,347,349,353,372,380,403,439,471,493,501,521,527,531,533,535,539,541,542,546,547,550,556:561,563,566,568,570,573:575,579,590:593,599:601,607,612,614,620,643:648,657,673,674,691,764,775,777,789,793,797:799,804:806,809,811,815,817,819,822,824,833,839,841,843,846,849,857,861,868,874,876,879,885,914,918,932,967,972,980,985,991:996,998,1001,1008,1010,1011,1013,1014,1016:1018,1020,1029,1040,1047,1058,1078,1115,1116,1127,1140,1143,1145,1149,1154,1172,1173,1177:1180,1185,1205,1232,1235,1237,1239,1256,1258,1261,1263,1269,1272,1291,1292,1301,1322,1327,1336,1337,1339,1342,1345:1347,1355,1358,1367,1374,1377,1381,1476,1485,1490)] = MS_DF_filter6$HumanGeneID[c(6,7,9,10,11,13,15:20,22,24,26,27,29:31,36,38,40,42,43,47,48,51,56,63,69,71,74,76,77,80,83:91,96,97,98,100,102,105,107,108,109,113:117,123,124,126,127,132,134,136,138,141,143,144,145,149,152,154:156,158,159,166,167,172:177,183,184,186,190,193:196,198,203,204,209,212,219,222:227,230,231,236:240,244,247,259,297,308,312,313,316,340,343,347,349,353,372,380,403,439,471,493,501,521,527,531,533,535,539,541,542,546,547,550,556:561,563,566,568,570,573:575,579,590:593,599:601,607,612,614,620,643:648,657,673,674,691,764,775,777,789,793,797:799,804:806,809,811,815,817,819,822,824,833,839,841,843,846,849,857,861,868,874,876,879,885,914,918,932,967,972,980,985,991:996,998,1001,1008,1010,1011,1013,1014,1016:1018,1020,1029,1040,1047,1058,1078,1115,1116,1127,1140,1143,1145,1149,1154,1172,1173,1177:1180,1185,1205,1232,1235,1237,1239,1256,1258,1261,1263,1269,1272,1291,1292,1301,1322,1327,1336,1337,1339,1342,1344:1347,1355,1358,1367,1374,1377,1381,1476,1485,1490)]




BP_Candidate_List_Mouse = MS_DF_filter6$MouseGeneID[c(6,7,9,10,11,13,15:20,22,24,26,27,29:31,36,38,40,42,43,47,48,51,56,63,69,71,74,76,77,80,83:91,96,97,98,100,102,105,107,108,109,113:117,123,124,126,127,132,134,136,138,141,143,144,145,149,152,154:156,158,159,166,167,172:177,183,184,186,190,193:196,198,203,204,209,212,219,222:227,230,231,236:240,244,247,259,297,308,312,313,316,340,343,347,349,353,372,380,403,439,471,493,501,521,527,531,533,535,539,541,542,546,547,550,556:561,563,566,568,570,573:575,579,590:593,599:601,607,612,614,620,643:648,657,673,674,691,764,775,777,789,793,797:798,804:806,809,811,815,817,819,822,824,833,839,841,843,846,849,854,857,861,868,874,876,879,885,914,918,932,967,972,980,985,991:996,998,1001,1008,1010,1011,1013,1014,1016:1018,1020,1029,1040,1047,1058,1078,1115,1116,1127,1140,1143,1145,1149,1154,1172,1173,1177:1180,1185,1232,1235,1237,1239,1256,1258,1261,1263,1269,1272,1291,1292,1301,1322,1327,1336,1337,1339,1342,1344:1347,1355,1358,1367,1374,1377,1381,1476,1485,1490)]

BP_Candidate_List_Human = MS_DF_filter6$HumanGeneID[c(6,7,9,10,11,13,15:20,22,24,26,27,29:31,36,38,40,42,43,47,48,51,56,63,69,71,74,76,77,80,83:91,96,97,98,100,102,105,107,108,109,113:117,123,124,126,127,132,134,136,138,141,143,144,145,149,152,154:156,158,159,166,167,172:177,183,184,186,190,193:196,198,203,204,209,212,219,222:227,230,231,236:240,244,247,259,297,308,312,313,316,340,343,347,349,353,372,380,403,439,471,493,501,521,527,531,533,535,539,541,542,546,547,550,556:561,563,566,568,570,573:575,579,590:593,599:601,607,612,614,620,643:648,657,673,674,691,764,775,777,789,793,797:798,804:806,809,811,815,817,819,822,824,833,839,841,843,846,849,854,857,861,868,874,876,879,885,914,918,932,967,972,980,985,991:996,998,1001,1008,1010,1011,1013,1014,1016:1018,1020,1029,1040,1047,1058,1078,1115,1116,1127,1140,1143,1145,1149,1154,1172,1173,1177:1180,1185,1232,1235,1237,1239,1256,1258,1261,1263,1269,1272,1291,1292,1301,1322,1327,1336,1337,1339,1342,1344:1347,1355,1358,1367,1374,1377,1381,1476,1485,1490)]


##### Functional Analysis Prep #####
library(PANTHER.db)
pthOrganisms(PANTHER.db) = "MOUSE"

library(enrichR)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(GO.db)

  ## Connect live to the Enrichr server/master database (website) and store the "space" as a variable. This contains a vast array of bioinformatics databases/websites.
    ## Store it as a variable for better visualization and eventual subsetting
DBs = listEnrichrDbs()

  ## View the variable in the editor to find the relevant databases and their indeces for subsetting (click on the dataframe in the Global Environment and manually peruse)

  ## e.g. "GO_Biological_Process_2018" (index # 130), and "Panther_2016" (index # 102)
    ## Concatenate them in a list for subsetting, or slice directly
DBs$libraryName[109]
diff_DBs = DBs$libraryName[c(50,127,91,92,109)]
DBs = DBs$libraryName[c(130,131,132,102,14,110)]

  ## Find human Entrez IDs (Ensembl IDs) from the list
Human_Ensembl_BP_IDs = as.integer(
  mapIds(org.Hs.eg.db,
         keys = as.character(
           c(BP_Candidate_List_Human[])
         ),
         keytype = "SYMBOL",
         column = "ENTREZID")
)
  ## Extract GO Terms for all human genes in the list
Human_BP_GO_IDs = c(
  mapIds(org.Hs.eg.db, 
         keys = as.character(
           c(BP_Candidate_List_Human[])
         ),
         keytype = "SYMBOL",
         column = "GOALL",
         multiVals = "list")
)

  ## Find mouse Entrez IDs (Ensembl IDs) from the list
Mouse_Ensembl_BP_IDs = as.integer(
  mapIds(org.Mm.eg.db,
         keys = as.character(
           c(BP_Candidate_List_Mouse[])
         ),
         keytype = "SYMBOL",
         column = "ENTREZID")
)
  ## Extract GO Terms for all mouse genes in the list
Mouse_BP_GO_IDs = c(
  mapIds(org.Mm.eg.db, 
         keys = as.character(
           c(BP_Candidate_List_Mouse[])
         ),
         keytype = "SYMBOL",
         column = "GOALL",
         multiVals = "list")
)


## Find the human genes that don't have Ensembl/Entrez IDs
Human_BP_indeces_without_Ensembl_IDs = c(
  which(is.na(Human_Ensembl_BP_IDs[]) == TRUE)
)
## Find the human genes that don't have GO IDs
Human_BP_indeces_without_GO_IDs = as.integer(c(
  which(is.na(Human_BP_GO_IDs[]) == TRUE)
))


  ## Find the genes without Ensembl/Entrez IDs
which(is.na(Mouse_Ensembl_BP_IDs[]) == TRUE)
Mouse_BP_indeces_without_Ensembl_IDs = c(
  which(is.na(Mouse_Ensembl_BP_IDs[]) == TRUE)
)
  ## Repeat for GO IDs
Mouse_BP_indeces_without_GO_IDs = as.integer(c(
  which(is.na(Mouse_BP_GO_IDs[]) == TRUE)
))



MS_DF_filter6$HumanGeneID[c(Human_BP_indeces_without_Ensembl_IDs, Human_BP_indeces_without_GO_IDs)]
  ## Concatenate the proteins' gene names into a list and rename the variable of the MS indeces with their gene names
Human_BPs_without_IDs = MS_DF_filter6$HumanGeneID[c(Human_BP_indeces_without_GO_IDs)]
  ## Use the indeces of the BP_indeces_without_GO_IDs for looping if statement in the creation of the data frame below
#BP_indeces_without_GO_IDs

#which(BP_indeces_without_Ensembl_IDs == BP_indeces_without_GO_IDs)
#match(BP_indeces_without_Ensembl_IDs, BP_indeces_without_GO_IDs)

  ## Create a matrix to house GO information for each gene
GO_INFO = matrix(nrow = length(c(BP_Candidate_List_Human)), ncol = 3)
rownames(GO_INFO) = c(BP_Candidate_List_Human)
colnames(GO_INFO) = c("# of GO Terms", "GO Term ID", "GO Term")

  ## Fill the first column with zeros for genes without GO Terms
GO_INFO[c(Human_BP_indeces_without_GO_IDs), 1] = 0
  ## Fill the second and third columns with "None"s for genes without GO Terms
GO_INFO[c(Human_BP_indeces_without_GO_IDs), c(2,3)] = "None"
  ## Fill the second column with the concatenated GO Term IDs associated with each gene, skipping over the genes with 0
for (i in 1:length(Human_BP_GO_IDs)){
  if (i == 119) {
    next
  }
  if (i == 225) {
    next
  }
  if (i == 246) {
    next
  }
  GO_INFO[i,2] = paste(names(mapIds(GO.db, Human_BP_GO_IDs[[i]][], "TERM", "GOID")), collapse = ";")
}
  ## Extract the GO IDs and remove the redundant ones;
    ## Fill the first column with the numbers of GO Terms associated with each gene, skipping over the genes with 0
for (i in 1:length(Human_BP_GO_IDs)){
  if (i == 119) {
    next
  }
  if (i == 225) {
    next
  }
  if (i == 246) {
    next
  }
  x = str_split(as.vector(GO_INFO[i,2]), pattern = ";", simplify = TRUE)
  x = c(as.character(x[which(x %in% x == TRUE)]))
  GO_INFO[i,2] = paste(unique(x), collapse = ";")
  GO_INFO[i,1] = paste(as.numeric(length(unique(x))))
}
  ## Fill the third column with the concatenated GO Terms associated with each gene, skipping over the genes with 0
for (i in 1:length(Human_BP_GO_IDs)){
  if (i == 119) {
    next
  }
  if (i == 225) {
    next
  }
  if (i == 246) {
    next
  }
  x = str_split(as.vector(GO_INFO[i,2]), pattern = ";", simplify = TRUE)
  GO_INFO[i,3] = paste(as.character(mapIds(GO.db, keys = x, keytype = "GOID", "TERM")), collapse = ";")
}

  ## Extract all of the terms into one vector to find unique GO terms
for (i in 1:length(Human_BP_GO_IDs)){
  if (i == 119) {
    next
  }
  if (i == 225) {
    next
  }
  if (i == 246) {
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
Functional_Analysis(Genes = BP_Candidate_List_Human)

#Enrichr_Functional_Output = PPIs[,c(1,4:8,10)]
                                  #GO_Cell_Comps[1:206,c(1,4:8,10)],
                                  #GO_Mol_Funcs[1:206,c(1,4:8,10)],
                                  #PATHWAYS[1:206,c(1,4:8,10)],
                                  #PPIs[1:206,c(1,4:8,10)])
#colnames(Enrichr_Functional_Output) = c("Protein-Protein Interaction Hub Proteins","Weighted Overlap Ratio","Pval","Qval*EnrichrZ","Qval","EnrichrZ","Proteins Involved")

#write.csv(Enrichr_Functional_Output, file = "M:\\Erik\\Data\\Omics\\MS\\TMEM_mycIP_BPCandidates_PPIs.csv", sep = ",", row.names = F, col.names = T)


##### Find terms in other genes/proteins #####
  ## Find the GO Terms of a desired gene/protein
Find_GOs_of_Gene_X = function(GENE_ID) {
  return(str_split(
    as.vector(GO_INFO[ which(GO_INFO_df$GeneID == GENE_ID) ,   3]), pattern = ";", simplify = TRUE)[
      which(str_split(
        as.vector(GO_INFO[ which(GO_INFO_df$GeneID == GENE_ID) ,    3]), pattern = ";", simplify = TRUE) %in% Unique_GOs)
    ]
  )
}
Find_GOs_of_Gene_X(GENE_ID = "JUN")

  ## Search for whether GO terms of one gene/protein of interest are also of another
Find_GOs_of_Two_Genes = function(GENE_ID1, GENE_ID2) {
  TermIdx = c(which(  
    str_split(
      as.vector(GO_INFO_df[ which(GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
    )[which(str_split(
      as.vector(GO_INFO_df[ which(GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
    ) %in% Unique_GOs)] %in%
      
      str_split(
        as.vector(GO_INFO_df[ which(GO_INFO_df$GeneID == GENE_ID2), 4]), pattern = ";", simplify = TRUE
      )[which(str_split(
        as.vector(GO_INFO_df[ which(GO_INFO_df$GeneID == GENE_ID2), 4]), pattern = ";", simplify = TRUE
      ) %in% Unique_GOs)]
  ))
  ## Find which terms those are
  str_split(
    as.vector(GO_INFO_df[ which(GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
  )[which(str_split(
    as.vector(GO_INFO_df[ which(GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
  ) %in% Unique_GOs)][TermIdx]
  return(str_split(
    as.vector(GO_INFO_df[ which(GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
  )[which(str_split(
    as.vector(GO_INFO_df[ which(GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
  ) %in% Unique_GOs)][TermIdx])
}
Find_GOs_of_Two_Genes(GENE_ID1 = "MTOR", GENE_ID2 = "JUN")
  
  ## Find the genes/proteins that contain a specific GO term(s)
Find_Genes_With_X_GO_Term = function(GO_Term){
  return(GO_INFO_df$GeneID[which(
    apply(GO_INFO_df, 1, function(x) any(grepl(GO_Term, x)))
  )])
}
Find_Genes_With_X_GO_Term(GO_Term = "generation of neurons")

  ## Find genes/proteins overlapping in multiple GO Terms
Find_Genes_Related_By_GO_Term = function(GO_Term1, GO_Term2){
  MatchIdx = c(which(GO_INFO_df$GeneID[which(
    apply(GO_INFO_df, 1, function(x) any(grepl(GO_Term1, x)))
  )] %in% GO_INFO_df$GeneID[which(
    apply(GO_INFO_df, 1, function(x) any(grepl(GO_Term2, x)))
  )]))
  ## Find which genes/proteins those are
  GO_INFO_df$GeneID[which(
    apply(GO_INFO_df, 1, function(x) any(grepl(GO_Term1, x)))
  )][MatchIdx]
  return(GO_INFO_df$GeneID[which(
    apply(GO_INFO_df, 1, function(x) any(grepl(GO_Term1, x)))
  )][MatchIdx])
}
Find_Genes_Related_By_GO_Term(GO_Term1 = "generation of neurons", GO_Term2 = "response to axon injury")



  ##### Repeat Analysis in Mouse #####

GO_INFO = matrix(nrow = length(c(BP_Candidate_List_Mouse)), ncol = 3)
rownames(GO_INFO) = c(BP_Candidate_List_Mouse)
colnames(GO_INFO) = c("# of GO Terms", "GO Term ID", "GO Term")

  ## Fill the first column with zeros for genes without GO Terms
GO_INFO[c(Mouse_BP_indeces_without_GO_IDs), 1] = 0
  ## Fill the second and third columns with "None"s for genes without GO Terms
GO_INFO[c(Mouse_BP_indeces_without_GO_IDs), c(2,3)] = "None"
  ## Fill the second column with the concatenated GO Term IDs associated with each gene, skipping over the genes with 0
for (i in 1:length(Mouse_BP_GO_IDs)){
  if (i == 3) {
    next
  }
  if (i == 6) {
    next
  }
  if (i == 28) {
    next
  }
  if (i == 77) {
    next
  }
  if (i == 84) {
    next
  }
  if (i == 135) {
    next
  }
  if (i == 172) {
    next
  }
  if (i == 246) {
    next
  }
  if (i == 268) {
    next
  }
  GO_INFO[i,2] = paste(names(mapIds(GO.db, BP_GO_IDs[[i]][], "TERM", "GOID")), collapse = ";")
}
  ## Extract the GO IDs and remove the redundant ones;
    ## Fill the first column with the numbers of GO Terms associated with each gene, skipping over the genes with 0
for (i in 1:length(Mouse_BP_GO_IDs)){
  if (i == 3) {
    next
  }
  if (i == 6) {
    next
  }
  if (i == 28) {
    next
  }
  if (i == 77) {
    next
  }
  if (i == 84) {
    next
  }
  if (i == 135) {
    next
  }
  if (i == 172) {
    next
  }
  if (i == 246) {
    next
  }
  if (i == 268) {
    next
  }
  x = str_split(as.vector(GO_INFO[i,2]), pattern = ";", simplify = TRUE)
  x = c(as.character(x[which(x %in% x == TRUE)]))
  GO_INFO[i,2] = paste(unique(x), collapse = ";")
  GO_INFO[i,1] = paste(as.numeric(length(unique(x))))
}
  ## Fill the third column with the concatenated GO Terms associated with each gene, skipping over the genes with 0
for (i in 1:length(Mouse_BP_GO_IDs)){
  if (i == 3) {
    next
  }
  if (i == 6) {
    next
  }
  if (i == 28) {
    next
  }
  if (i == 77) {
    next
  }
  if (i == 84) {
    next
  }
  if (i == 135) {
    next
  }
  if (i == 172) {
    next
  }
  if (i == 246) {
    next
  }
  if (i == 268) {
    next
  }
  x = str_split(as.vector(GO_INFO[i,2]), pattern = ";", simplify = TRUE)
  GO_INFO[i,3] = paste(as.character(mapIds(GO.db, keys = x, keytype = "GOID", "TERM")), collapse = ";")
}

  ## Extract all of the terms into one vector to find unique GO terms
for (i in 1:length(Mouse_BP_GO_IDs)){
  if (i == 3) {
    next
  }
  if (i == 6) {
    next
  }
  if (i == 28) {
    next
  }
  if (i == 77) {
    next
  }
  if (i == 84) {
    next
  }
  if (i == 135) {
    next
  }
  if (i == 172) {
    next
  }
  if (i == 246) {
    next
  }
  if (i == 268) {
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


Functional_Analysis(Genes = BP_Candidate_List_Mouse)

Find_GOs_of_Gene_X(GENE_ID = "Jun")
Find_GOs_of_Two_Genes(GENE_ID1 = "Mtor", GENE_ID2 = "Jun")
Find_Genes_With_X_GO_Term(GO_Term = "generation of neurons")
Find_Genes_Related_By_GO_Term(GO_Term1 = "generation of neurons", GO_Term2 = "response to axon injury")
