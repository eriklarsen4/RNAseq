


  ## Script for obtaining Enrichr-based bioinformatics data
    ## Developed by Erik Larsen

    ## Condense lines of code by clicking all arrows next to code line #. Click again to open, step-wise to maintain organization

##### Environment Prep #####

  ## Install packages to use the "enrichr", "panther" packages
    ## At least for the first installation, use the following code; THIS REQUIRES R version 4.0.0+ !!!

  ## Uncomment the following if/install statements if these packages have not previously been installed, and install them.
    #if (!requireNamespace("BiocManager", quietly = TRUE))
    #  install.packages("BiocManager")
    #BiocManager::install(version = "3.11")

    #if (!requireNamespace("BiocManager", quietly = TRUE))
    #  install.packages("BiocManager")
    #BiocManager::install("PANTHER.db")

    #if (!requireNamespace("BiocManager", quietly = TRUE))
    #  install.packages("BiocManager")
    #BiocManager::install("EnhancedVolcano")

  ## Load the Panther 2016 Pathway Database package and use "mouse" as reference genome/specified species for downstream bioinformatics analysis
library(PANTHER.db)
pthOrganisms(PANTHER.db) = "MOUSE"
  ## (Don't update, if the option pops up
# n)

## Load the enrichr package
library(enrichR)
## Load previously installed packages from the library (if not previously installed, install!)
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
library(ggsci)
library(viridis)
library(dendextend)
library(ggdendro)
library(grid)
#library(EnhancedVolcano)
library(pheatmap)
#library(RColorBrewer)

  ## Connect live to the Enrichr server/master database (website) and store the "space" as a variable. This contains a vast array of bioinformatics databases/websites
    ## Store it as a variable for better visualization and eventual subsetting
DBs = listEnrichrDbs()

  ## View the variable in the editor to find the relevant databases and their indeces for subsetting (click on the dataframe in the global environment and manually peruse)

  ## In this case, "GO_Biological_Process_2018" (index # 130), and "Panther_2016" (index # 102)
    ## Concatenate them in a list for subsetting, or slice directly (we'll slice directly in the next line)
DBs = DBs$libraryName[c(130,102)]

##### Young Genotype Data prep #####
  ## Import the dataset you want to analyze
DESeq2_Hippo_YWvYM = read.csv("M:/Erik/Data/Omics/Hippocampus/Genotype Comparisons/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results for Young Hippocampus without ymu 4141.csv")
  ## Filter (subset) ribosomal RNAs, Riks, and genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns.
DESeq2_Hippo_YWvYM3 = subset(DESeq2_Hippo_YWvYM, (!is.na(DESeq2_Hippo_YWvYM[,"AdjP"])))

  ## Filter the rRNAs
DESeq2_Hippo_YWvYM3 = DESeq2_Hippo_YWvYM3 %>% filter(!grepl(DESeq2_Hippo_YWvYM3$GeneID, pattern = "Rpl.+.?$"))
  ## Again
DESeq2_Hippo_YWvYM3 = DESeq2_Hippo_YWvYM3 %>% filter(!grepl(DESeq2_Hippo_YWvYM3$GeneID, pattern = "Rps.+.?$"))
  ## Again
DESeq2_Hippo_YWvYM3 = DESeq2_Hippo_YWvYM3 %>% filter(!grepl(DESeq2_Hippo_YWvYM3$GeneID, pattern = "Mrpl.+.?$"))
  ## Again
DESeq2_Hippo_YWvYM3 = DESeq2_Hippo_YWvYM3 %>% filter(!grepl(DESeq2_Hippo_YWvYM3$GeneID, pattern = "Mrps.+.?$"))
  ## Filter the unrecognized genes
DESeq2_Hippo_YWvYM3 = DESeq2_Hippo_YWvYM3 %>% filter(!grepl(DESeq2_Hippo_YWvYM3$GeneID, pattern = ".*Rik$"))

  ## Export these Genes to the server for Enrichr upload to obtain all the genes in the relevant pathways/terms
#write.table(DESeq2_Hippo_YWvYM3, "M:\\Erik\\Data\\Omics\\Hippocampus\\Genotype Comparisons\\Enrichr DESeq2 df Young wo ymu4141.csv", sep = ",", col.names = T)


  ## Perform relevant analysis through Enrichr; in this case, FDR = 0.05
    ## first: GO Bio Processes; this analyzes genes from your dataset in Enrichr, which draws terms from geneontology.org
Hippo_YWvYM_GOs = as.data.frame(
  enrichr(
    c(DESeq2_Hippo_YWvYM3$GeneID[which(DESeq2_Hippo_YWvYM3$AdjP < 0.05)]), DBs[1])
)
    ## second: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
Hippo_YWvYM_PATHWAYS = as.data.frame(
  enrichr(
    c(DESeq2_Hippo_YWvYM3$GeneID[which(DESeq2_Hippo_YWvYM3$AdjP < 0.05)]), DBs[2])
)

  ## Clean up both dataframes and export both dataframes.
  ## Dataframes for binary heatmap visualization (presence in pathway/term or not) need additional cleaning in Python;
  ## Following that manipulation, data will be re-imported for visualization
Hippo_YWvYM_GOs = Hippo_YWvYM_GOs %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Hippo_YWvYM_GOs$Overlap = as.numeric(Hippo_YWvYM_GOs$Overlap)
Hippo_YWvYM_GOs$`Process Size` = as.numeric(Hippo_YWvYM_GOs$`Process Size`)
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Hippo_YWvYM_GOs = add_column(Hippo_YWvYM_GOs, EnrichrZscore = Hippo_YWvYM_GOs$GO_Biological_Process_2018.Combined.Score/log(Hippo_YWvYM_GOs$GO_Biological_Process_2018.P.value), .before = 4)


  ## Vectorize and make the DEGs within GO Biological Processes compatible with the gene names output from DESeq2 (EntrezID; capital first letter, lowercase all rest)
    ## Use for subsetting in downstream visualizations (indexing and labeling)
Hippo_YWvYM_GO_GENES = str_split(Hippo_YWvYM_GOs$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Hippo_YWvYM_GO_GENES)){
  Hippo_YWvYM_GO_GENES[i,] = paste(substr(Hippo_YWvYM_GO_GENES[i,], 1, 1),
                           tolower(substr(Hippo_YWvYM_GO_GENES[i,], 2, 7)), sep = "") 
}

  ## Repeat same process for PATHWAYS
Hippo_YWvYM_PATHWAYS = Hippo_YWvYM_PATHWAYS %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Hippo_YWvYM_PATHWAYS$Overlap = as.numeric(Hippo_YWvYM_PATHWAYS$Overlap)
Hippo_YWvYM_PATHWAYS$`Process Size` = as.numeric(Hippo_YWvYM_PATHWAYS$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Hippo_YWvYM_PATHWAYS = add_column(Hippo_YWvYM_PATHWAYS, EnrichrZscore = Hippo_YWvYM_PATHWAYS$Panther_2016.Combined.Score/log(Hippo_YWvYM_PATHWAYS$Panther_2016.P.value), .before = 4)
  ## Re-order the dataframe by Q value
#Hippo_YWvYM_PATHWAYS = Hippo_YWvYM_PATHWAYS[order(Hippo_YWvYM_PATHWAYS$Panther_2016.Adjusted.P.value, decreasing = FALSE),]
  ## Re-set the row index to make sure indeces are correct
#row.names(Hippo_YWvYM_PATHWAYS) = NULL

  ## Vectorize and make the DEGs within GO Biological Processes compatible with the gene names output from DESeq2 (EntrezID; capital first letter, lowercase all rest)
    ## Use for subsetting in downstream visualizations (indexing and labeling)
Hippo_YWvYM_PATHWAY_GENES = str_split(Hippo_YWvYM_PATHWAYS$Panther_2016.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Hippo_YWvYM_PATHWAY_GENES)){
  Hippo_YWvYM_PATHWAY_GENES[i,] = paste(substr(Hippo_YWvYM_PATHWAY_GENES[i,], 1, 1),
                                tolower(substr(Hippo_YWvYM_PATHWAY_GENES[i,], 2, 7)), sep = "") 
}

  ## Remove unused columns ("old p-values" and "odds ratio")
Hippo_YWvYM_GOs = Hippo_YWvYM_GOs[,-c(7,8,9)]
Hippo_YWvYM_PATHWAYS = Hippo_YWvYM_PATHWAYS[,-c(7,8,9)]

  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ##  --->>> The fraction of DEGs belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Hippo_YWvYM_GOs = add_column(Hippo_YWvYM_GOs, Weighted_Overlap_Ratio = Hippo_YWvYM_GOs$Overlap*(Hippo_YWvYM_GOs$Overlap/Hippo_YWvYM_GOs$`Process Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
Hippo_YWvYM_GOs = add_column(Hippo_YWvYM_GOs, Modified_Combined_Score = (Hippo_YWvYM_GOs$EnrichrZscore)*(log10(Hippo_YWvYM_GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)
  ## Re-order the dataframe by Modified Combined Score
Hippo_YWvYM_GOs = Hippo_YWvYM_GOs[order(Hippo_YWvYM_GOs$Modified_Combined_Score, decreasing = TRUE),]
  ## Re-set the row index to make sure indeces are correct
row.names(Hippo_YWvYM_GOs) = NULL


  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ##  --->>> The fraction of DEGs belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Hippo_YWvYM_PATHWAYS = add_column(Hippo_YWvYM_PATHWAYS, Weighted_Overlap_Ratio = Hippo_YWvYM_PATHWAYS$Overlap*(Hippo_YWvYM_PATHWAYS$Overlap/Hippo_YWvYM_PATHWAYS$`Process Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
Hippo_YWvYM_PATHWAYS = add_column(Hippo_YWvYM_PATHWAYS, Modified_Combined_Score = (Hippo_YWvYM_PATHWAYS$EnrichrZscore)*(log10(Hippo_YWvYM_PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)
  ## Re-order the dataframe by Modified Combined Score
Hippo_YWvYM_PATHWAYS = Hippo_YWvYM_PATHWAYS[order(Hippo_YWvYM_PATHWAYS$Modified_Combined_Score, decreasing = TRUE),]
  ## Re-set the row index to make sure indeces are correct
row.names(Hippo_YWvYM_PATHWAYS) = NULL



  ## Find genes within a Process/Pathway of interest that are DE (down) in mutants
    ## For example, find which genes that are differentially downregulated in young hippocampus mutants
  ## Concatenate a list of all the overlapping DEGs within the process
Postsyn_memb_org_Overlap = c(which(DESeq2_Hippo_YWvYM3$GeneID %in% Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]]))
  ## Concatenate a different list of the indeces from the overlapping genes that are up/down in mutants
Up_Genes_indeces_in_Postsyn_memb_org = c(which(DESeq2_Hippo_YWvYM3[c(Postsyn_memb_org_Overlap),3] > 0))
  ## Find which genes they are to manually confirm
DESeq2_Hippo_YWvYM3$GeneID[Postsyn_memb_org_Overlap[Up_Genes_indeces_in_Postsyn_memb_org]]


## Export data as .CSVs to the server
#write.table(Hippo_YWvYM_GOs, "M:\\Erik\\Data\\Omics\\Hippocampus\\Genotype Comparisons\\5mo Hippocampus GT woymu4141 GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(Hippo_YWvYM_PATHWAYS, "M:\\Erik\\Data\\Omics\\Hippocampus\\Genotype Comparisons\\5mo Hippocampus GT woymu4141 Panther Pathways.csv", sep = ",", col.names = T)


  ## Create a column in the DESeq2 dataframe that scales the Adjusted P-value by log-base 10
DESeq2_Hippo_YWvYM3$log10Pval = -log10(DESeq2_Hippo_YWvYM3$AdjP)

  ## Create the data points of interest for functional display.. ggplot struggles to plot subsets of dataframes, so creating a variable that directly labels all the points accordingly helps.
    ## Create a new column (genes of interest); Using "GeneID" is irrelevant; the values will be replaced in the next step
DESeq2_Hippo_YWvYM3$g.o.i. = DESeq2_Hippo_YWvYM3$GeneID

## Turn the concatenated list of genes, COPIED AND PASTED FROM ENRICHR, into compatible Entrez GeneIDs (lowercase all non-first letters)
## This reference list can be used for future plotting/labeling in volcanoes, MAs, and heatmaps.
Postsynaptic_Membrane_Organization = c("NLGN2","NLGN4X","NLGN3","NLGN4Y","NLGN1","NRXN2","NRXN1","CEL","LRP4","GPHN","DLG4","DNAJA3","APOE","MUSK","CHRNB1","SHANK3")

for (i in 1:length(Postsynaptic_Membrane_Organization)){
  Postsynaptic_Membrane_Organization[i] = paste(substr(Postsynaptic_Membrane_Organization[i], 1, 1),
                                    tolower(substr(Postsynaptic_Membrane_Organization[i], 2, 8)), sep = "") 
}
  ## Concatenate the strings of the DEGs in the Young Hippocampus DESeq2 dataset that are also in this process
c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]])
  ## Find the indeces of these genes in the DESeq2 dataframe for plotting purposes
which(DESeq2_Hippo_YWvYM3$GeneID %in% c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]]))
## Label these points in the volcano/MA plots

## Create a "not in" operator
"%ni%" = Negate("%in%")
which(c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]]) %ni% c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$`Process Size`[1]]))

  ## Indeces from the DESeq2 dataframe of genes NOT overlapping in the GO Bio Process, "Postsynaptic Membrane Organization"
which(DESeq2_Hippo_YWvYM3$GeneID %ni% c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]]))
  ## (the genes themselves)
DESeq2_Hippo_YWvYM3[which(DESeq2_Hippo_YWvYM3$GeneID %ni% c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]])),1]
  ## Concatenated list of the overlapping DEGs from the DESeq2 dataset and the "Postsynaptic Membrane Organization" GO
DEG_Overlap = c(66,281,290,521,692,816,917)
  ## Subsetted list of the DEGs of all the genes in the GO
which(DESeq2_Hippo_YWvYM3$GeneID %in% c(Postsynaptic_Membrane_Organization))
  ## Subsetted list of the genes in the GO not DE
Postsynaptic_Membrane_Organization[-DEG_Overlap]




##### Volcano Plot of 5mo GT Hippocampus #####

  ## Fill the points with appropriately indexed data
    ## "Non-differentially Expressed Genes" are defined as being "genes of interest" with Adjusted P-values > 0.01
DESeq2_Hippo_YWvYM3$g.o.i.[which(DESeq2_Hippo_YWvYM3$AdjP > 0.05)] = "Non-DEGs"
    ## "Non-differentially expressed genes" identified in the "Postsynaptic Membrane Organization" GO Biological Process
DESeq2_Hippo_YWvYM3$g.o.i.[which(DESeq2_Hippo_YWvYM3$GeneID %in% c(Postsynaptic_Membrane_Organization))] = "Genes in\nPostsynaptic\nMembrane\nOrganization GO"
    ## "Differentially Expressed Genes" also identified in Postsynaptic Membrane Organization GO ("overlapping" genes)
#DESeq2_Hippo_YWvYM3$g.o.i.[c(which(DESeq2_Hippo_YWvYM3$GeneID %in% c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]])))] = "DEGs in\nPostsynaptic\nMembrane\nOrganization"
  ## "Differentially Expressed Genes" are defined as the genes differentially expressed (AdjP < 0.01) but not overlapping with the "Postsynaptic Membrane Organization" GO genes
DESeq2_Hippo_YWvYM3$g.o.i.[c(which(DESeq2_Hippo_YWvYM3$AdjP <= 0.05))][-c(which(DESeq2_Hippo_YWvYM3$GeneID %in% c(Postsynaptic_Membrane_Organization)))] = "DEGs"

  ## Turn the plotting column into a factor vector
DESeq2_Hippo_YWvYM3$g.o.i. = as.factor(DESeq2_Hippo_YWvYM3$g.o.i.)
  ## Check to make sure
class(DESeq2_Hippo_YWvYM3$g.o.i.)
  ## Check to make sure all the levels exist and didn't get overwritten by another/were properly indexed
levels(DESeq2_Hippo_YWvYM3$g.o.i.)

  ## Create labels on points we're interested in labeling in volcano/MA plots
    ## Create the column, "labs"; hide all of the text labels with: ""
DESeq2_Hippo_YWvYM3$labs = ""

  ## Figure out the indeces of genes we're interested in labeling
which(DESeq2_Hippo_YWvYM3$GeneID %in% c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]]))
Hippo_YWvYM_GO_GENES[1,1:6]
which(DESeq2_Hippo_YWvYM3$GeneID == "Nlgn3") # 816
which(DESeq2_Hippo_YWvYM3$GeneID == "Nlgn2") # 917
which(DESeq2_Hippo_YWvYM3$GeneID == "Dlg4") # 281
which(DESeq2_Hippo_YWvYM3$GeneID == "Nrxn2") # 290
which(DESeq2_Hippo_YWvYM3$GeneID == "Lrp4") # 66
which(DESeq2_Hippo_YWvYM3$GeneID == "Apoe") # 692
which(DESeq2_Hippo_YWvYM3$GeneID == "Shank3") # 521


  ## Fill the appropriate rows of the "labs" column with the corresponding Gene names for labeling those points
DESeq2_Hippo_YWvYM3$labs[c(which(DESeq2_Hippo_YWvYM3$GeneID %in% c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]])))] = DESeq2_Hippo_YWvYM3$GeneID[c(which(DESeq2_Hippo_YWvYM3$GeneID %in% c(Hippo_YWvYM_GO_GENES[1,1:Hippo_YWvYM_GOs$Overlap[1]])))]

  ## Create another "labels" column and fill with cherry-picked genes
DESeq2_Hippo_YWvYM3$labs2 = ""
DESeq2_Hippo_YWvYM3$labs2[c(2,816,917,281,290,66,692,521)] = DESeq2_Hippo_YWvYM3$GeneID[c(2,816,917,281,290,66,692,521)]

  ## Notes on ggplot commands:
  ## geom-point = alpha is a ggplot add-on function that controls point color transparency
  ## coord_cartesian is a ggplot add-on function that determines axes sizes
  ## labs is a ggplot add-on function that determines the plot labels; x-axis, y-axis, title, legend colors
  ## scale_color_manual is a ggplot add-on function that sets the colors of the plot, corresponding to the groups determined in the the base function ("col" in aes)
  ## Theme is used to manipulate the physical location of the plot title, and the appearance of the entire plot

volyhippo = ggplot(DESeq2_Hippo_YWvYM3) +
  geom_point(data = subset(DESeq2_Hippo_YWvYM3, `g.o.i.` == "Non-DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.2) +
  #geom_point(data = subset(DESeq2_e133, `g.o.i.` == "Non-DEGs in Neuron Differentiation"),
  #           aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.1) +
  geom_point(data = subset(DESeq2_Hippo_YWvYM3, `g.o.i.` == "DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.1) +
  geom_point(data = subset(DESeq2_Hippo_YWvYM3, `g.o.i.` == "Genes in\nPostsynaptic\nMembrane\nOrganization GO"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.8) +
  coord_cartesian(xlim = c(-1,1), ylim = c(0,20)) +
  labs(x = "log2(FC)", y = "-log10(Adj.P-val)", title = "5mo Tmem184b GT Hippocampus Gene Exp. Changes", col = "Gene Data Type") +
  geom_hline(yintercept = 1.3, type = "solid", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("darkgoldenrod4", "gray48", "navy"))

  ## geom_text_repel labels the selected points in the color commanded
  ## Re-save the variable with the added functionality; for unlabeled plot
volyhippo = volyhippo  + geom_text_repel(data = DESeq2_Hippo_YWvYM3, x = DESeq2_Hippo_YWvYM3$log2.FC., y = DESeq2_Hippo_YWvYM3$log10Pval, color = "black", aes(label = DESeq2_Hippo_YWvYM3$labs2)) 

## Plot it
volyhippo


##### 5mo hippoc MA plot #####
## Transform the BaseMean (Mean of normalized counts across all samples) column
DESeq2_Hippo_YWvYM3$BaseMean = log(DESeq2_Hippo_YWvYM3$BaseMean)

## Layer the different classes of genes, containing the entire dataset. Subset the dataset for each layer
MAyhippo = ggplot(DESeq2_Hippo_YWvYM3) +
  geom_point(data = subset(DESeq2_Hippo_YWvYM3, `g.o.i.` == "Non-DEGs"),
             aes(x = `BaseMean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.2) +
  #geom_point(data = subset(DESeq2_e133, `g.o.i.` == "Non-DEGs in Neuron Differentiation"),
  #           aes(x = `BaseMean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.2) +
  geom_point(data = subset(DESeq2_Hippo_YWvYM3, `g.o.i.` == "DEGs"),
             aes(x = `BaseMean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.2) +
  geom_point(data = subset(DESeq2_Hippo_YWvYM3, `g.o.i.` == "Genes in\nPostsynaptic\nMembrane\nOrganization GO"),
             aes(x = `BaseMean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.5) +
  coord_cartesian(xlim = c(0,13), ylim = c(-5,1.5)) +
  labs(x = "ln(Mean of Normalized TPM)", y = "log2(FC)", title = "5mo Tmem184b GT Hippocampus Gene Exp. Changes", col = "Gene Data Type") +
  geom_hline(yintercept = 0, type = "solid", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("darkgoldenrod4", "navy", "gray48"))

  ## Add labels
MAyhippo = MAyhippo + geom_text_repel(data = DESeq2_Hippo_YWvYM3, x = DESeq2_Hippo_YWvYM3$BaseMean, y = DESeq2_Hippo_YWvYM3$log2.FC., color = "black", aes(label = DESeq2_Hippo_YWvYM3$labs2))
  ## Plot it
MAyhippo

##### 5mo hippo Bar Plots of GO Terms Using ggplot2 #####
  ## Sort by Mod Comb Score and filter some of the Processes out, including the top 19 Processes
yhippo_GO_Bar = Hippo_YWvYM_GOs[order(Hippo_YWvYM_GOs$Modified_Combined_Score, decreasing = TRUE),]
  ## Re-set the row indeces
row.names(yhippo_GO_Bar) = NULL
  ## Tranform Q
yhippo_GO_Bar$neglog10Q = -log10(yhippo_GO_Bar$GO_Biological_Process_2018.Adjusted.P.value)
  ## Filter
yhippo_GO_Bar = yhippo_GO_Bar[1:10,]
  ## Remove GO terms
yhippo_GO_Bar$GO_Biological_Process_2018.Term = yhippo_GO_Bar$GO_Biological_Process_2018.Term %>% gsub(x = yhippo_GO_Bar$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")

  ## Condense the wordy Process strings
yhippo_GO_Bar$GO_Biological_Process_2018.Term[4] = "regulation of\nneurotransmitter\nreceptor activity"
yhippo_GO_Bar$GO_Biological_Process_2018.Term[5] = "regulation of\nchaperone-mediated autophagy"
yhippo_GO_Bar$GO_Biological_Process_2018.Term[6] = "mRNA cis splicing,\nvia spliceosome"
yhippo_GO_Bar$GO_Biological_Process_2018.Term[7] = "positive regulation of\nsynaptic transmission,\nglutamatergic"
yhippo_GO_Bar$GO_Biological_Process_2018.Term[8] = "regulation of\noxidative phosphorylation"
yhippo_GO_Bar$GO_Biological_Process_2018.Term[9] = "positive regulation of\nexcitatory postsynaptic potential"


  ## Re-scale the dataframe if using terms along x-axis and Q vals on y-axis
#yhippo_GO_Bar = yhippo_GO_Bar[order(yhippo_GO_Bar$GO_Biological_Process_2018.Adjusted.P.value, decreasing = TRUE),]

  ## Plot the bar plot
  ## Store as a variable
  ## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
## scale_fill_continuous() creates the color gradient and legend
## coord_flip() switches axes for better visualization of the long GO terms


yHippo_GO_BAR = ggplot(yhippo_GO_Bar, aes(x = reorder(`GO_Biological_Process_2018.Term` , `Modified_Combined_Score`), `Modified_Combined_Score`)) + 
  geom_col(stat = "identity" , aes(fill = `Modified_Combined_Score`)) + 
  scale_fill_continuous(name = "Modified\nCombined\nScore") + 
  labs(x = "2018 GO Biological Process", y = "Modified Enrichr Combined Score", title = "5mo Tmem184b GT Hippocampus GO Analysis") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5)) +
  coord_flip()


## Plot it.
yHippo_GO_BAR



##### 5mo hippo Bar Plots of Panther Pathways Using ggplot2 #####
  ## Sort by Modified Combined Score (if not already done), add transformed Q, and filter by the top 10
yhippo_Pathway_Bar = Hippo_YWvYM_PATHWAYS[order(Hippo_YWvYM_PATHWAYS$Modified_Combined_Score, decreasing = TRUE),]

yhippo_Pathway_Bar$neglogQ = -log10(yhippo_Pathway_Bar$Panther_2016.Adjusted.P.value)

yhippo_Pathway_Bar = yhippo_Pathway_Bar[1:10,]
  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
yhippo_Pathway_Bar$Panther_2016.Term = yhippo_Pathway_Bar$Panther_2016.Term %>% gsub(x = yhippo_Pathway_Bar$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")

yhippo_Pathway_Bar$Panther_2016.Term[2] = "Opioid\nproopiomelanocortin\npathway"
yhippo_Pathway_Bar$Panther_2016.Term[3] = "Opioid\nproenkaphalin\npathway"
yhippo_Pathway_Bar$Panther_2016.Term[4] = "Nicotine\npharmacodynamics\npathway"
yhippo_Pathway_Bar$Panther_2016.Term[6] = "Alzheimer disease-\namyloid secretase pathway"
yhippo_Pathway_Bar$Panther_2016.Term[7] = "Alpha adrenergic\nreceptor signaling pathway"
yhippo_Pathway_Bar$Panther_2016.Term[8] = "Formyltetrahydroformate\nbiosynthesis"
yhippo_Pathway_Bar$Panther_2016.Term[10] = "Dopamine receptor\nmediated signaling pathway"


  ## Plot the bar plot
  ## Store as a variable
  ## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long pathway names
yhippo_PATHWAY_BAR_Q = ggplot(yhippo_Pathway_Bar, aes(x = reorder(`Panther_2016.Term`, `Modified_Combined_Score`), y = `Modified_Combined_Score`)) + 
  geom_col(stat = "identity" , aes(fill = `Modified_Combined_Score`)) + 
  scale_fill_continuous(name = "Modified\nCombined\nScore") + 
  labs(x = "Panther 2016 Pathways", y = "Modified Combined Score", title = "5mo Tmem184b GT Hippocampus Pathway Analysis") + 
  theme_light() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  coord_flip()

  ## Plot it.
yhippo_PATHWAY_BAR_Q



##### 5mo hippo Heatmaps with labels using ggplot #####

  ## Import the relevant data
yhippo_DEGs_and_GOs =  read_csv("M:/Erik/Data/Omics/TimeCourse/Custom Python Enrichr GO Clustergram e13 GT.csv")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
yhippo_GO_Heatmap = melt(yhippo_DEGs_and_GOs, id = "DEGs")
colnames(yhippo_GO_Heatmap) = c("DEGs", "GOs", "Presence In Process")

  ## Remove the "GO...." strings in the GO Process column
yhippo_GO_Heatmap$GOs = yhippo_GO_Heatmap$GOs %>% gsub(x = yhippo_GO_Heatmap$GOs, pattern = " \\(.*\\)$", replacement = "")

yhippo_GO_Heatmap$GOs[343:513] = "regulation of\ncalcium ion-dependent exocytosis"
yhippo_GO_Heatmap$GOs[685:855] = "establishment of\nchromosome localization"
yhippo_GO_Heatmap$GOs[1198:1368] = "G1/S transition of\nmitotic cell cycle"


  ## Make the digital values factors so it's compatible to graph
yhippo_GO_Heatmap$`Presence In Process` = as.factor(yhippo_GO_Heatmap$`Presence In Process`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

yhippo_GO_HEATlabs = ggplot(yhippo_GO_Heatmap, aes(`GOs`, `DEGs`)) + geom_tile(aes(fill = `Presence In Process`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "2018 GO Biological Processes", y = "Differentially Expressed Genes (FDR < 0.01)", title = "Tmem184b-GT-affected GO Processes at e13") + guides(fill = guide_legend(title = "Presence \nin Process")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


  ## Plot the heatmap
yhippo_GO_HEATlabs

  ## Import the relevant data
yhippo_DEGs_and_Pathways =  read_csv("M:/Erik/Data/Omics/TimeCourse/Custom Python Enrichr Pathway Clustergram e13 GT.csv")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
yhippo_Pathways_Heatmap = melt(yhippo_DEGs_and_Pathways, id = "DEGs")
colnames(yhippo_Pathways_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
yhippo_Pathways_Heatmap$Pathways = yhippo_Pathways_Heatmap$Pathways %>% gsub(x = yhippo_Pathways_Heatmap$Pathways, pattern = " Homo sapiens .+.?$", replacement = "")

  ## Make the digital values factors so it's compatible to graph
yhippo_Pathways_Heatmap$`Presence In Pathway` = as.factor(yhippo_Pathways_Heatmap$`Presence In Pathway`)
yhippo_Pathways_Heatmap$Pathways = as.factor(yhippo_Pathways_Heatmap$Pathways)

yhippo_Pathways_Heatmap$Pathways[757:840] = "Muscarinic acetylcholine receptor\n2 and 4 signaling pathway"
yhippo_Pathways_Heatmap$Pathways[1:84] = "Axon guidance\nmediated by Slit/Robo"
yhippo_Pathways_Heatmap$Pathways[85:168] = "De novo pyrimidine\ndeoxyribonucleotide biosynthesis"
yhippo_Pathways_Heatmap$Pathways[589:672] = "Dopamine receptor-meidated\nsignaling pathway"

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

yhippo_Pathway_HEATlabs = ggplot(yhippo_Pathways_Heatmap, aes(x = `Pathways`, y = `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Pathway`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "Panther 2016 Processes", y = "Differentially Expressed Genes (FDR < 0.01)", title = "Tmem184b-GT-affected Panther Pathways at e13") + guides(fill = guide_legend(title = "Presence\nin Pathway")) + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


  ## Plot the heatmap
yhippo_Pathway_HEATlabs


##### Create a pathway/GO term-specific heatmap across replicates #####
  ## Import the e13 normalized counts file. These are expression estimates for each gene, for each sample/replicate, where each gene's value is scaled by its sample's effect size
yhippo_Normalized_Counts = read_csv("M:/Erik/Data/Omics/Hippocampus/Genotype Comparisons/Processed Galaxy Output/Counts Files to Upload/Galaxy111_rLog-Normalized_counts_file_on_young_Hippocampus_wo_ymu4141.csv", col_names = TRUE)
  ## Rename columns
colnames(yhippo_Normalized_Counts) = c("GeneID", "WT1", "WT2", "WT3", "WT4", "Mut1", "Mut2", "Mut3")
  ## Make sure the Enrichr reference results are correctly ordered when subsetting the "2)" iteration.
Hippo_YWvYM_GOs = Hippo_YWvYM_GOs[order(Hippo_YWvYM_GOs$Modified_Combined_Score, decreasing = TRUE), ]
  ## Re-set the index to make sure indeces are correct
row.names(Hippo_YWvYM_GOs) = NULL

  ## Subset the dataframe (and run the subsequent code all at once) for one of the 3 following parameters:
  ## By the results obtained in the Enrichr search, identifying "neuron differentiation" as a significantly affected biological process:
  ## 1) For all genes in the term
#RNASeqRepResultsNeurDiff = e13_Normalized_Counts[e13_Normalized_Counts$GeneID %in% c(neuron_differentiation),]
  ## 2) For all the genes in the term + Neurog1, Neurog2, Tlx3
RNASeqRepResultsPMO = yhippo_Normalized_Counts[yhippo_Normalized_Counts$GeneID %in% c(Postsynaptic_Membrane_Organization),]
  ## 3) For only the overlapping (DE) genes
#RNASeqRepResultsNeurDiff = e13_Normalized_Counts[e13_Normalized_Counts$GeneID %in% c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]]),]
  ## 4) For only DE downregulated overlapping genes
#RNASeqRepResultsNeurDiff = e13_Normalized_Counts[e13_Normalized_Counts$GeneID %in% c(DESeq2_e133$GeneID[Neuron_Diff_Overlap[Down_Genes_indeces_in_Neuron_Diff]]),]

  ## 5) For all DEGs in the DESeq2 dataset
#RNASeqRepResultsNeurDiff = e13_Normalized_Counts[e13_Normalized_Counts$GeneID %in% c(DESeq2_e133$GeneID), ]


  ## Create dataframes to calculate Z-scores across replicates of each gene.
  ## These will hold the original data frame values until filled in later commands
  ## "Big" will store the means and sds of each gene
RNASeqRepResultsPMOBIG = RNASeqRepResultsPMO
## "3" will store only the means by genes for all replicates
RNASeqRepResultsPMO3 = RNASeqRepResultsPMO
## "4" will store only the sds by gene for all replicates
RNASeqRepResultsPMO4 = RNASeqRepResultsPMO
## "Z" will store only the Z-scores, which will be used directly to create the heatmaps
RNASeqRepResultsPMOZ = RNASeqRepResultsPMO

## Create a new column to fill with the means and standard deviations of each gene's transcript counts per million ("TPM")
RNASeqRepResultsPMOBIG$mean = 0
RNASeqRepResultsPMOBIG$sd = 0

## Loop through the data frame and fill the "mean" and "sd" columns with their appropriate values
for (i in 1:nrow(RNASeqRepResultsPMOBIG)){
  RNASeqRepResultsPMOBIG$mean[i] = (sum(RNASeqRepResultsPMO[i,c(2:8)])/ncol(RNASeqRepResultsPMO[,c(2:8)]))
}
for (j in 1:nrow(RNASeqRepResultsPMOBIG)){
  RNASeqRepResultsPMOBIG$sd[j] = sd(RNASeqRepResultsPMO[j,c(2:8)])
}
## Create a dataframe, storing the gene-specific mean of normalized TPM in all columns/replicates for Z-score calculating
RNASeqRepResultsPMO3[,c(2:8)] = RNASeqRepResultsPMOBIG$mean
## Create a dataframe, storing the gene-specific normalized TPM standard deviations in all columns/replicates for Z-score calculating
RNASeqRepResultsPMO4[,c(2:8)] = RNASeqRepResultsPMOBIG$sd
## Create the Z-score dataframe
RNASeqRepResultsPMOZ[,c(2:8)] = (RNASeqRepResultsPMO[,c(2:8)] - RNASeqRepResultsPMO3[,c(2:8)])/RNASeqRepResultsPMO4[,c(2:8)]
## Remove the genes that were detected to have 0 TPMs across all samples
## Create a function that evaluates a vector/dataframe's (x's) numerical values. Returns equal length vector with T/F bools.
## Subset the dataframe that filters those genes.
row_has_na = apply(RNASeqRepResultsPMOZ, 1, function(x){any(is.na(x))})
RNASeqRepResultsPMOZ = RNASeqRepResultsPMOZ[!row_has_na,]

## Carefully subset the list of interest; this changes depending on what is desired to be shown
#RNASeqRepResultsPMOZA = RNASeqRepResultsPMOZ %>%
#  filter(GeneID %in% DESeq2_Hippo_YWvYM3$GeneID[c(1:1635,2087)])

RNASeqRepResultsPMOZA = RNASeqRepResultsPMOZ

## Create a list of gene names **ordered by Euclidean distance** by which to re-order the dataframe/heatmap
Euclid_dist_order = hclust(dist(RNASeqRepResultsPMOZA[,c(2:8)], method = "euclidean"))$order
## The names (not the numbers)
Euclid_dist_ord_Genes = c(RNASeqRepResultsPMOZA$GeneID[Euclid_dist_order])

## Transform again to order by clusters (Euclidean distances)
#RNASeqRepResultsPMOZA = RNASeqRepResultsPMOZA %>%
# mutate(GeneID =  factor(GeneID, levels = Euclid_dist_ord_Genes)) %>%
#  arrange(GeneID)

## Re-arrange the data so that columns and rows can be run appropriately in a heatmap; give appropriate column names to this new dataframe
RNASeqRepResultsPMO2 = melt(RNASeqRepResultsPMOZ, id = "GeneID")
colnames(RNASeqRepResultsPMO2) = c("GeneID", "Genotype", "TPM Z-score")

## Create and store the heatmap core as a variable, with appropriate scale gradient limits (look up the max/min Z-score!)
## Use geom_tile to map the transcript counts by Z-score
## Fill the scale gradient, with red = downregulation, green = upregulation; Title the legend
## Add labels

## Customize by removing the filler surrounding the graph and the tickmarks
## Center the title and adjust the genes/genotype text sizes and locations
yHippoc_HEATER = ggplot(RNASeqRepResultsPMO2, aes(`Genotype`, `GeneID`)) + 
  geom_tile(aes(fill = `TPM Z-score`), color = "black") + 
  scale_fill_gradient2(low = "navy", high = "gold3", name = "TPM\nZ-score", limits = c(-1.26,1.76)) + 
  labs(x = "Genotype", y = "Downregulated Genes within Process", title = "e13 Tmem184b GT Trx Profile of\n'neuron differentiation' (GO: 0030182)") + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(plot.title = element_text(hjust = 0.5), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), axis.text.y = element_text(vjust = 0.5, size = 8), legend.title = element_text(hjust = 0))

## Plot the heatmap
yHippoc_HEATER



## Find where the "Postsynaptic Membrane Organization" DEGs are in the clustered matrix for subsetting the heatmap
#which(Euclid_dist_ord_Genes == DESeq2_Hippo_YWvYM3$GeneID[Postsyn_memb_org_Overlap[1]]) # 8
## Check 8 in the the re-ordered list containing gene names; should be "Lrp4"
#Euclid_dist_ord_Genes[8]
#which(Euclid_dist_ord_Genes == DESeq2_Hippo_YWvYM3$GeneID[Postsyn_memb_org_Overlap[2]]) # 12
#which(Euclid_dist_ord_Genes == DESeq2_Hippo_YWvYM3$GeneID[Postsyn_memb_org_Overlap[3]]) # 11
#which(Euclid_dist_ord_Genes == DESeq2_Hippo_YWvYM3$GeneID[Postsyn_memb_org_Overlap[4]]) # 14
#which(Euclid_dist_ord_Genes == DESeq2_Hippo_YWvYM3$GeneID[Postsyn_memb_org_Overlap[5]]) # 10
#which(Euclid_dist_ord_Genes == DESeq2_Hippo_YWvYM3$GeneID[Postsyn_memb_org_Overlap[6]]) # 9
#which(Euclid_dist_ord_Genes == DESeq2_Hippo_YWvYM3$GeneID[Postsyn_memb_org_Overlap[7]]) # 13


## Concatenate the list of indeces that refer to genes in the overlap list (Postsynaptic_Membrane_Organization) within the re-ordered-by-clustered-Euclidean-distance list of genes. 
## This should read out the list of genes that are DE within the Postsynaptic Membrane Organization GO Bio Process
Euclid_dist_ord_Genes[c(8,12,11,14,10,9,13)]
## Concatenate the list of indeces that refer to the indeces of the Zscore dataframe (RNASeqRepResultsNeurDiffZA) where those overlap genes are.
## This should read out the indeces within the Zscore dataframe of each gene in the overlap list.
Euclid_dist_order[c(8,12,11,14,10,9,13)]
## Check by subsetting the Zscore dataframe
RNASeqRepResultsPMOZA[Euclid_dist_order[c(8,12,11,14,10,9,13)],]

## Create the pheatmap for the DEGs involved in "neuron differentiation"
pheatmap(mat = RNASeqRepResultsPMOZA[Euclid_dist_order[c(8,12,11,14,10,9,13)] ,2:8], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, labels_row = c(Euclid_dist_ord_Genes[c(8,12,11,14,10,9,13)]))



##### Intersect the DIV14 Down DEGs with the Rescue Up DEGs #####

## Subset the DEGs in mutant embryonic cultures at DIV14
MUT_IV14 = c(which(DESeq2_In_vitro_DIV143$AdjP < 0.01))
MUT_IV14 = DESeq2_In_vitro_DIV143[MUT_IV14,]
## Save the gene names as a list
Mut_DIV14 = c(MUT_IV14$GeneID)
## Subset DEGs down in this subset
DOWN_MUT_IV14 = c(which(DESeq2_In_vitro_DIV143$AdjP < 0.01 & DESeq2_In_vitro_DIV143$log2.FC. < 0))
DOWN_MUT_IV14 = DESeq2_In_vitro_DIV143[DOWN_MUT_IV14,]
## Save the gene names as a list
Down_Mut_DIV14 = c(DOWN_MUT_IV14$GeneID)

## Subset the DEGs in rescues
RESCUE = c(which(DESeq2_In_vitro_Rescue3$AdjP < 0.01))
RESCUE = DESeq2_In_vitro_Rescue3[RESCUE,]
## Save the gene names as a list
Rescue = c(RESCUE$GeneID)
## Subset DEGs up in this subset
UP_RESCUE = c(which(DESeq2_In_vitro_Rescue3$AdjP < 0.01 & DESeq2_In_vitro_Rescue3$log2.FC. > 0))
UP_RESCUE = DESeq2_In_vitro_Rescue3[UP_RESCUE,]
## Save the gene names as a list
Up_Rescue = c(UP_RESCUE$GeneID)

## Create a list of the genes down in mutants and up in rescues
## Combine the  DEGs
Tmem_Direct_Targets = c(Down_Mut_DIV14, Up_Rescue)
## Find how many genes are significantly different along with Tmem (both directions)
length(which(duplicated(Tmem_Direct_Targets)))
## Find which genes
which(duplicated(Tmem_Direct_Targets))

## Create an index of those genes
Tmem_Direct_Target_Candidates_Index = c(which(duplicated(Tmem_Direct_Targets)))
## Confirm there are as many as those who were duplicated in the concatenated list
length(Tmem_Direct_Targets[Tmem_Direct_Target_Candidates_Index])
## Create a list of only those genes by name
Tmem_Direct_Target_Candidates = c(Tmem_Direct_Targets[Tmem_Direct_Target_Candidates_Index])

## Subset the DEGs in mutant embryonic cultures at DIV14
MUT_IV14 = c(which(DESeq2_In_vitro_DIV143$AdjP < 0.01))
MUT_IV14 = DESeq2_In_vitro_DIV143[MUT_IV14,]
## Save the gene names as a list
Mut_DIV14 = c(MUT_IV14$GeneID)
## Subset DEGs up in this subset
UP_MUT_IV14 = c(which(DESeq2_In_vitro_DIV143$AdjP < 0.01 & DESeq2_In_vitro_DIV143$log2.FC. > 0))
UP_MUT_IV14 = DESeq2_In_vitro_DIV143[UP_MUT_IV14,]
## Save the gene names as a list
Up_Mut_DIV14 = c(UP_MUT_IV14$GeneID)

## Subset the DEGs in rescues
RESCUE = c(which(DESeq2_In_vitro_Rescue3$AdjP < 0.01))
RESCUE = DESeq2_In_vitro_Rescue3[RESCUE,]
## Save the gene names as a list
Rescue = c(RESCUE$GeneID)
## Subset DEGs up in this subset
DOWN_RESCUE = c(which(DESeq2_In_vitro_Rescue3$AdjP < 0.01 & DESeq2_In_vitro_Rescue3$log2.FC. < 0))
DOWN_RESCUE = DESeq2_In_vitro_Rescue3[DOWN_RESCUE,]
## Save the gene names as a list
Down_Rescue = c(DOWN_RESCUE$GeneID)

## Create a list of the genes up in mutants and down in rescues
## Combine the  DEGs
Tmem_Indirect_Targets = c(Up_Mut_DIV14, Down_Rescue)
## Find how many genes are significantly different along with Tmem (both directions)
length(which(duplicated(Tmem_Indirect_Targets)))
## Find which genes
which(duplicated(Tmem_Indirect_Targets))

## Create an index of those genes
Tmem_Indirect_Target_Candidates_Index = c(which(duplicated(Tmem_Indirect_Targets)))
## Confirm there are as many as those who were duplicated in the concatenated list
length(Tmem_Indirect_Targets[Tmem_Indirect_Target_Candidates_Index])
## Create a list of only those genes by name
Tmem_Indirect_Target_Candidates = c(Tmem_Indirect_Targets[Tmem_Indirect_Target_Candidates_Index])




IX_GOs = as.data.frame(
  enrichr(
    c(Tmem_Direct_Target_Candidates), DBs[1])
)

## second: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
IX_PATHWAYS = as.data.frame(
  enrichr(
    c(Tmem_Direct_Target_Candidates), DBs[2])
)

IX_GOs = IX_GOs %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
## Convert the values into integers for computation
IX_GOs$Overlap = as.numeric(IX_GOs$Overlap)
IX_GOs$`Process Size` = as.numeric(IX_GOs$`Process Size`)

## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
IX_GOs = add_column(IX_GOs, EnrichrZscore = IX_GOs$GO_Biological_Process_2018.Combined.Score/log(IX_GOs$GO_Biological_Process_2018.P.value), .before = 6)

## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
IX_GOs = add_column(IX_GOs, Weighted_Overlap_Ratio = IX_GOs$Overlap*(IX_GOs$Overlap/IX_GOs$`Process Size`), .before = 4)
## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
IX_GOs = add_column(IX_GOs, Modified_Combined_Score = (abs(IX_GOs$EnrichrZscore))*(-log(IX_GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)



IX_PATHWAYS = IX_PATHWAYS %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
## Convert the values into integers for computation
IX_PATHWAYS$Overlap = as.numeric(IX_PATHWAYS$Overlap)
IX_PATHWAYS$`Process Size` = as.numeric(IX_PATHWAYS$`Process Size`)
## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
IX_PATHWAYS = add_column(IX_PATHWAYS, EnrichrZscore = IX_PATHWAYS$Panther_2016.Combined.Score/log(IX_PATHWAYS$Panther_2016.P.value), .before = 6)
## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
IX_PATHWAYS = add_column(IX_PATHWAYS, Weighted_Overlap_Ratio = IX_PATHWAYS$Overlap*(IX_PATHWAYS$Overlap/IX_PATHWAYS$`Process Size`), .before = 4)
## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
IX_PATHWAYS = add_column(IX_PATHWAYS, Modified_Combined_Score = (abs(IX_PATHWAYS$EnrichrZscore))*(-log(IX_PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)



## Clean up the dataframes
IX_GOs = IX_GOs[,-c(9,10,11)]
IX_PATHWAYS = IX_PATHWAYS[,-c(9,10,11)]


##### Intersect the DIV14 Up DEGs with the Rescue Down DEGs #####
## Repeat for supplemental material (i.e. DEGs opposite to TMEM)
IX_supp_GOs = as.data.frame(
  enrichr(
    c(Tmem_Indirect_Target_Candidates), DBs[1])
)
## second: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
IX_supp_PATHWAYS = as.data.frame(
  enrichr(
    c(Tmem_Indirect_Target_Candidates), DBs[2])
)

IX_supp_GOs = IX_supp_GOs %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
## Convert the values into integers for computation
IX_supp_GOs$Overlap = as.numeric(IX_supp_GOs$Overlap)
IX_supp_GOs$`Process Size` = as.numeric(IX_supp_GOs$`Process Size`)

## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
IX_supp_GOs = add_column(IX_supp_GOs, EnrichrZscore = IX_supp_GOs$GO_Biological_Process_2018.Combined.Score/log(IX_supp_GOs$GO_Biological_Process_2018.P.value), .before = 6)

## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
IX_supp_GOs = add_column(IX_supp_GOs, Weighted_Overlap_Ratio = IX_supp_GOs$Overlap*(IX_supp_GOs$Overlap/IX_supp_GOs$`Process Size`), .before = 4)
## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
IX_supp_GOs = add_column(IX_supp_GOs, Modified_Combined_Score = (abs(IX_supp_GOs$EnrichrZscore))*(-log(IX_supp_GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)



IX_supp_PATHWAYS = IX_supp_PATHWAYS %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
## Convert the values into integers for computation
IX_supp_PATHWAYS$Overlap = as.numeric(IX_supp_PATHWAYS$Overlap)
IX_supp_PATHWAYS$`Process Size` = as.numeric(IX_supp_PATHWAYS$`Process Size`)
## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
IX_supp_PATHWAYS = add_column(IX_supp_PATHWAYS, EnrichrZscore = IX_supp_PATHWAYS$Panther_2016.Combined.Score/log(IX_supp_PATHWAYS$Panther_2016.P.value), .before = 6)
## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
IX_supp_PATHWAYS = add_column(IX_supp_PATHWAYS, Weighted_Overlap_Ratio = IX_supp_PATHWAYS$Overlap*(IX_supp_PATHWAYS$Overlap/IX_supp_PATHWAYS$`Process Size`), .before = 4)
## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
IX_supp_PATHWAYS = add_column(IX_supp_PATHWAYS, Modified_Combined_Score = (abs(IX_supp_PATHWAYS$EnrichrZscore))*(-log(IX_supp_PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)



## Clean up the dataframes
IX_supp_GOs = IX_supp_GOs[,-c(9,10,11)]
IX_supp_PATHWAYS = IX_supp_PATHWAYS[,-c(9,10,11)]



##### Intersection Rescue Bar Plots of GO Processes and Pathways Using ggplot2 #####

## Sort by -log10(Q). Retain the top 10 Processes 
IX_GOs$neglog10Q = -log10(IX_GOs$GO_Biological_Process_2018.Adjusted.P.value)
IX_GOs = IX_GOs[order(IX_GOs$neglog10Q, decreasing = TRUE),]
IX_GOs_Bar_ADJP = IX_GOs[1:10,]
## Remove the "GO...." strings in the GO Process column
IX_GOs_Bar_ADJP$GO_Biological_Process_2018.Term = IX_GOs_Bar_ADJP$GO_Biological_Process_2018.Term %>% gsub(x = IX_GOs_Bar_ADJP$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")

## Shorten the term strings
IX_GOs_Bar_ADJP$GO_Biological_Process_2018.Term[4] = "positive regulation of transcription\nfrom RNA polymerase II promoter\nin response to endoplasmic reticulum stress"
IX_GOs_Bar_ADJP$GO_Biological_Process_2018.Term[9] = "regulation of aspartic-type endopeptidase activity\ninvolved in amyloid precursor protein catabolic process"
IX_GOs_Bar_ADJP$GO_Biological_Process_2018.Term[10] = "intrinsic apoptotic signaling pathway in response to\nendoplasmic reticulum stress"


## Plot the bar plot
## Store as a variable
## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
## scale_fill_continuous() creates the color gradient and legend
## coord_flip() switches axes for better visualization of the long GO terms
IX_GOs_BAR_ADJP = ggplot(IX_GOs_Bar_ADJP, aes(x = reorder(`GO_Biological_Process_2018.Term`, `GO_Biological_Process_2018.Adjusted.P.value`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "2018 GO Biological Process", y = "Enrichr -log10(Adj. P-val)", title = "GO Analysis of DEGs Affected by Tmem184b in vitro") +
  theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Plot it.
IX_GOs_BAR_ADJP



## Sort by -log10(Q). Retain the top 10 Processes 
IX_PATHWAYS$neglog10Q = -log10(IX_PATHWAYS$Panther_2016.Adjusted.P.value)
IX_PATHWAYS = IX_PATHWAYS[order(IX_PATHWAYS$neglog10Q, decreasing = TRUE),]
IX_Pathway_Bar_ADJP = IX_PATHWAYS[1:10,]
## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
IX_Pathway_Bar_ADJP$Panther_2016.Term = IX_Pathway_Bar_ADJP$Panther_2016.Term %>% gsub(x = IX_Pathway_Bar_ADJP$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")

## Shorten the term strings
IX_Pathway_Bar_ADJP$Panther_2016.Term[2] = "Heterotrimeric G-protein signaling pathway-\nGq alpha and Go alpha mediated pathway"

## Plot the bar plot
## Store as a variable
## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
## scale_fill_continuous() creates the color gradient and legend
## coord_flip() switches axes for better visualization of the long pathway names
IX_PATHWAY_BAR_ADJP = ggplot(IX_Pathway_Bar_ADJP, aes(x = reorder(`Panther_2016.Term`, `Panther_2016.Adjusted.P.value`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "Panther 2016 Pathways", y = "Enrichr -log10(Adj. P-val)", title = "Pathway Analysis of DEGs Affected by Tmem184b in vitro") +
  theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Plot it.
IX_PATHWAY_BAR_ADJP




## Export data as .CSVs to the server
#write.table(IX_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\IX DIV14 GT DRG GO Bio Processes 0 LFC.csv", sep = ",", col.names = T)
#write.table(IX_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\IX DIV14 GT DRG Panther Pathways 0 LFC.csv", sep = ",", col.names = T)


##### Intersection Rescue GO Bio Processes and Panther Pathways Heatmaps with labels using ggplot2 #####
## Import the relevant data
IX_Pathways =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr Pathway Clustergram IX.csv")
## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
IX_Pathway_Heatmap = melt(IX_Pathways, id = "DEGs")
colnames(IX_Pathway_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

## Remove excessive pathway name strings
IX_Pathway_Heatmap$Pathways = IX_Pathway_Heatmap$Pathways %>% gsub(x = IX_Pathway_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

## Make the digital values factors so it's compatible to graph
IX_Pathway_Heatmap$`Presence In Pathway` = as.factor(IX_Pathway_Heatmap$`Presence In Pathway`)

## Create and store the heatmap core as a variable
## Use geom_tile to map the binary values
## Customize by removing the filler surrounding the graph and the tickmarks
## Center the Title

IX_HEATlabs = ggplot(IX_Pathway_Heatmap, aes(`Pathways`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Pathway`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.01)", title = "Pathway Analysis of DEGs Affected by Tmem184b in vitro") + 
  guides(fill = guide_legend(title = "Presence in Pathway")) + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()

## Plot the heatmap
IX_HEATlabs


## Import the relevant data
IX_GOs =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr GO Clustergram IX.csv")
## Export the gene list for photoshopping
#write.csv(DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\aDRG FDR 05 Heatmap DEGs.csv", sep = ",")


## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
IX_GO_Heatmap = melt(IX_GOs, id = "DEGs")
colnames(IX_GO_Heatmap) = c("DEGs", "GOs", "Presence In Process")

## Remove the "GO...." strings in the GO Process column
IX_GO_Heatmap$GOs = IX_GO_Heatmap$GOs %>% gsub(x = IX_GO_Heatmap$GOs, pattern = " \\(.*\\)$", replacement = "")

## Make the digital values factors so it's compatible to graph
IX_GO_Heatmap$`Presence In Process` = as.factor(IX_GO_Heatmap$`Presence In Process`)

## Create and store the heatmap core as a variable
## Use geom_tile to map the binary values
## Customize by removing the filler surrounding the graph and the tickmarks
## Center the Title

IX_HEATlabs = ggplot(IX_GO_Heatmap, aes(`GOs`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Process`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "2018 GO Biological Process", y = "Differentially Expressed Genes (FDR < 0.05)", title = "GO Analysis of DEGs Affected by Tmem184b in vitro") + 
  guides(fill = guide_legend(title = "Presence in Process")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


## Plot the heatmap
IX_HEATlabs



