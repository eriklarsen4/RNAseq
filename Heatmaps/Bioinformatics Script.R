

## Script for obtaining Enrichr-based bioinformatics data
## Developed by Erik Larsen

  ## Condense lines of code by command: Alt+o.

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
## (Don't update if the option pops up
## "n")

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
#library(dplyr)
library(readr)
#library(readxl)
#library(tibble)
library(AnnotationDbi)
library(stats4)
library(BiocGenerics)
library(parallel)
#library(ggsci)
#library(viridis)
library(dendextend)
library(ggdendro)
library(grid)
#library(EnhancedVolcano)
library(pheatmap)
#library(RColorBrewer)
  ## Create Martha's heatmap palette for downstream use
  #Marthas_palette = colorRampPalette(c("blue", "black", "yellow"))(n = 1000)
    # ("darkgoldenrod4", "white", and "navy" with n = 50 is preferable)

  ## Connect live to the Enrichr server/master database (website) and store the "space" as a variable. This contains a vast array of bioinformatics databases/websites
    ## Store it as a variable for better visualization and eventual subsetting
DBs = listEnrichrDbs()

  ## View the variable in the editor to find the relevant databases and their indeces for subsetting (click on the dataframe in the global environment and manually peruse)

  ## In this case, "GO_Biological_Process_2018" (index # 130), and "Panther_2016" (index # 102)
    ## Concatenate them in a list for subsetting, or slice directly (we'll slice directly in the next line)
DBs = DBs$libraryName[c(130,102)]

##### e13 Data prep #####
  ## Import the dataset you want to analyze
DESeq2_e13 = read.csv("M:/Erik/Data/Omics/TimeCourse/Processed Galaxy Output/Test Results to Upload/DESeq2_result_file_on_GT_e13_PARA.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
DESeq2_e133 = subset(DESeq2_e13, (!is.na(DESeq2_e13[,"AdjP"])))

  ## Perform relevant analysis through Enrichr; in this case, FDR = 0.01
    ## first: GO Bio Processes; this analyzes genes from your dataset in Enrichr, which draws terms from geneontology.org
e13_GOs = as.data.frame(
  enrichr(
    c(DESeq2_e133$GeneID[which(DESeq2_e133$AdjP < 0.01)]), DBs[1])
  )
  ## second: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
e13_PATHWAYS = as.data.frame(
  enrichr(
    c(DESeq2_e133$GeneID[which(DESeq2_e133$AdjP < 0.01)]), DBs[2])
  )

  ## Clean up both dataframes and export both dataframes.
      ## Dataframes for binary heatmap visualization (presence in pathway/term or not) need additional cleaning in Python;
        ## Following that manipulation, data will be re-imported for visualization
e13_GOs = e13_GOs %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
    ## Convert the values into integers for computation
e13_GOs$Overlap = as.numeric(e13_GOs$Overlap)
e13_GOs$`Process Size` = as.numeric(e13_GOs$`Process Size`)
    ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
e13_GOs = add_column(e13_GOs, EnrichrZscore = e13_GOs$GO_Biological_Process_2018.Combined.Score/log(e13_GOs$GO_Biological_Process_2018.P.value), .before = 4)
    ## Re-order the dataframe by Q value
e13_GOs = e13_GOs[order(e13_GOs$GO_Biological_Process_2018.Adjusted.P.value, decreasing = FALSE),]
    ## Re-set the index to make sure indeces are correct
row.names(e13_GOs) = NULL


    ## Vectorize and make the DEGs within GO Biological Processes compatible with the gene names output from DESeq2 (EntrezID; capital first letter, lowercase all rest)
    ## Use for subsetting in downstream visualizations (indexing and labeling)
e13_GO_GENES = str_split(e13_GOs$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(e13_GO_GENES)){
  e13_GO_GENES[i,] = paste(substr(e13_GO_GENES[i,], 1, 1),
                    tolower(substr(e13_GO_GENES[i,], 2, 7)), sep = "") 
}

    ## Repeat same process for PATHWAYS
e13_PATHWAYS = e13_PATHWAYS %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
    ## Convert the values into integers for computation
e13_PATHWAYS$Overlap = as.numeric(e13_PATHWAYS$Overlap)
e13_PATHWAYS$`Process Size` = as.numeric(e13_PATHWAYS$`Process Size`)

    ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
e13_PATHWAYS = add_column(e13_PATHWAYS, EnrichrZscore = e13_PATHWAYS$Panther_2016.Combined.Score/log(e13_PATHWAYS$Panther_2016.P.value), .before = 4)
    ## Re-order the dataframe by Q value
e13_PATHWAYS = e13_PATHWAYS[order(e13_PATHWAYS$Panther_2016.Adjusted.P.value, decreasing = FALSE),]
    ## Re-set the row index to make sure indeces are correct
row.names(e13_PATHWAYS) = NULL

    ## Vectorize and make the DEGs within GO Biological Processes compatible with the gene names output from DESeq2 (EntrezID; capital first letter, lowercase all rest)
    ## Use for subsetting in downstream visualizations (indexing and labeling)
e13_PATHWAY_GENES = str_split(e13_PATHWAYS$Panther_2016.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(e13_PATHWAY_GENES)){
  e13_PATHWAY_GENES[i,] = paste(substr(e13_PATHWAY_GENES[i,], 1, 1),
                    tolower(substr(e13_PATHWAY_GENES[i,], 2, 7)), sep = "") 
}

    ## Remove unused columns ("old p-values" and "odds ratio")
e13_GOs = e13_GOs[,-c(7,8,9)]
e13_PATHWAYS = e13_PATHWAYS[,-c(7,8,9)]

    ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
      ##  --->>> The fraction of DEGs belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
e13_GOs = add_column(e13_GOs, Weighted_Overlap_Ratio = e13_GOs$Overlap*(e13_GOs$Overlap/e13_GOs$`Process Size`), .before = 4)
    ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
e13_GOs = add_column(e13_GOs, Modified_Combined_Score = (e13_GOs$EnrichrZscore)*(log10(e13_GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)

    ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
      ##  --->>> The fraction of DEGs belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
e13_PATHWAYS = add_column(e13_PATHWAYS, Weighted_Overlap_Ratio = e13_PATHWAYS$Overlap*(e13_PATHWAYS$Overlap/e13_PATHWAYS$`Process Size`), .before = 4)
    ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
e13_PATHWAYS = add_column(e13_PATHWAYS, Modified_Combined_Score = (e13_PATHWAYS$EnrichrZscore)*(log10(e13_PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)


    ## Find genes within a Process/Pathway of interest that are DE (down) in mutants
      ## For example, find which genes that are differentially downregulated in e13 mutants are in the "Neuron Differentiation" GO Bio Process; include the Neurogenins and Tlx3
      ## Concatenate a list of all the overlapping DEGs within the process
Neuron_Diff_Overlap = c(which(DESeq2_e133$GeneID %in% e13_GO_GENES[6,1:e13_GOs$Overlap[6]]), which(DESeq2_e133$GeneID == "Neurog1"), which(DESeq2_e133$GeneID == "Neurog2"), which(DESeq2_e133$GeneID == "Tlx3"))
      ## Concatenate a different list of the indeces from the overlapping genes that are down in mutants
#Down_Genes_indeces_in_Neuron_Diff = c(which(DESeq2_e133[c(Neuron_Diff_Overlap),3] < 0))
      ## Find which genes they are to manually confirm
#DESeq2_e133$GeneID[Neuron_Diff_Overlap[Down_Genes_indeces_in_Neuron_Diff]]


    ## Export data as .CSVs to the server
#write.table(e13_GOs, "M:\\Erik\\Data\\Omics\\Timecourse\\e13 GT DRG GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(e13_PATHWAYS, "M:\\Erik\\Data\\Omics\\Timecourse\\e13 GT DRG Panther Pathways.csv", sep = ",", col.names = T)


  ## Create a column in the DESeq2 dataframe that scales the Adjusted P-value by log-base 10
DESeq2_e133$log10Pval = -log10(DESeq2_e133$AdjP)

  ## Create the data points of interest for functional display.. ggplot struggles to plot subsets of dataframes, so creating a variable that directly labels all the points accordingly helps.
  ## Create a new column (genes of interest); Using "GeneID" is irrelevant; the values will be replaced in the next step
DESeq2_e133$g.o.i. = DESeq2_e133$GeneID

  ## Turn the concatenated list of genes, COPIED AND PASTED FROM ENRICHR, into compatible Entrez GeneIDs (lowercase all non-first letters)
    ## This reference list can be used for future plotting/labeling in volcanoes, MAs, and heatmaps.
neuron_differentiation = c("SPINK5","EPHA5","RYK","TGFB2","SMARCA1","WNT5A","WNT5B","RUNX1","SOX11","ISL2","POU4F2","NRBP2","RUNX3","DCLK2","RUNX2","SKI","NIF3L1","MAPK10","PTF1A","MTPN","RB1","LRP6","FZD10","HS6ST1","ASCL1","TUBB3","SPOCK1","PITX3","PITX2","FZD2","FZD1","WNT10B","FZD4","WNT10A","XBP1","FZD3","WNT3A","MAP2K1","FZD5","FZD8","FZD7","INHBA","MANF","FOXN4","POU3F2","CCSAP","NR4A2","WNK1","TDP2","UNC119","IL1RAPL1","PPT1","VSX1","PIN1","ITM2C","UNC119B","BRSK1","BRSK2","WNT2B","ATP8A2","FMR1","CEND1","RORB","SOX4","FGF8","ADNP2","IER2","WNT9A","WNT9B","PAX6","USH2A","TGFBR1","VEGFA","PHOX2B","WNT16","PHOX2A","PTPRD","ID1","CDK5","ID3","ID2","ID4","ALK","TENM1","BTG4","TENM2","TENM3","TENM4","WNT8B","PCSK9","WNT8A","WNT6","FXR2","FXR1","NEUROD4","WNT11","SHH","CASP3","OTX2","WNT1","SKIL","WNT2","WNT3","WNT4","NEUROD1","NTRK2","SKOR1","MYEF2","NLGN4X","SKOR2","MEF2C","EN1","THOC2","PROX1","WNT7A","WNT7B","PROX2","ATP2B2","EPOP","LMX1B","LMX1A","DDIT4","ERCC2","HAND2","ERCC3","CDK5R1","FOXA2","FOXA1","MYT1L","LDB1","CTNNB1","CDNF","SRF","ID2B","GLI2","NGRN","MAPK9","MAPK8","OPHN1")
for (i in 1:length(neuron_differentiation)){
  neuron_differentiation[i] = paste(substr(neuron_differentiation[i], 1, 1),
                                    tolower(substr(neuron_differentiation[i], 2, 8)), sep = "") 
}
  ## Include Tlx3, Neurogenins
neuron_differentiation2 = c("SPINK5","EPHA5","RYK","TGFB2","SMARCA1","WNT5A","WNT5B","RUNX1","SOX11","ISL2","POU4F2","NRBP2","RUNX3","DCLK2","RUNX2","SKI","NIF3L1","MAPK10","PTF1A","MTPN","RB1","LRP6","FZD10","HS6ST1","ASCL1","TUBB3","SPOCK1","PITX3","PITX2","FZD2","FZD1","WNT10B","FZD4","WNT10A","XBP1","FZD3","WNT3A","MAP2K1","FZD5","FZD8","FZD7","INHBA","MANF","FOXN4","POU3F2","CCSAP","NR4A2","WNK1","TDP2","UNC119","IL1RAPL1","PPT1","VSX1","PIN1","ITM2C","UNC119B","BRSK1","BRSK2","WNT2B","ATP8A2","FMR1","CEND1","RORB","SOX4","FGF8","ADNP2","IER2","WNT9A","WNT9B","PAX6","USH2A","TGFBR1","VEGFA","PHOX2B","WNT16","PHOX2A","PTPRD","ID1","CDK5","ID3","ID2","ID4","ALK","TENM1","BTG4","TENM2","TENM3","TENM4","WNT8B","PCSK9","WNT8A","WNT6","FXR2","FXR1","NEUROD4","WNT11","SHH","CASP3","OTX2","WNT1","SKIL","WNT2","WNT3","WNT4","NEUROD1","NTRK2","SKOR1","MYEF2","NLGN4X","SKOR2","MEF2C","EN1","THOC2","PROX1","WNT7A","WNT7B","PROX2","ATP2B2","EPOP","LMX1B","LMX1A","DDIT4","ERCC2","HAND2","ERCC3","CDK5R1","FOXA2","FOXA1","MYT1L","LDB1","CTNNB1","CDNF","SRF","ID2B","GLI2","NGRN","MAPK9","MAPK8","OPHN1","TLX3","NEUROG1","NEUROG2")
for (i in 1:length(neuron_differentiation2)){
  neuron_differentiation2[i] = paste(substr(neuron_differentiation2[i], 1, 1),
                                     tolower(substr(neuron_differentiation2[i], 2, 8)), sep = "") 
}

  ## Concatenate the strings of the DEGs in the e13 DESeq2 dataset that are also in this process
#c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]])
  ## Find the indeces of these genes in the DESeq2 dataframe for plotting purposes
#which(DESeq2_e133$GeneID %in% c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]]))
  ## Label these points in the volcano/MA plots

  ## Create a "not in" operator
"%ni%" = Negate("%in%")
#which(c(e13_GO_GENES[1,1:e13_GOs$Overlap[1]]) %ni% c(e13_GO_GENES[1,1:e13_GOs$`Process Size`[1]]))

  ## Indeces from the DESeq2 dataframe of genes NOT overlapping in the GO Bio Process, "Neuron Differentiation"
which(DESeq2_e133$GeneID %ni% c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]]))
  ## (the genes themselves)
DESeq2_e133$GeneID[which(neuron_differentiation %ni% c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]]))]
  ## Concatenated list of the overlapping DEGs from the DESeq2 dataset and the "Neuron Differentiation" GO
DEG_Overlap = c(1,9,11,12,17,18,23,25,27,30,37,40,43,47,54,62,67,73,82,86,95,104,105,114,119,137)
  ## Subsetted list of the DEGs of all the genes in the GO
which(DESeq2_e133$GeneID %in% neuron_differentiation[DEG_Overlap])
  ## Subsetted list of the genes in the GO not DE
neuron_differentiation[-DEG_Overlap]




##### Volcano Plot of e13 GT DRG #####

  ## Fill the points with appropriately indexed data
  ## "Non-differentially Expressed Genes" are defined as being "genes of interest" with Adjusted P-values > 0.01
DESeq2_e133$g.o.i.[which(DESeq2_e133$AdjP > 0.01)] = "Non-DEGs"
  ## "Non-differentially expressed genes" identified in the "Neuron Differentiation" GO Biological Process
#DESeq2_e133$g.o.i.[which(DESeq2_e133$GeneID %in% c(neuron_differentiation[-DEG_Overlap]))] = "Non-DEGs in Neuron Differentiation"
  ## "Differentially Expressed Genes" also identified in Neuron Differentiation GO ("overlapping" genes)
DESeq2_e133$g.o.i.[c(which(DESeq2_e133$GeneID %in% c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]])))] = "DEGs in\nNeuron Differentiation"
  ## "Differentially Expressed Genes" are defined as the genes differentially expressed (AdjP < 0.01) but not overlapping with the "Neuron Differentiation" GO genes
DESeq2_e133$g.o.i.[c(which(DESeq2_e133$AdjP <= 0.01))][-c(which(DESeq2_e133$GeneID %in% neuron_differentiation[DEG_Overlap]))] = "DEGs"

  ## Turn the plotting column into a factor vector
DESeq2_e133$g.o.i. = as.factor(DESeq2_e133$g.o.i.)
  ## Check to make sure
class(DESeq2_e133$g.o.i.)
  ## Check to make sure all the levels exist and didn't get overwritten by another/were properly indexed
levels(DESeq2_e133$g.o.i.)

  ## Create labels on points we're interested in labeling in volcano/MA plots
    ## Create the column, "labs"; hide all of the text labels with: ""
DESeq2_e133$labs = ""

## Figure out the indeces of genes we're interested in labeling
#which(DESeq2_e133$GeneID %in% c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]]))
#e13_GO_GENES[6,1:26]
#which(DESeq2_e133$GeneID == "Fzd10") # 829
#which(DESeq2_e133$GeneID == "Neurod1") # 790
#which(DESeq2_e133$GeneID == "Mapk9") # 1457
#which(DESeq2_e133$GeneID == "Neurod4") # 278
#which(DESeq2_e133$GeneID == "Wnt4") # 811
#which(DESeq2_e133$GeneID == "Fzd2") # 1431
#which(DESeq2_e133$GeneID == "Wnt3a") # 97
#which(DESeq2_e133$GeneID == "Fzd8") # 678
#which(DESeq2_e133$GeneID == "Pou4f2") # 1186
#which(DESeq2_e133$GeneID == "Mapk10") # 803

  ## Fill the appropriate rows of the "labs" column with the corresponding Gene names for labeling those points
DESeq2_e133$labs[c(which(DESeq2_e133$GeneID %in% c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]])))] = DESeq2_e133$GeneID[c(which(DESeq2_e133$GeneID %in% c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]])))]

  ## Create another "labels" column and fill with cherry-picked genes
DESeq2_e133$labs2 = ""
DESeq2_e133$labs2[c(34,97,678,790,811,829,1186,1431)] = DESeq2_e133$GeneID[c(34,97,678,790,811,829,1186,1431)]

    ## Notes on ggplot commands:
    ## geom-point = alpha is a ggplot add-on function that controls point color transparency
    ## coord_cartesian is a ggplot add-on function that determines axes sizes
    ## labs is a ggplot add-on function that determines the plot labels; x-axis, y-axis, title, legend colors
    ## scale_color_manual is a ggplot add-on function that sets the colors of the plot, corresponding to the groups determined in the the base function ("col" in aes)
    ## Theme is used to manipulate the physical location of the plot title, and the appearance of the entire plot

vole13 = ggplot(DESeq2_e133) +
  geom_point(data = subset(DESeq2_e133, `g.o.i.` == "Non-DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.2) +
  #geom_point(data = subset(DESeq2_e133, `g.o.i.` == "Non-DEGs in Neuron Differentiation"),
  #           aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.1) +
  geom_point(data = subset(DESeq2_e133, `g.o.i.` == "DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.1) +
  geom_point(data = subset(DESeq2_e133, `g.o.i.` == "DEGs in\nNeuron Differentiation"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.8) +
  coord_cartesian(xlim = c(-8,8), ylim = c(0,40)) +
  labs(x = "log2(FC)", y = "-log10(Adj.P-val)", title = "e13 Tmem184b GT DRG Gene Exp. Changes", col = "Gene Data Type") +
  geom_hline(yintercept = 2.0, type = "solid", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("darkgoldenrod4", "navy", "gray48"))

    ## geom_text_repel labels the selected points in the color commanded
      ## Re-save the variable with the added functionality; for unlabeled plot
vole13 = vole13  + geom_text_repel(data = DESeq2_e133, x = DESeq2_e133$log2.FC., y = DESeq2_e133$log10Pval, color = "black", aes(label = DESeq2_e133$labs2)) 

  ## Plot it
vole13


##### e13 MA plot #####
  ## Transform the BaseMean (Mean of normalized counts across all samples) column
DESeq2_e133$BaseMean = log(DESeq2_e133$BaseMean)

  ## Layer the different classes of genes, containing the entire dataset. Subset the dataset for each layer
MAe13 = ggplot(DESeq2_e133) +
  geom_point(data = subset(DESeq2_e133, `g.o.i.` == "Non-DEGs"),
             aes(x = `BaseMean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.2) +
  #geom_point(data = subset(DESeq2_e133, `g.o.i.` == "Non-DEGs in Neuron Differentiation"),
  #           aes(x = `BaseMean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.2) +
  geom_point(data = subset(DESeq2_e133, `g.o.i.` == "DEGs"),
             aes(x = `BaseMean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.2) +
  geom_point(data = subset(DESeq2_e133, `g.o.i.` == "DEGs in\nNeuron Differentiation"),
             aes(x = `BaseMean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.5) +
  coord_cartesian(xlim = c(0,14), ylim = c(-5,5)) +
  labs(x = "ln (Mean of TPM)", y = "log2 (FC)", title = "e13 Tmem184b GT DRG Gene Exp. Changes", col = "Gene Data Type") +
  geom_hline(yintercept = 0, type = "solid", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") + 
  scale_color_manual(values = c("darkgoldenrod4", "navy", "gray48"))

  ## Add labels
MAe13 = MAe13 + geom_text_repel(data = DESeq2_e133, x = DESeq2_e133$BaseMean, y = DESeq2_e133$log2.FC., color = "black", aes(label = DESeq2_e133$labs2))
  ## Plot it
MAe13

##### e13 Bar Plots of GO Terms Using ggplot2 #####
  ## Sort by Q and filter some of the Processes out, including the top 19 Processes
e13_GO_Bar = e13_GOs[order(e13_GOs$GO_Biological_Process_2018.Adjusted.P.value, decreasing = FALSE),]
  ## Re-set the row indeces
row.names(e13_GO_Bar) = NULL
  ## Tranform Q
e13_GO_Bar$neglog10Q = -log10(e13_GO_Bar$GO_Biological_Process_2018.Adjusted.P.value)
  ## Filter
e13_GO_Bar = e13_GO_Bar[1:10,]
  ## Remove GO terms
e13_GO_Bar$GO_Biological_Process_2018.Term = e13_GO_Bar$GO_Biological_Process_2018.Term %>% gsub(x = e13_GO_Bar$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")

  ## Condense the wordy Process strings
e13_GO_Bar$GO_Biological_Process_2018.Term[3] = "regulation of\n calcium ion-dependent exocytosis"
e13_GO_Bar$GO_Biological_Process_2018.Term[5] = "establishment of\n chromosome localization"
e13_GO_Bar$GO_Biological_Process_2018.Term[8] = "G1/S transition of\n mitotic cell cycle"

  ## Re-scale the dataframe if using terms along x-axis and Q vals on y-axis
e13_GO_Bar = e13_GO_Bar[order(e13_GO_Bar$GO_Biological_Process_2018.Adjusted.P.value, decreasing = TRUE),]

  ## Plot the bar plot
  ## Store as a variable
  ## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long GO terms


e13_GO_BAR_Q = ggplot(e13_GO_Bar, aes(x = reorder(`GO_Biological_Process_2018.Term` , `GO_Biological_Process_2018.Adjusted.P.value`), `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "2018 GO Biological Process", y = "Enrichr -log10(Adj. P-val)", title = "e13 Tmem184b GT DRG GO Analysis") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), axis.text.x = element_blank())


  ## Plot it.
e13_GO_BAR_Q



##### e13 Bar Plots of Panther Pathways Using ggplot2 #####
  ## Sort by Q and filter some of the pathways out for graphing
e13_Pathway_Bar = e13_PATHWAYS[order(e13_PATHWAYS$Panther_2016.Adjusted.P.value, decreasing = FALSE),]

e13_Pathway_Bar$Panther_2016.Adjusted.P.value = -log10(e13_Pathway_Bar$Panther_2016.Adjusted.P.value)

e13_Pathway_Bar = e13_Pathway_Bar[1:10,]
  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
e13_Pathway_Bar$Panther_2016.Term = e13_Pathway_Bar$Panther_2016.Term %>% gsub(x = e13_Pathway_Bar$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")

e13_Pathway_Bar$Panther_2016.Term[1] = "Axon guidance\nmediated by Slit/Robo"
e13_Pathway_Bar$Panther_2016.Term[2] = "De novo pyrimidine\ndeoxyribonucleotide biosynthesis"
e13_Pathway_Bar$Panther_2016.Term[8] = "Dopamine receptor-mediated\nsignaling pathway"
e13_Pathway_Bar$Panther_2016.Term[10] = "Muscarinic acetylcholine receptor\n2 and 4 signaling pathway"

  ## Plot the bar plot
  ## Store as a variable
  ## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long pathway names
e13_PATHWAY_BAR_Q = ggplot(e13_Pathway_Bar, aes(x = reorder(`Panther_2016.Term`, `Panther_2016.Adjusted.P.value`), y = `Panther_2016.Adjusted.P.value`)) + 
  geom_col(stat = "identity" , aes(fill = `Panther_2016.Adjusted.P.value`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "Panther 2016 Pathways", y = "Enrichr -log10(Adj. P-val)", title = "e13 Tmem184b GT/GT DRG Pathway Analysis") + 
  theme_light() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + 
  coord_flip()

  ## Plot it.
e13_PATHWAY_BAR_Q



##### e13 Heatmaps with labels using ggplot #####

  ## Import the relevant data
e13_DEGs_and_GOs =  read_csv("M:/Erik/Data/Omics/TimeCourse/Custom Python Enrichr GO Clustergram e13 GT.csv")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
e13_GO_Heatmap = melt(e13_DEGs_and_GOs, id = "DEGs")
colnames(e13_GO_Heatmap) = c("DEGs", "GOs", "Presence In Process")

  ## Remove the "GO...." strings in the GO Process column
e13_GO_Heatmap$GOs = e13_GO_Heatmap$GOs %>% gsub(x = e13_GO_Heatmap$GOs, pattern = " \\(.*\\)$", replacement = "")

e13_GO_Heatmap$GOs[343:513] = "regulation of\ncalcium ion-dependent exocytosis"
e13_GO_Heatmap$GOs[685:855] = "establishment of\nchromosome localization"
e13_GO_Heatmap$GOs[1198:1368] = "G1/S transition of\nmitotic cell cycle"


  ## Make the digital values factors so it's compatible to graph
e13_GO_Heatmap$`Presence In Process` = as.factor(e13_GO_Heatmap$`Presence In Process`)


  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

e13_GO_HEATlabs = ggplot(e13_GO_Heatmap, aes(`GOs`, `DEGs`)) + geom_tile(aes(fill = `Presence In Process`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "2018 GO Biological Processes", y = "Differentially Expressed Genes (FDR < 0.01)", title = "Tmem184b-GT-affected GO Processes at e13") + guides(fill = guide_legend(title = "Presence \nin Process")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


  ## Plot the heatmap
e13_GO_HEATlabs

  ## Import the relevant data
e13_DEGs_and_Pathways =  read_csv("M:/Erik/Data/Omics/TimeCourse/Custom Python Enrichr Pathway Clustergram e13 GT.csv")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
e13_Pathways_Heatmap = melt(e13_DEGs_and_Pathways, id = "DEGs")
colnames(e13_Pathways_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
e13_Pathways_Heatmap$Pathways = e13_Pathways_Heatmap$Pathways %>% gsub(x = e13_Pathways_Heatmap$Pathways, pattern = " Homo sapiens .+.?$", replacement = "")

  ## Make the digital values factors so it's compatible to graph
e13_Pathways_Heatmap$`Presence In Pathway` = as.factor(e13_Pathways_Heatmap$`Presence In Pathway`)
e13_Pathways_Heatmap$Pathways = as.factor(e13_Pathways_Heatmap$Pathways)

e13_Pathways_Heatmap$Pathways[757:840] = "Muscarinic acetylcholine receptor\n2 and 4 signaling pathway"
e13_Pathways_Heatmap$Pathways[1:84] = "Axon guidance\nmediated by Slit/Robo"
e13_Pathways_Heatmap$Pathways[85:168] = "De novo pyrimidine\ndeoxyribonucleotide biosynthesis"
e13_Pathways_Heatmap$Pathways[589:672] = "Dopamine receptor-meidated\nsignaling pathway"

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

e13_Pathway_HEATlabs = ggplot(e13_Pathways_Heatmap, aes(x = `Pathways`, y = `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Pathway`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "Panther 2016 Processes", y = "Differentially Expressed Genes (FDR < 0.01)", title = "Tmem184b-GT-affected Panther Pathways at e13") + guides(fill = guide_legend(title = "Presence\nin Pathway")) + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


## Plot the heatmap
e13_Pathway_HEATlabs


##### Create a pathway/GO term-specific heatmap across replicates #####
  ## Import the e13 normalized counts file. These are expression estimates for each gene, for each sample/replicate, where each gene's value is scaled by its sample's effect size
e13_Normalized_Counts = read_csv("M:/Erik/Data/Omics/TimeCourse/Raw Galaxy Output/Normalized Counts/DESeq2_Normalized_Counts_on_GT_e13_PARA.csv", col_names = TRUE)
  ## Rename columns
colnames(e13_Normalized_Counts) = c("GeneID", "WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3")
  ## Make sure the Enrichr reference results are correctly ordered when subsetting the "2)" iteration.
e13_GOs = e13_GOs[order(e13_GOs$GO_Biological_Process_2018.Adjusted.P.value, decreasing = FALSE), ]
    ## Re-set the index to make sure indeces are correct
row.names(e13_GOs) = NULL

  ## Subset the dataframe (and run the subsequent code all at once) for one of the 3 following parameters:
    ## By the results obtained in the Enrichr search, identifying "neuron differentiation" as a significantly affected biological process:
      ## 1) For all genes in the term
#RNASeqRepResultsNeurDiff = e13_Normalized_Counts[e13_Normalized_Counts$GeneID %in% c(neuron_differentiation),]
      ## 2) For all the genes in the term + Neurog1, Neurog2, Tlx3
RNASeqRepResultsNeurDiff = e13_Normalized_Counts[e13_Normalized_Counts$GeneID %in% c(neuron_differentiation2),]
      ## 3) For only the overlapping (DE) genes
#RNASeqRepResultsNeurDiff = e13_Normalized_Counts[e13_Normalized_Counts$GeneID %in% c(e13_GO_GENES[6,1:e13_GOs$Overlap[6]]),]
      ## 4) For only DE downregulated overlapping genes
#RNASeqRepResultsNeurDiff = e13_Normalized_Counts[e13_Normalized_Counts$GeneID %in% c(DESeq2_e133$GeneID[Neuron_Diff_Overlap[Down_Genes_indeces_in_Neuron_Diff]]),]

      ## 5) For all DEGs in the DESeq2 dataset
#RNASeqRepResultsNeurDiff = e13_Normalized_Counts[e13_Normalized_Counts$GeneID %in% c(DESeq2_e133$GeneID), ]


  ## Create dataframes to calculate Z-scores across replicates of each gene.
    ## These will hold the original data frame values until filled in later commands
  ## "Big" will store the means and sds of each gene
RNASeqRepResultsNeurDiffBIG = RNASeqRepResultsNeurDiff
  ## "3" will store only the means by genes for all replicates
RNASeqRepResultsNeurDiff3 = RNASeqRepResultsNeurDiff
  ## "4" will store only the sds by gene for all replicates
RNASeqRepResultsNeurDiff4 = RNASeqRepResultsNeurDiff
  ## "Z" will store only the Z-scores, which will be used directly to create the heatmaps
RNASeqRepResultsNeurDiffZ = RNASeqRepResultsNeurDiff

  ## Create a new column to fill with the means and standard deviations of each gene's transcript counts per million ("TPM")
RNASeqRepResultsNeurDiffBIG$mean = 0
RNASeqRepResultsNeurDiffBIG$sd = 0

  ## Loop through the data frame and fill the "mean" and "sd" columns with their appropriate values
for (i in 1:nrow(RNASeqRepResultsNeurDiffBIG)){
  RNASeqRepResultsNeurDiffBIG$mean[i] = (sum(RNASeqRepResultsNeurDiff[i,c(2:7)])/ncol(RNASeqRepResultsNeurDiff[,c(2:7)]))
}
  for (j in 1:nrow(RNASeqRepResultsNeurDiffBIG)){
    RNASeqRepResultsNeurDiffBIG$sd[j] = sd(RNASeqRepResultsNeurDiff[j,c(2:7)])
}
  ## Create a dataframe, storing the gene-specific mean of normalized TPM in all columns/replicates for Z-score calculating
RNASeqRepResultsNeurDiff3[,c(2:7)] = RNASeqRepResultsNeurDiffBIG$mean
  ## Create a dataframe, storing the gene-specific normalized TPM standard deviations in all columns/replicates for Z-score calculating
RNASeqRepResultsNeurDiff4[,c(2:7)] = RNASeqRepResultsNeurDiffBIG$sd
  ## Create the Z-score dataframe
RNASeqRepResultsNeurDiffZ[,c(2:7)] = (RNASeqRepResultsNeurDiff[,c(2:7)] - RNASeqRepResultsNeurDiff3[,c(2:7)])/RNASeqRepResultsNeurDiff4[,c(2:7)]
  ## Remove the genes that were detected to have 0 TPMs across all samples
    ## Create a function that evaluates a vector/dataframe's (x's) numerical values. Returns equal length vector with T/F bools.
      ## Subset the dataframe that filters those genes.
row_has_na = apply(RNASeqRepResultsNeurDiffZ, 1, function(x){any(is.na(x))})
RNASeqRepResultsNeurDiffZ = RNASeqRepResultsNeurDiffZ[!row_has_na,]

  ## Carefully subset the list of interest; this changes depending on what is desired to be shown
#RNASeqRepResultsNeurDiffZA = RNASeqRepResultsNeurDiffZ %>%
#  filter(GeneID %in% DESeq2_e133$GeneID[c(1:1635,2087)])

RNASeqRepResultsNeurDiffZA = RNASeqRepResultsNeurDiffZ
                                    
  ## Create a list of gene names **ordered by Euclidean distance** by which to re-order the dataframe/heatmap
Euclid_dist_order = hclust(dist(RNASeqRepResultsNeurDiffZA[,c(2:7)], method = "euclidean"))$order
  ## The names (not the numbers)
Euclid_dist_ord_Genes = c(RNASeqRepResultsNeurDiffZA$GeneID[Euclid_dist_order])

  ## Transform again to order by clusters (Euclidean distances)
#RNASeqRepResultsNeurDiffZA = RNASeqRepResultsNeurDiffZA %>%
 # mutate(GeneID =  factor(GeneID, levels = Euclid_dist_ord_Genes)) %>%
#  arrange(GeneID)

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap; give appropriate column names to this new dataframe
RNASeqRepResultsNeurDiff2 = melt(RNASeqRepResultsNeurDiffZ, id = "GeneID")
colnames(RNASeqRepResultsNeurDiff2) = c("GeneID", "Genotype", "TPM Z-score")

  ## Create and store the heatmap core as a variable, with appropriate scale gradient limits (look up the max/min Z-score!)
  ## Use geom_tile to map the transcript counts by Z-score
  ## Fill the scale gradient, with red = downregulation, green = upregulation; Title the legend
  ## Add labels

  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the title and adjust the genes/genotype text sizes and locations
e13_HEATER = ggplot(RNASeqRepResultsNeurDiff2, aes(`Genotype`, `GeneID`)) + 
  geom_tile(aes(fill = `TPM Z-score`), color = "black") + 
  scale_fill_gradient2(low = "navy", high = "gold3", name = "TPM\nZ-score", limits = c(-1.26,1.76)) + 
  labs(x = "Genotype", y = "Downregulated Genes within Process", title = "e13 Tmem184b GT Trx Profile of\n'neuron differentiation' (GO: 0030182)") + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(plot.title = element_text(hjust = 0.5), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), axis.text.y = element_text(vjust = 0.5, size = 8), legend.title = element_text(hjust = 0))

  ## Plot the heatmap
e13_HEATER



  ## Find where the "Neuron Differentiation" DEGs are in the clustered matrix for subsetting the heatmap
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[1]]) # 28
    ## Check 28 in the the re-ordered list containing gene names; should be "Nrbp2"
#Euclid_dist_ord_Genes[28]
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[2]]) # 98
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[3]]) # 61
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[4]]) # 71
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[5]]) # 72
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[6]]) # 76
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[7]]) # 26
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[8]]) # 30
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[9]]) # 32
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[10]]) # 79
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[11]]) # 24
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[12]]) # 80
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[13]]) # 101
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[14]]) # 74
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[15]]) # 23
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[16]]) # 33
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[17]]) # 78
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[18]]) # 34
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[19]]) # 29
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[20]]) # 31
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[21]]) # 58
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[22]]) # 99
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[23]]) # 20
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[24]]) # 52
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[25]]) # 36
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[26]]) # 97
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[27]]) # 73
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[28]]) # 77
#which(Euclid_dist_ord_Genes == DESeq2_e133$GeneID[Neuron_Diff_Overlap[29]]) # 62

  ## Concatenate the list of indeces that refer to genes in the overlap list (Neuron_Diff_Overlap) within the re-ordered-by-clustered-Euclidean-distance list of genes. 
    ## This should read out the list of genes that are DE within the Neuron Differentiation GO Bio Process + Tlx3, Neurog1, Neurog2
Euclid_dist_ord_Genes[c(28,98,61,71,72,76,26,30,32,79,24,80,101,74,23,33,78,34,29,31,58,99,20,52,36,97,73,77,62)]
  ## Concatenate the list of indeces that refer to the indeces of the Zscore dataframe (RNASeqRepResultsNeurDiffZA) where those overlap genes are.
    ## This should read out the indeces within the Zscore dataframe of each gene in the overlap list.
Euclid_dist_order[c(28,98,61,71,72,76,26,30,32,79,24,80,101,74,23,33,78,34,29,31,58,99,20,52,36,97,73,77,62)]
      ## Check by subsetting the Zscore dataframe
RNASeqRepResultsNeurDiffZA[Euclid_dist_order[c(28,98,61,71,72,76,26,30,32,79,24,80,101,74,23,33,78,34,29,31,58,99,20,52,36,97,73,77,62)],]

  ## Create the pheatmap for the DEGs involved in "neuron differentiation"
pheatmap(mat = RNASeqRepResultsNeurDiffZA[Euclid_dist_order[c(28,98,61,71,72,76,26,30,32,79,24,80,101,74,23,33,78,34,29,31,58,99,20,52,36,97,73,77,62)] ,2:7], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, cutree_rows = 2, labels_row = c(Euclid_dist_ord_Genes[c(28,98,61,71,72,76,26,30,32,79,24,80,101,74,23,33,78,34,29,31,58,99,20,52,36,97,73,77,62)]))



##### Adult Data prep #####
  ## Import the dataset you want to analyze
DESeq2_Adults = read.csv("M:/Erik/Data/Omics/RNAseq/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
DESeq2_Adults3 = subset(DESeq2_Adults, (!is.na(DESeq2_Adults[,"AdjP"])))
  ## Include Nppb
DESeq2_Adults9 = DESeq2_Adults[c(1:12955,25308),]
  ## Put a placeholder Q value in for now
DESeq2_Adults9[12956,7] = 0.9999999
  ## Re-order the index
rownames(DESeq2_Adults9) = NULL

  ## Perform relevant analysis through enrichr; in this case, FDR <= 0.01
    ## first: GO Bio Processes; this analyzes genes from your dataset in enrichr, which draws terms from geneontology.org
Adult_GOs = as.data.frame(
  enrichr(
    c(DESeq2_Adults3$GeneID[which(DESeq2_Adults3$AdjP <= 0.05)]), DBs[1])
)
    ## second: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
Adult_PATHWAYS = as.data.frame(
  enrichr(
    c(DESeq2_Adults3$GeneID[which(DESeq2_Adults3$AdjP <= 0.05)]), DBs[2])
)

  ## Remove unnecessary columns
Adult_GOs = Adult_GOs[,-c(5,6,7)]
Adult_PATHWAYS = Adult_PATHWAYS[,-c(5,6,7)]

  ## Clean up both dataframes and export both dataframes.
    ## The pathway analysis dataframe needs additional cleaning in Python;
    ## Following that manipulation, data will be re-imported in other scripts for visualization using ggplot
Adult_GOs = Adult_GOs %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Adult_GOs$Overlap = as.numeric(Adult_GOs$Overlap)
Adult_GOs$`Process Size` = as.numeric(Adult_GOs$`Process Size`)

    ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Adult_GOs = add_column(Adult_GOs, EnrichrZscore = Adult_GOs$GO_Biological_Process_2018.Combined.Score/log(Adult_GOs$GO_Biological_Process_2018.P.value), .before = 6)
    ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
      ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Adult_GOs = add_column(Adult_GOs, Weighted_Overlap_Ratio = Adult_GOs$Overlap*(Adult_GOs$Overlap/Adult_GOs$`Process Size`), .before = 4)
    ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
Adult_GOs = add_column(Adult_GOs, Modified_Combined_Score = (abs(Adult_GOs$EnrichrZscore))*(-log(Adult_GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)


## Repeat the above code for Panther Pathways
Adult_PATHWAYS = Adult_PATHWAYS %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Adult_PATHWAYS$Overlap = as.numeric(Adult_PATHWAYS$Overlap)
Adult_PATHWAYS$`Process Size` = as.numeric(Adult_PATHWAYS$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Adult_PATHWAYS = add_column(Adult_PATHWAYS, EnrichrZscore = Adult_PATHWAYS$Panther_2016.Combined.Score/log(Adult_PATHWAYS$Panther_2016.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Adult_PATHWAYS = add_column(Adult_PATHWAYS, Weighted_Overlap_Ratio = Adult_PATHWAYS$Overlap*(Adult_PATHWAYS$Overlap/Adult_PATHWAYS$`Process Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
Adult_PATHWAYS = add_column(Adult_PATHWAYS, Modified_Combined_Score = (abs(Adult_PATHWAYS$EnrichrZscore))*(-log(Adult_PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)


Adult_GO_GENES = str_split(Adult_GOs$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Adult_GO_GENES)){
  Adult_GO_GENES[i,] = paste(substr(Adult_GO_GENES[i,], 1, 1),
                    tolower(substr(Adult_GO_GENES[i,], 2, 7)), sep = "") 
}
Adult_PATHWAY_GENES = str_split(Adult_PATHWAYS$Panther_2016.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Adult_PATHWAY_GENES)){
  Adult_PATHWAY_GENES[i,] = paste(substr(Adult_PATHWAY_GENES[i,], 1, 1),
                             tolower(substr(Adult_PATHWAY_GENES[i,], 2, 7)), sep = "") 
}

## ...In progress... Gene - Term associations..
#Genes_and_Terms = data.frame(1:nrow(DESeq2_Expression_Results3), 2)
#colnames(Genes_and_Terms) = c("GeneID", "Terms")
#Genes_and_Terms$GeneID = ""
#Genes_and_Terms$GeneID = as.character(DESeq2_Expression_Results3$GeneID)
#Genes_and_Terms$Terms = ""

  ## Export data as .CSVs to the server
#write.table(Adult_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\Adult GT DRG GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(Adult_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\Adult GT DRG Panther Pathways.csv", sep = ",", col.names = T)

##### Volcano Plot of Adult GT DRG #####

  ## Create a column in the DESeq2 dataframe that scales the Adjusted P-value by log-base 10
DESeq2_Adults9$log10Pval = -log10(DESeq2_Adults9$AdjP)

  ## Create the data points of interest for functional display.. ggplot struggles to plot subsets of dataframes, so creating a variable that directly labels all the points accordingly helps.
  ## Create a new column (genes of interest); Using "GeneID" is irrelevant; the values will be replaced in the next step
DESeq2_Adults9$g.o.i. = DESeq2_Adults9$GeneID

  ## Fill the points with appropriately indexed data;
  ## "Differentially Expressed Genes" are defined as being "genes of interest" with Adjusted P-Values of < 0.05
DESeq2_Adults9$g.o.i.[which(DESeq2_Adults9$AdjP <= 0.05)]= "DEGs"

  ## Optional third color labeling
#DESeq2_Adults3$g.o.i.[c(1:405)][-c()] = ""

  ## "Non-differentially Expressed Genes" are defined as being "genes of interest" with Adjusted P-values > 0.05
DESeq2_Adults9$g.o.i.[which(DESeq2_Adults9$AdjP >= 0.05)]= "Non-DEGs"

  ## Subset of genes with well-known functions and effects on DRG pain/itch behavior, according to literature
#DESeq2_Adults9$g.o.i.[c(2,5,8,9,13,14,21,31,34,46,54,74,80,101,118,128,145,208,250,268,352,353,354,356)] = "Itch-related DEGs"

which(DESeq2_Adults9$GeneID == "Trpm6")

differentiation = c(Find_Genes_With_X_GO_Term(GO_Term = "differentiation"))

DESeq2_Adults9$g.o.i.[diff_indx] = "Differentiation"

  ## Create labels on points we're interested in labeling in volcano plots
  ## Create the column, "labs"; hide all of the text labels with: ""
DESeq2_Adults9$labs = ""


  ## Label only these selected items; adjust for each experiment; for the current example and a more global volcano plot, highlight the most extreme data points on the plot
#ix_label1 = c(DESeq2_Adults3$GeneID[DESeq2_Adults3$log10Pval > 28])

  ## Fill the appropriate rows of the "labs" column with the corresponding gene names for labeling those points
DESeq2_Adults9$labs[c(c(1,2,5,8,9,13,21,31,54,118,352,356))] = DESeq2_Adults9$GeneID[c(c(1,2,5,8,9,13,21,31,54,118,352,356))]


  ## Notes on ggplot commands:
  ## geom-point = alpha is a ggplot add-on function that controls point color transparency
  ## coord_cartesian is a ggplot add-on function that determines axes sizes
  ## labs is a ggplot add-on function that determines the plot labels; x-axis, y-axis, title, legend colors
  ## scale_color_manual is a ggplot add-on function that sets the colors of the plot, corresponding to the groups determined in the the base function ("col" in aes)
  ## Theme is used to manipulate the physical location of the plot title


  ## Pay attention to the exact strings of column names.
voladult = ggplot(DESeq2_Adults9) +
  geom_point(data = subset(DESeq2_Adults9, `g.o.i.` == "Non-DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.2) +
  #geom_point(data = subset(DESeq2_e133, `g.o.i.` == "Non-DEGs in Neuron Differentiation"),
  #           aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.1) +
  geom_point(data = subset(DESeq2_Adults9, `g.o.i.` == "DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.3) +
  geom_point(data = subset(DESeq2_Adults9, `g.o.i.` == "Stimulus Response\nDEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.9) +
  coord_cartesian(xlim = c(-3.09,1.6), ylim = c(0,52)) +
  labs(title = expression(paste(italic("Tmem184b")^"GT/GT", " DRG mRNA Exp. Changes")), x = "log2 (FC) rel. to WT", y = "-log10 (Adj.P-val)", col = "Gene Data Type") +
  #labs(x = "log2 (FC) rel. to WT", y = "-log10(Adj.P-val)", title = "Adult Tmem184b GT DRG Gene Exp. Changes", col = "Gene Data Type") +
  geom_hline(yintercept = 1.3, type = "solid", color = "black") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.80,0.87)) + 
  scale_color_manual(values = c("darkgoldenrod4", "gray48", "navy"))

  ## geom_text_repel labels the selected points in the color commanded
  ## Re-save the variable with the added functionality; for unlabeled plot
voladult = voladult  + geom_text_repel(data = DESeq2_Adults9, x = DESeq2_Adults9$log2.FC., y = DESeq2_Adults9$log10Pval, color = "black", aes(label = DESeq2_Adults9$labs)) 

  ## Plot it
voladult


##### Adult MA plot #####
  ## Labels for MA
DESeq2_Adults3$labs[c(c(1,2,13,54,118,352,356,965,3166))] = DESeq2_Adults3$GeneID[c(c(1,2,13,54,118,352,356,965,3166))]
  ## Transform the mean of normalized counts
DESeq2_Adults3$Base.Mean = log(DESeq2_Adults3$Base.Mean)

Adult_MA = ggplot(DESeq2_Adults3) +
  geom_point(data = subset(DESeq2_Adults3, `g.o.i.` == "Non-DEGs"),
             aes(x = `Base.Mean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.2) +
  geom_point(data = subset(DESeq2_Adults3, `g.o.i.` == "DEGs"), 
             aes(x = `Base.Mean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.3) +
  geom_point(data = subset(DESeq2_Adults3, `g.o.i.` == "Itch-related DEGs"),
             aes(x = `Base.Mean`, y = `log2.FC.`, color = `g.o.i.`), alpha = 0.8) +
  coord_cartesian(xlim = c(0,13), ylim = c(-4,1.5)) +
  labs(x = "ln (Mean of TPM)", y = "log2 (FC) mRNA Rel. to WT", title = "Adult DRG Tmem184b GT mRNA Exp.", col = "Gene Data Type") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("darkgoldenrod4", "navy", "gray48"))

  ## Add labels
Adult_MA = Adult_MA + geom_text_repel(data = DESeq2_Adults3, x = DESeq2_Adults3$Base.Mean, y = DESeq2_Adults3$log2.FC., color = "black", aes(label = `labs`))

  ## Plot it
Adult_MA


##### Adult Bar Plots of GO Terms Using ggplot2 #####
  ## Sort by -log10(Q). Retain the top 10 Processes
Adult_GOs$neglog10Q = -log10(Adult_GOs$GO_Biological_Process_2018.Adjusted.P.value)

Adult_GOs = Adult_GOs[order(Adult_GOs$neglog10Q, decreasing = TRUE),]

Adult_GO_Bar_ADJP = Adult_GOs[1:10,]
  ## Remove the "GO...." strings in the GO Process column
Adult_GO_Bar_ADJP$GO_Biological_Process_2018.Term = Adult_GO_Bar_ADJP$GO_Biological_Process_2018.Term %>% gsub(x = Adult_GO_Bar_ADJP$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")

  ## Shorten strings of terms
Adult_GO_Bar_ADJP$GO_Biological_Process_2018.Term[1] = "positive regulation of\nneuron projection development"
Adult_GO_Bar_ADJP$GO_Biological_Process_2018.Term[2] = "positive regulation of CD4-positive,\nalpha-beta T cell activation"
Adult_GO_Bar_ADJP$GO_Biological_Process_2018.Term[3] = "G-protein coupled receptor\n signaling pathway, coupled to \ncyclic nucleotide second messenger"
Adult_GO_Bar_ADJP$GO_Biological_Process_2018.Term[7] = "negative regulation of calcium\n ion-dependent exocytosis"


  ## Plot the bar plot
  ## Store as a variable
  ## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long GO terms
Adult_GO_BAR_ADJP = ggplot(Adult_GO_Bar_ADJP, aes(x = reorder(`GO_Biological_Process_2018.Term`, `neglog10Q`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "2018 GO Biological Process", y = "Enrichr -log10(Adj. P-val)", title = "Adult Tmem184b GT/GT DRG GO Analysis") + 
  theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  coord_flip()

## Plot it.
Adult_GO_BAR_ADJP


##### Adult Bar Plots of Panther Pathways Using ggplot2 #####
  ## Convert the Q-values to -log10 scale
Adult_PATHWAYS$neglog10Q = -log10(Adult_PATHWAYS$Panther_2016.Adjusted.P.value)
  ## Remove all but the top 15 pathways
Adult_PATHWAYS = Adult_PATHWAYS[order(Adult_PATHWAYS$neglog10Q, decreasing = TRUE),]

Adult_Pathway_Bar_ADJP = Adult_PATHWAYS[1:10,]
  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
Adult_Pathway_Bar_ADJP$Panther_2016.Term = Adult_Pathway_Bar_ADJP$Panther_2016.Term %>% gsub(x = Adult_Pathway_Bar_ADJP$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")
  
  ## Shorten term strings
Adult_Pathway_Bar_ADJP$Panther_2016.Term[1] = "Heterotrimeric G-protein signaling pathway-\nGq alpha and Go alpha mediated pathway"
Adult_Pathway_Bar_ADJP$Panther_2016.Term[2] = "Histamine H1 receptor\nmediated signaling pathway"
Adult_Pathway_Bar_ADJP$Panther_2016.Term[4] = "Thyrotropin-releasing hormone\nreceptor signaling pathway"
Adult_Pathway_Bar_ADJP$Panther_2016.Term[5] = "Oxytocin receptor\nmediated signaling pathway"
Adult_Pathway_Bar_ADJP$Panther_2016.Term[6] = "Muscarinic acetylcholine receptor\n1 and 3 signaling pathway"
Adult_Pathway_Bar_ADJP$Panther_2016.Term[7] = "5HT2 type receptor\nmediated signaling pathway"
Adult_Pathway_Bar_ADJP$Panther_2016.Term[8] = "Angiotensin II-stimulated signaling\nthrough G proteins and beta-arrestin"
Adult_Pathway_Bar_ADJP$Panther_2016.Term[10] = "Alzheimer disease-\namyloid secretase pathway"


  ## Plot the bar plot
  ## Store as a variable
  ## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long pathway names
Adult_PATHWAY_BAR_ADJP = ggplot(Adult_Pathway_Bar_ADJP, aes(x = reorder(`Panther_2016.Term`, `neglog10Q`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "Panther 2016 Pathways", y = "Enrichr -log10(Adj. P-val)", title = "Adult Tmem184b GT/GT DRG Pathway Analysis") + 
  theme_classic() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")  + 
  coord_flip()

  ## Plot it.
Adult_PATHWAY_BAR_ADJP



##### Adult Panther Pathways Heatmap with labels using ggplot2 #####
  ## Import the relevant data
DEGs_and_Pathways =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr Pathway Clustergram Adult GT.csv")
  ## Export the gene list for photoshopping
#write.csv(DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\aDRG FDR 05 Heatmap DEGs.csv", sep = ",")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
Adult_Binom_Heatmap = melt(DEGs_and_Pathways, id = "DEGs")
colnames(Adult_Binom_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

  ## Remove excessive pathway name strings
Pathway_Names = Adult_Binom_Heatmap$Pathways %>% gsub(x = Adult_Binom_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

  ## Put the new strings back into the dataframe
Adult_Binom_Heatmap[,2] = Pathway_Names

  ## Make the digital values factors so it's compatible to graph
Adult_Binom_Heatmap$`Presence In Pathway` = as.factor(Adult_Binom_Heatmap$`Presence In Pathway`)



  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

HEATlabs = ggplot(Adult_Binom_Heatmap, aes(`Pathways`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Pathway`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.05)", title = "D.E.G.s of Affected Pathways in Tmem184b GT/GT aDRG") + 
  guides(fill = guide_legend(title = "Presence in Pathway")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


## Plot the heatmap
HEATlabs




## Import the relevant data
DEGs_and_GOs =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr GO Clustergram Adult GT.csv")
  ## Export the gene list for photoshopping
#write.csv(DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\aDRG FDR 05 Heatmap DEGs.csv", sep = ",")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
Adult_Binom_Heatmap = melt(DEGs_and_GOs, id = "DEGs")
colnames(Adult_Binom_Heatmap) = c("DEGs", "GOs", "Presence In Process")
  ## Remove the "GO...." strings in the GO Process column
Adult_Binom_Heatmap$GOs = Adult_Binom_Heatmap$GOs %>% gsub(x = Adult_Binom_Heatmap$GOs, pattern = " \\(.*\\)$", replacement = "")
  ## Make the digital values factors so it's compatible to graph
Adult_Binom_Heatmap$`Presence In Process` = as.factor(Adult_Binom_Heatmap$`Presence In Process`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

HEATlabs = ggplot(Adult_Binom_Heatmap, aes(`GOs`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Process`), color = "black") + scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.05)", title = "D.E.G.s of Affected GO Processes in Tmem184b GT/GT aDRG") + 
  guides(fill = guide_legend(title = "Presence in Process")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


## Plot the heatmap
HEATlabs


##### Adult Panther Pathways Heatmap without labels using ggplot2 #####

## Without labels
#HEATnolabs = ggplot(Adult_Binom_Heatmap, aes(`Pathways`, `DEGs`)) + geom_tile(aes(fill = `Presence In Pathway`), color = "black") + scale_fill_manual(values = c("navy", "gold3")) + labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.05)", title = "D.E.G.s of Affected Pathways in Tmem184b GT/GT aDRG") + guides(fill = guide_legend(title = "Presence in Pathway")) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0), axis.text= element_blank()) + coord_flip()

## Plot the heatmap
#HEATnolabs



##### Create a heatmap across replicates #####
  ## Import the adult normalized counts file. These are expression estimates for each gene, for each sample/replicate, where each gene's value is normalized to its sample's effect size
Adult_Normalized_Counts = read_csv("M:/Erik/Data/Omics/RNAseq/Processed Galaxy Output/Counts Files to Upload/RNASeqRepResults.csv", col_names = TRUE)
  ## Rename columns
colnames(Adult_Normalized_Counts) = c("GeneID", "WT1", "WT2", "WT3", "WT4", "Mut1", "Mut2", "Mut3", "Mut4")

RNASeqRepResultsAdultAll = Adult_Normalized_Counts[Adult_Normalized_Counts$GeneID %in% c(DESeq2_Adults9$GeneID),]

  ## Create dataframes to calculate Z-scores across replicates of each gene.
    ## These will hold the original data frame values until filled in later commands
      ## "Big" will store the means and sds of each gene
RNASeqRepResultsAdultAllBIG = RNASeqRepResultsAdultAll
      ## "3" will store only the means by genes for all replicates
RNASeqRepResultsAdultAll3 = RNASeqRepResultsAdultAll
      ## "4" will store only the sds by gene for all replicates
RNASeqRepResultsAdultAll4 = RNASeqRepResultsAdultAll
      ## "Z" will store only the Z-scores, which will be used directly to create the heatmaps
RNASeqRepResultsAdultAllZ = RNASeqRepResultsAdultAll

  ## Create a new column to fill with the means and standard deviations of each gene's transcript counts per million ("TPM")
RNASeqRepResultsAdultAllBIG$mean = 0
RNASeqRepResultsAdultAllBIG$sd = 0

  ## Loop through the dataframe and fill the "mean" and "sd" columns with their appropriate values
for (i in 1:nrow(RNASeqRepResultsAdultAllBIG)){
  RNASeqRepResultsAdultAllBIG$mean[i] = (sum(RNASeqRepResultsAdultAll[i,c(2:9)])/ncol(RNASeqRepResultsAdultAll[,c(2:9)]))
}
  for (j in 1:nrow(RNASeqRepResultsAdultAllBIG)){
    RNASeqRepResultsAdultAllBIG$sd[j] = sd(RNASeqRepResultsAdultAll[j,c(2:9)], na.rm = TRUE)
}
    ## Create a dataframe storing the gene-specific mean normalized TPM in all columns/replicates for Z-score calculating
RNASeqRepResultsAdultAll3[,c(2:9)] = RNASeqRepResultsAdultAllBIG$mean
    ## Create a dataframe storing the gene-specific normalized TPM standard deviationsin all columns/replicates for Z-score calculating
RNASeqRepResultsAdultAll4[,c(2:9)] = RNASeqRepResultsAdultAllBIG$sd
    ## Create the Z-score dataframe
RNASeqRepResultsAdultAllZ[,c(2:9)] = (RNASeqRepResultsAdultAll[,c(2:9)] - RNASeqRepResultsAdultAll3[,c(2:9)])/RNASeqRepResultsAdultAll4[,c(2:9)]
    ## Remove the genes that were detected to have 0 TPMs across all samples
    ## Create a function that evaluates a vector/dataframe's (x's) numerical values. Returns equal length vector with T/F bools.
    ## Subset the dataframe that filters those genes.
row_has_na = apply(RNASeqRepResultsAdultAllZ, 1, function(x){any(is.na(x))})
RNASeqRepResultsAdultAllZ = RNASeqRepResultsAdultAllZ[!row_has_na,]

  ## Transform the profile to include only DEGs
RNASeqRepResultsAdultAllZA = RNASeqRepResultsAdultAllZ %>%
 filter(GeneID %in% DESeq2_Adults9$GeneID[c(1:405,12956)])


  ## Create a list of gene names ordered by euclidean distance by which to re-order the dataframe
Euclid_dist_order = hclust(dist(RNASeqRepResultsAdultAllZA[,c(2:9)], method = "euclidean"))$order
    ## The names (not the numbers)
Euclid_dist_ord_Genes = c(RNASeqRepResultsAdultAllZA$GeneID[Euclid_dist_order])

  ## Transform again to order by clusters (Euclidean distances)
RNASeqRepResultsAdultAllZA = RNASeqRepResultsAdultAllZA %>%
  mutate(GeneID =  factor(GeneID, levels = Euclid_dist_ord_Genes)) %>%
  arrange(GeneID)

    ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap; give appropriate column names to this new dataframe
RNASeqRepResultsAdultAll2 = melt(RNASeqRepResultsAdultAllZA, id = "GeneID")
colnames(RNASeqRepResultsAdultAll2) = c("GeneID", "Genotype", "TPM Z-score")

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the transcript counts by Z-score
  ## Fill the scale gradient, with red = downregulation, green = upregulation; Title the legend
  ## Add labels
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the title and adjust the genes/genotype text sizes and locations
AdultHEATER = ggplot(RNASeqRepResultsAdultAll2, aes(`Genotype`, `GeneID`)) + 
  geom_tile(aes(fill = `TPM Z-score`)) + 
  scale_fill_gradient2(low = "navy", mid = "white", high = "darkgoldenrod4", name = "TPM\nZ-score") + # limits = c(-1.850,2.251)) + 
  labs(x = "Genotype", y = "DEGs", title = "Adult Tmem184b GT DRG Trx Profile of DEGs") + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(plot.title = element_text(hjust = 0.5), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(vjust = 0.5, size = 3), legend.title = element_text(hjust = 0))

#, panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1)

## Plot the heatmap
AdultHEATER


  ## Find where the "Itch-related DEGs" are in the clustered matrix subsetted to DEGs for the heatmap

which(Euclid_dist_ord_Genes == "Il31ra") # 1534
which(Euclid_dist_ord_Genes == "Cysltr2") # 1533
which(Euclid_dist_ord_Genes == "Npy2r") # 1428
which(Euclid_dist_ord_Genes == "Sst") # 1421
which(Euclid_dist_ord_Genes == "Htr1a") # 1536
which(Euclid_dist_ord_Genes == "P2rx3") # 1565
which(Euclid_dist_ord_Genes == "Lpar3") # 1911
which(Euclid_dist_ord_Genes == "Lpar5") # 1341
which(Euclid_dist_ord_Genes == "Scn11a") # 1923
which(Euclid_dist_ord_Genes == "Scn10a") # 1795
which(Euclid_dist_ord_Genes == "Mrgprd") # 1922
which(Euclid_dist_ord_Genes == "Trpc6") # 905
which(Euclid_dist_ord_Genes == "Trpc3") # 1752
which(Euclid_dist_ord_Genes == "F2rl2") # 1939
which(Euclid_dist_ord_Genes == "Htr1f") # 1530
which(Euclid_dist_ord_Genes == "Osmr") # 3193
which(Euclid_dist_ord_Genes == "Fxyd2") # 9568
which(Euclid_dist_ord_Genes == "Htr4") # 9577
which(Euclid_dist_ord_Genes == "Mrgprx1") # 9571
which(Euclid_dist_ord_Genes == "Ptgdr") # 909
which(Euclid_dist_ord_Genes == "Trpa1") # 3993
which(Euclid_dist_ord_Genes == "Trpm6") # 1705
which(Euclid_dist_ord_Genes == "Hrh1") # 1543
which(Euclid_dist_ord_Genes == "Mrgpra3") # 1413
which(Euclid_dist_ord_Genes == "Nppb") # 1422
which(Euclid_dist_ord_Genes == "Tmem184b") # 1537

  ## Confirm by indexing; copy and paste to customize the clustered heatmap
Euclid_dist_ord_Genes[c(1534,1533,1428,1421,1536,1565,1911,1341,1923,1795,1922,905,1752,1939,1530,3193,9568,9577,9571,909,3993,1705,1543,1413,1422)]
  ## Confirm the location/order of the transformed matrix/dataframe
    ## Should be "Il31ra"
RNASeqRepResultsAdultAllZA[1534,]


  ## Visualize the Z-scored, Euclidean-clustered and ordered gene TPMs across replicates with "pheatmap"
pheatmap(mat = RNASeqRepResultsAdultAllZA[,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, show_rownames = F)
  ## Visualize around the Tmem184b cluster
pheatmap(mat = RNASeqRepResultsAdultAllZA[1530:1544,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, labels_row = Euclid_dist_ord_Genes[1530:1544])
  ## Visualize the "Itch-related DEGs"
pheatmap(mat = RNASeqRepResultsAdultAllZA[c(1534,1533,1428,1421,1536,1565,1911,1341,1923,1795,1922,905,1752,1939,1530,3193,9568,9577,9571,909,3993,1705,1543,1413,1422),2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, cutree_rows = 5, labels_row = Euclid_dist_ord_Genes[c(1534,1533,1428,1421,1536,1565,1911,1341,1923,1795,1922,905,1752,1939,1530,3193,9568,9577,9571,909,3993,1705,1543,1413,1422)])
  ## Visualize the "Itch-related DEGs"
pheatmap(mat = RNASeqRepResultsAdultAllZA[c(1420:1560),2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, cutree_rows = 5, labels_row = Euclid_dist_ord_Genes[c(1400:1600)])

## Visualize the "Itch-related DEGs"
pheatmap(mat = RNASeqRepResultsAdultAllZA[c(1534,1533,1428,1421,1536,1565,1911,1341,1923,1795,1922,905,1752,1939,1530,3193,9568,9577,9571,909,3993,1705,1543,1413,1422,1537),2:9], clustering_distance_rows = "euclidean", color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), angle_col = 0, treeheight_row = 35, treeheight_col = 9, cutree_rows = 4, labels_row = Euclid_dist_ord_Genes[c(1534,1533,1428,1421,1536,1565,1911,1341,1923,1795,1922,905,1752,1939,1530,3193,9568,9577,9571,909,3993,1705,1543,1413,1422,1537)])




  ## Find where the "Itch-related DEGs" are in the clustered matrix subsetted to DEGs for the heatmap
which(Euclid_dist_ord_Genes == "Il31ra") # 132
which(Euclid_dist_ord_Genes == "Cysltr2") # 131
which(Euclid_dist_ord_Genes == "Npy2r") # 168
which(Euclid_dist_ord_Genes == "Sst") # 120
which(Euclid_dist_ord_Genes == "Htr1a") # 357
which(Euclid_dist_ord_Genes == "P2rx3") # 157
which(Euclid_dist_ord_Genes == "Lpar3") # 167
which(Euclid_dist_ord_Genes == "Lpar5") # 361
which(Euclid_dist_ord_Genes == "Scn11a") # 387
which(Euclid_dist_ord_Genes == "Scn10a") # 395
which(Euclid_dist_ord_Genes == "Mrgprd") # 386
which(Euclid_dist_ord_Genes == "Trpc6") # 403
which(Euclid_dist_ord_Genes == "Trpc3") # 262
which(Euclid_dist_ord_Genes == "F2rl2") # 244
which(Euclid_dist_ord_Genes == "Htr1f") # 187
which(Euclid_dist_ord_Genes == "Osmr") # 123
which(Euclid_dist_ord_Genes == "Fxyd2") # 336
which(Euclid_dist_ord_Genes == "Htr4") # 404
which(Euclid_dist_ord_Genes == "Mrgprx1") # 196
which(Euclid_dist_ord_Genes == "Ptgdr") # 406
which(Euclid_dist_ord_Genes == "Trpa1") # 299
which(Euclid_dist_ord_Genes == "Trpm6") # 347
which(Euclid_dist_ord_Genes == "Hrh1") # 259
which(Euclid_dist_ord_Genes == "Mrgpra3") # 180
which(Euclid_dist_ord_Genes == "Nppb") # 339
which(Euclid_dist_ord_Genes == "Tmem184b") # 358

  ## Confirm by indexing; copy and paste to customize the clustered heatmap
Euclid_dist_ord_Genes[c(132,131,168,120,357,157,167,361,387,395,386,403,262,244,187,123,336,404,196,406,299,347,259,180,339)]
  ## Confirm the location/order of the transformed matrix/dataframe
    ## Should be "Il31ra"
RNASeqRepResultsAdultAllZA[132,]

  ## Visualize the Z-scored, Euclidean-clustered and ordered gene TPMs across replicates with "pheatmap"
pheatmap(mat = RNASeqRepResultsAdultAllZA[,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, cutree_rows = 2, treeheight_row = 35, treeheight_col = 7, show_rownames = F)
  ## Visualize around the Tmem184b cluster
pheatmap(mat = RNASeqRepResultsAdultAllZA[353:363,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, labels_row = Euclid_dist_ord_Genes[353:363])
  ## Visualize the "Itch-related DEGs"
pheatmap(mat = RNASeqRepResultsAdultAllZA[c(132,131,168,120,357,157,167,361,387,395,386,403,262,244,187,123,336,404,196,406,299,347,259,180,339,358),2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustsering_distance_rows = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 9, cutree_rows = 4, labels_row = Euclid_dist_ord_Genes[c(132,131,168,120,357,157,167,361,387,395,386,403,262,244,187,123,336,404,196,406,299,347,259,180,339,358)])




max(RNASeqRepResultsAdultAllZA[,2:9])
min(RNASeqRepResultsAdultAllZA[,2:9])
max(RNASeqRepResultsAdultAll2[,3])
min(RNASeqRepResultsAdultAll2[,3])

##### In vitro Data prep #####
  ## Import the dataset you want to analyze
DESeq2_In_vitro_DIV14 = read.csv("M:/Erik/Data/Omics/TimeCourse/Processed Galaxy Output/Test Results to Upload/Cultured Embryos for Genotype at DIV14.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
DESeq2_In_vitro_DIV143 = subset(DESeq2_In_vitro_DIV14, (!is.na(DESeq2_In_vitro_DIV14[,"AdjP"])))

  ## Create a column in the DESeq2 dataframe that scales the Adjusted P-value by log-base 10
DESeq2_In_vitro_DIV143$log10Pval = -log10(DESeq2_In_vitro_DIV143$AdjP)

  ## Find the mean 
#av = subset(DESeq2_In_vitro_DIV143, DESeq2_In_vitro_DIV143$AdjP < 0.01 & DESeq2_In_vitro_DIV143$log2.FC. < 0)
#mean(av$log2.FC.)

  ## Create the data points of interest for functional display.. ggplot struggles to plot subsets of dataframes, so creating a variable that directly labels all the points accordingly helps.
    ## Create a new column (genes of interest); Using "GeneID" is irrelevant; the values will be replaced in the next step
DESeq2_In_vitro_DIV143$g.o.i. = DESeq2_In_vitro_DIV143$GeneID

  ## Fill the points with appropriately indexed data
DESeq2_In_vitro_DIV143$g.o.i.[which(c(DESeq2_In_vitro_DIV143$AdjP < 0.01 & DESeq2_In_vitro_DIV143$log2.FC. > 0))]= "Up DEGs"
DESeq2_In_vitro_DIV143$g.o.i.[which(c(DESeq2_In_vitro_DIV143$AdjP < 0.01 & DESeq2_In_vitro_DIV143$log2.FC. < 0))]= "Down DEGs"

#DESeq2_In_vitro_DIV143$g.o.i.[which(c(DESeq2_In_vitro_DIV143$AdjP < 0.01 & (DESeq2_In_vitro_DIV143$log2.FC. > -0.7325703 & DESeq2_In_vitro_DIV143$log2.FC. < 0.7325703)))] = "Stable DEGs"

DESeq2_In_vitro_DIV143$g.o.i.[which(c(DESeq2_In_vitro_DIV143$AdjP > 0.01))] = "Non-DEGs"

##### In vitro volcano plot of GT DRG #####
  
  ## Notes on ggplot commands:
    ## geom-point = alpha is a ggplot add-on function that controls point color transparency
    ## coord_cartesian is a ggplot add-on function that determines axes sizes
    ## labs is a ggplot add-on function that determines the plot labels; x-axis, y-axis, title, legend colors
    ## scale_color_manual is a ggplot add-on function that sets the colors of the plot, corresponding to the groups determined in the the base function ("col" in aes)
    ## Theme is used to manipulate the physical location of the plot title

  ## Pay attention to the exact strings of column names.
volinvitro = ggplot(DESeq2_In_vitro_DIV143) +
  geom_point(data = subset(DESeq2_In_vitro_DIV143, `g.o.i.` == "Up DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.3) +
  geom_point(data = subset(DESeq2_In_vitro_DIV143, `g.o.i.` == "Down DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.3) +
  #geom_point(data = subset(DESeq2_In_vitro_DIV143, `g.o.i.` == "Stable DEGs"),
   #          aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.8) +
  geom_point(data = subset(DESeq2_In_vitro_DIV143, `g.o.i.` == "Non-DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.2) +
  coord_cartesian(xlim = c(-5,5), ylim = c(0,40)) +
  labs(x = "log2(FC)", y = "-log10(Adj.P-val)", title = "In vitro e13 Tmem184b GT DRG Gene Exp. Changes", col = "Gene Data Type") +
  scale_color_manual(values = c("navy", "gray48", "darkgoldenrod4")) +
  #geom_hline(yintercept = 2.0, type = "solid", color = "black") +
  #geom_vline(xintercept = -0.7325703, linetype = "dashed", color = "black") +
  #geom_vline(xintercept = 0.7325703, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none")

  ## Plot it
  volinvitro
  
  
  
##### In vitro MA plot of  GT DRG #####
  ## Transform the mean normalized counts
DESeq2_In_vitro_DIV143$BaseMean = log(DESeq2_In_vitro_DIV143$BaseMean)

  ## Desgin the MA plot
MA = ggplot(data = DESeq2_In_vitro_DIV143, aes(x = `BaseMean`, y = `log2.FC.`, col = `g.o.i.`)) + 
  geom_point(alpha = 0.4) + coord_cartesian(xlim = c(0,13), ylim = c(-5.0,5.0)) + 
  labs(x = "ln(Mean of Normalized Counts", y = "log2(FC)", title = "In vitro e13 Tmem184b GT DRG Gene Expression Changes", col = "Gene Data Type") + 
  scale_color_manual(values = c("navy", "gray87", "gray48", "darkgoldenrod4")) + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

  ## Plot it
MA


##### In vitro data prep for bar plots and heatmaps #####
  ## Perform relevant analysis through enrichr; in this case, FDR <= 0.01
  ## first: GO Bio Processes; this analyzes genes from your dataset in enrichr, which draws terms from geneontology.org
Up_DIV14_GOs = as.data.frame(
  enrichr(
    c(DESeq2_In_vitro_DIV143$GeneID[which(DESeq2_In_vitro_DIV143$g.o.i. == "Up DEGs")]), DBs[1])
)
Down_DIV14_GOs = as.data.frame(
  enrichr(
    c(DESeq2_In_vitro_DIV143$GeneID[which(DESeq2_In_vitro_DIV143$g.o.i. == "Down DEGs")]), DBs[1])
)

  ## second: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
Up_DIV14_PATHWAYS = as.data.frame(
  enrichr(
    c(DESeq2_In_vitro_DIV143$GeneID[which(DESeq2_In_vitro_DIV143$g.o.i == "Up DEGs")]), DBs[2])
)
Down_DIV14_PATHWAYS = as.data.frame(
  enrichr(
    c(DESeq2_In_vitro_DIV143$GeneID[which(DESeq2_In_vitro_DIV143$g.o.i == "Down DEGs")]), DBs[2])
)


  ## Remove unnecessary columns
Up_DIV14_GOs = Up_DIV14_GOs[,-c(5,6,7)]
Down_DIV14_GOs = Down_DIV14_GOs[,-c(5,6,7)]
Up_DIV14_PATHWAYS = Up_DIV14_PATHWAYS[,-c(5,6,7)]
Down_DIV14_PATHWAYS = Down_DIV14_PATHWAYS[,-c(5,6,7)]


## Clean up both dataframes and export both dataframes.
  ## Dataframes need additional cleaning in Python;
    ## Following that manipulation, data will be re-imported in other scripts for visualization using ggplot
Up_DIV14_GOs = Up_DIV14_GOs %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Up_DIV14_GOs$Overlap = as.numeric(Up_DIV14_GOs$Overlap)
Up_DIV14_GOs$`Process Size` = as.numeric(Up_DIV14_GOs$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Up_DIV14_GOs = add_column(Up_DIV14_GOs, EnrichrZscore = Up_DIV14_GOs$GO_Biological_Process_2018.Combined.Score/log(Up_DIV14_GOs$GO_Biological_Process_2018.P.value), .before = 6)

  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Up_DIV14_GOs = add_column(Up_DIV14_GOs, Weighted_Overlap_Ratio = Up_DIV14_GOs$Overlap*(Up_DIV14_GOs$Overlap/Up_DIV14_GOs$`Process Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
Up_DIV14_GOs = add_column(Up_DIV14_GOs, Modified_Combined_Score = (abs(Up_DIV14_GOs$EnrichrZscore))*(-log(Up_DIV14_GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)


Down_DIV14_GOs = Down_DIV14_GOs %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Down_DIV14_GOs$Overlap = as.numeric(Down_DIV14_GOs$Overlap)
Down_DIV14_GOs$`Process Size` = as.numeric(Down_DIV14_GOs$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Down_DIV14_GOs = add_column(Down_DIV14_GOs, EnrichrZscore = Down_DIV14_GOs$GO_Biological_Process_2018.Combined.Score/log(Down_DIV14_GOs$GO_Biological_Process_2018.P.value), .before = 6)

  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Down_DIV14_GOs = add_column(Down_DIV14_GOs, Weighted_Overlap_Ratio = Down_DIV14_GOs$Overlap*(Down_DIV14_GOs$Overlap/Down_DIV14_GOs$`Process Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
Down_DIV14_GOs = add_column(Down_DIV14_GOs, Modified_Combined_Score = (abs(Down_DIV14_GOs$EnrichrZscore))*(-log(Down_DIV14_GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)





  ## Repeat the above code for Panther Pathways analyzed from up DEGs
Up_DIV14_PATHWAYS = Up_DIV14_PATHWAYS %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Up_DIV14_PATHWAYS$Overlap = as.numeric(Up_DIV14_PATHWAYS$Overlap)
Up_DIV14_PATHWAYS$`Process Size` = as.numeric(Up_DIV14_PATHWAYS$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Up_DIV14_PATHWAYS = add_column(Up_DIV14_PATHWAYS, EnrichrZscore = Up_DIV14_PATHWAYS$Panther_2016.Combined.Score/log(Up_DIV14_PATHWAYS$Panther_2016.P.value), .before = 6)

  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Up_DIV14_PATHWAYS = add_column(Up_DIV14_PATHWAYS, Weighted_Overlap_Ratio = Up_DIV14_PATHWAYS$Overlap*(Up_DIV14_PATHWAYS$Overlap/Up_DIV14_PATHWAYS$`Process Size`), .before = 4)

  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
Up_DIV14_PATHWAYS = add_column(Up_DIV14_PATHWAYS, Modified_Combined_Score = (abs(Up_DIV14_PATHWAYS$EnrichrZscore))*(-log(Up_DIV14_PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)


  ## Repeat the above code for Panther Pathways analyzed from down DEGs
Down_DIV14_PATHWAYS = Down_DIV14_PATHWAYS %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Down_DIV14_PATHWAYS$Overlap = as.numeric(Down_DIV14_PATHWAYS$Overlap)
Down_DIV14_PATHWAYS$`Process Size` = as.numeric(Down_DIV14_PATHWAYS$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Down_DIV14_PATHWAYS = add_column(Down_DIV14_PATHWAYS, EnrichrZscore = Down_DIV14_PATHWAYS$Panther_2016.Combined.Score/log(Down_DIV14_PATHWAYS$Panther_2016.P.value), .before = 6)

  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Down_DIV14_PATHWAYS = add_column(Down_DIV14_PATHWAYS, Weighted_Overlap_Ratio = Down_DIV14_PATHWAYS$Overlap*(Down_DIV14_PATHWAYS$Overlap/Down_DIV14_PATHWAYS$`Process Size`), .before = 4)

  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
Down_DIV14_PATHWAYS = add_column(Down_DIV14_PATHWAYS, Modified_Combined_Score = (abs(Down_DIV14_PATHWAYS$EnrichrZscore))*(-log(Down_DIV14_PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)





Up_DIV14_GO_GENES = str_split(Up_DIV14_GOs$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Up_DIV14_GO_GENES)){
  Up_DIV14_GO_GENES[i,] = paste(substr(Up_DIV14_GO_GENES[i,], 1, 1),
                             tolower(substr(Up_DIV14_GO_GENES[i,], 2, 7)), sep = "") 
}
Down_DIV14_GO_GENES = str_split(Down_DIV14_GOs$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Down_DIV14_GO_GENES)){
  Down_DIV14_GO_GENES[i,] = paste(substr(Down_DIV14_GO_GENES[i,], 1, 1),
                                tolower(substr(Down_DIV14_GO_GENES[i,], 2, 7)), sep = "") 
}


Up_DIV14_PATHWAY_GENES = str_split(Up_DIV14_PATHWAYS$Panther_2016.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Up_DIV14_PATHWAY_GENES)){
  Up_DIV14_PATHWAY_GENES[i,] = paste(substr(Up_DIV14_PATHWAY_GENES[i,], 1, 1),
                                  tolower(substr(Up_DIV14_PATHWAY_GENES[i,], 2, 7)), sep = "") 
}
Down_DIV14_PATHWAY_GENES = str_split(Down_DIV14_PATHWAYS$Panther_2016.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Down_DIV14_PATHWAY_GENES)){
  Down_DIV14_PATHWAY_GENES[i,] = paste(substr(Down_DIV14_PATHWAY_GENES[i,], 1, 1),
                                  tolower(substr(Down_DIV14_PATHWAY_GENES[i,], 2, 7)), sep = "") 
}


## ...In progress... Gene - Term associations..
#Genes_and_Terms = data.frame(1:nrow(DESeq2_Expression_Results3), 2)
#colnames(Genes_and_Terms) = c("GeneID", "Terms")
#Genes_and_Terms$GeneID = ""
#Genes_and_Terms$GeneID = as.character(DESeq2_Expression_Results3$GeneID)
#Genes_and_Terms$Terms = ""

## Export data as .CSVs to the server
#write.table(Adult_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\Adult GT DRG GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(Adult_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\Adult GT DRG Panther Pathways.csv", sep = ",", col.names = T)



##### In vitro Bar Plots of GO Terms Using ggplot2 #####
  ## Sort by -log10(Q). Retain the top 10 Processes
Up_DIV14_GOs$neglog10Q = -log10(Up_DIV14_GOs$GO_Biological_Process_2018.Adjusted.P.value)

Up_DIV14_GOs = Up_DIV14_GOs[order(Up_DIV14_GOs$neglog10Q, decreasing = TRUE),]

Up_GO_Bar_ADJP = Up_DIV14_GOs[1:10,]
  ## Remove the "GO...." strings in the GO Process column
Up_GO_Bar_ADJP$GO_Biological_Process_2018.Term = Up_GO_Bar_ADJP$GO_Biological_Process_2018.Term %>% gsub(x = Up_GO_Bar_ADJP$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")

  ## Plot the bar plot
  ## Store as a variable
  ## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long GO terms
Up_GO_BAR_ADJP = ggplot(Up_GO_Bar_ADJP, aes(x = reorder(`GO_Biological_Process_2018.Term`, `neglog10Q`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "2018 GO Biological Process", y = "Enrichr -log10(Adj. P-val)", title = "e13 in vitro Tmem184b GT/GT DRG GO Analysis of Up DEGs") + theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + 
  coord_flip()

  ## Plot it.
Up_GO_BAR_ADJP

  ## Repeat for Down DEGs
  ## Sort by -log10(Q). Retain the top 10 Processes
Down_DIV14_GOs$neglog10Q = -log10(Down_DIV14_GOs$GO_Biological_Process_2018.Adjusted.P.value)

Down_DIV14_GOs = Down_DIV14_GOs[order(Down_DIV14_GOs$neglog10Q, decreasing = TRUE),]

Down_GO_Bar_ADJP = Down_DIV14_GOs[1:10,]
  ## Remove the "GO...." strings in the GO Process column
Down_GO_Bar_ADJP$GO_Biological_Process_2018.Term = Down_GO_Bar_ADJP$GO_Biological_Process_2018.Term %>% gsub(x = Down_GO_Bar_ADJP$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")

  ## Plot the bar plot
  ## Store as a variable
  ## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long GO terms
Down_GO_BAR_ADJP = ggplot(Down_GO_Bar_ADJP, aes(x = reorder(`GO_Biological_Process_2018.Term`, `neglog10Q`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "2018 GO Biological Process", y = "Enrichr -log10(Adj. P-val)", title = "e13 in vitro Tmem184b GT/GT DRG GO Analysis of Down DEGs") + theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + 
  coord_flip()

  ## Plot it.
Down_GO_BAR_ADJP


##### In vitro Bar Plots of Panther Pathways Using ggplot2 #####
  ## Convert the Q-values to -log10 scale
Up_DIV14_PATHWAYS$neglog10Q = -log10(Up_DIV14_PATHWAYS$Panther_2016.Adjusted.P.value)
  ## Remove all but the top 15 pathways
Up_DIV14_PATHWAYS = Up_DIV14_PATHWAYS[order(Up_DIV14_PATHWAYS$neglog10Q, decreasing = TRUE),]

Up_DIV14_Pathway_Bar_ADJP = Up_DIV14_PATHWAYS[1:10,]
  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
Up_DIV14_Pathway_Bar_ADJP$Panther_2016.Term = Up_DIV14_Pathway_Bar_ADJP$Panther_2016.Term %>% gsub(x = Up_DIV14_Pathway_Bar_ADJP$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")


  ## Plot the bar plot
  ## Store as a variable
  ## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long pathway names
Up_DIV14_PATHWAY_BAR_ADJP = ggplot(Up_DIV14_Pathway_Bar_ADJP, aes(x = reorder(`Panther_2016.Term`, `neglog10Q`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "Panther 2016 Pathways", y = "Enrichr -log10(Adj. P-val)", title = "e13 in vitro Tmem184b GT/GT DRG Pathway Analysis of Up DEGs") + theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + 
  coord_flip()

  ## Plot it.
Up_DIV14_PATHWAY_BAR_ADJP




  ## Convert the Q-values to -log10 scale
Down_DIV14_PATHWAYS$neglog10Q = -log10(Down_DIV14_PATHWAYS$Panther_2016.Adjusted.P.value)
  ## Remove all but the top 15 pathways
Down_DIV14_PATHWAYS = Down_DIV14_PATHWAYS[order(Down_DIV14_PATHWAYS$neglog10Q, decreasing = TRUE),]

Down_DIV14_Pathway_Bar_ADJP = Down_DIV14_PATHWAYS[1:10,]
  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
Down_DIV14_Pathway_Bar_ADJP$Panther_2016.Term = Down_DIV14_Pathway_Bar_ADJP$Panther_2016.Term %>% gsub(x = Down_DIV14_Pathway_Bar_ADJP$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")


  ## Plot the bar plot
  ## Store as a variable
  ## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long pathway names
Down_DIV14_PATHWAY_BAR_ADJP = ggplot(Down_DIV14_Pathway_Bar_ADJP, aes(x = reorder(`Panther_2016.Term`, `neglog10Q`), y = `neglog10Q`)) + geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + labs(x = "Panther 2016 Pathways", y = "Enrichr -log10(Adj. P-val)", title = "e13 in vitro Tmem184b GT/GT DRG Pathway Analysis of Down DEGs") + theme_light() + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + coord_flip()

  ## Plot it.
Down_DIV14_PATHWAY_BAR_ADJP




## Export data as .CSVs to the server
#write.table(Up_DIV14_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 DIV14 GT DRG Up GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(Up_DIV14_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 DIV14 GT DRG Up Panther Pathways.csv", sep = ",", col.names = T)
#write.table(Down_DIV14_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 DIV14 GT DRG Down GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(Down_DIV14_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 DIV14 GT DRG Down Panther Pathways.csv", sep = ",", col.names = T)


##### In vitro GO Bio Processes and Panther Pathways Heatmaps with labels using ggplot2 #####
  ## Import the relevant data
DIV14_Up_DEGs_and_Pathways =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr Pathway Clustergram e13 DIV14 GT DRG Up.csv")
## Export the gene list for photoshopping
#write.csv(DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\aDRG FDR 05 Heatmap DEGs.csv", sep = ",")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
DIV14_Up_Pathway_Heatmap = melt(DIV14_Up_DEGs_and_Pathways, id = "DEGs")
colnames(DIV14_Up_Pathway_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

  ## Remove excessive pathway name strings
DIV14_Up_Pathway_Heatmap$Pathways = DIV14_Up_Pathway_Heatmap$Pathways %>% gsub(x = DIV14_Up_Pathway_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

  ## Make the digital values factors so it's compatible to graph
DIV14_Up_Pathway_Heatmap$`Presence In Pathway` = as.factor(DIV14_Up_Pathway_Heatmap$`Presence In Pathway`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

DIV14_Up_HEATlabs = ggplot(DIV14_Up_Pathway_Heatmap, aes(`Pathways`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Pathway`), color = "black") + scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.01)", title = "D.E.G.s of Affected Pathways in Tmem184b GT/GT e13 DIV14") + 
  guides(fill = guide_legend(title = "Presence in Pathway")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()

## Plot the heatmap
DIV14_Up_HEATlabs


  ## Import the relevant data
DIV14_Down_DEGs_and_Pathways =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr Pathway Clustergram e13 DIV14 GT DRG Down.csv")
## Export the gene list for photoshopping
#write.csv(DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\aDRG FDR 05 Heatmap DEGs.csv", sep = ",")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
DIV14_Down_Pathway_Heatmap = melt(DIV14_Down_DEGs_and_Pathways, id = "DEGs")
colnames(DIV14_Down_Pathway_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

  ## Remove excessive pathway name strings
DIV14_Down_Pathway_Heatmap$Pathways = DIV14_Down_Pathway_Heatmap$Pathways %>% gsub(x = DIV14_Down_Pathway_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

  ## Make the digital values factors so it's compatible to graph
DIV14_Down_Pathway_Heatmap$`Presence In Pathway` = as.factor(DIV14_Down_Pathway_Heatmap$`Presence In Pathway`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

DIV14_Down_HEATlabs = ggplot(DIV14_Down_Pathway_Heatmap, aes(`Pathways`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Pathway`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.01)", title = "D.E.G.s of Affected Pathways in Tmem184b GT/GT e13 DIV14") + 
  guides(fill = guide_legend(title = "Presence in Pathway")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()

  ## Plot the heatmap
DIV14_Down_HEATlabs



  ## Import the relevant data
DIV14_Up_DEGs_and_GOs =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr GO Clustergram e13 DIV14 GT DRG Up.csv")
## Export the gene list for photoshopping
#write.csv(DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\aDRG FDR 05 Heatmap DEGs.csv", sep = ",")


  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
DIV14_Up_GO_Heatmap = melt(DIV14_Up_DEGs_and_GOs, id = "DEGs")
colnames(DIV14_Up_GO_Heatmap) = c("DEGs", "GOs", "Presence In Process")

  ## Remove the "GO...." strings in the GO Process column
DIV14_Up_GO_Heatmap$GOs = DIV14_Up_GO_Heatmap$GOs %>% gsub(x = DIV14_Up_GO_Heatmap$GOs, pattern = " \\(.*\\)$", replacement = "")

  ## Make the digital values factors so it's compatible to graph
DIV14_Up_GO_Heatmap$`Presence In Process` = as.factor(DIV14_Up_GO_Heatmap$`Presence In Process`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

DIV14_Up_HEATlabs = ggplot(DIV14_Up_GO_Heatmap, aes(`GOs`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Process`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "2018 GO Biological Process", y = "Differentially Expressed Genes (FDR < 0.05)", title = "D.E.G.s of Affected GO Processes in Tmem184b GT/GT e13 DIV14") + 
  guides(fill = guide_legend(title = "Presence in Process")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


  ## Plot the heatmap
DIV14_Up_HEATlabs



  ## Import the relevant data
DIV14_Down_DEGs_and_GOs =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr GO Clustergram e13 DIV14 GT DRG Down.csv")
## Export the gene list for photoshopping
#write.csv(DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\aDRG FDR 05 Heatmap DEGs.csv", sep = ",")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
DIV14_Down_GO_Heatmap = melt(DIV14_Down_DEGs_and_GOs, id = "DEGs")
colnames(DIV14_Down_GO_Heatmap) = c("DEGs", "GOs", "Presence In Process")

  ## Remove the "GO...." strings in the GO Process column
DIV14_Down_GO_Heatmap$GOs = DIV14_Down_GO_Heatmap$GOs %>% gsub(x = DIV14_Down_GO_Heatmap$GOs, pattern = " \\(.*\\)$", replacement = "")

  ## Make the digital values factors so it's compatible to graph
DIV14_Down_GO_Heatmap$`Presence In Process` = as.factor(DIV14_Down_GO_Heatmap$`Presence In Process`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

DIV14_Down_HEATlabs = ggplot(DIV14_Down_GO_Heatmap, aes(`GOs`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Process`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "2018 GO Biological Process", y = "Differentially Expressed Genes (FDR < 0.05)", title = "D.E.G.s of Affected GO Processes in Tmem184b GT/GT e13 DIV14") + 
  guides(fill = guide_legend(title = "Presence in Process")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


  ## Plot the heatmap
DIV14_Down_HEATlabs


##### In vitro DRG RESCUE Data Prep #####
  ## Import the dataset you want to analyze
DESeq2_In_vitro_Rescue = read.csv("M:/Erik/Data/Omics/TimeCourse/Processed Galaxy Output/Test Results to Upload/Cultured Mutant Embryos with Rescue.csv")
  ## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
DESeq2_In_vitro_Rescue3 = subset(DESeq2_In_vitro_Rescue, (!is.na(DESeq2_In_vitro_Rescue[,"AdjP"])))

  ## Create a column in the DESeq2 dataframe that scales the Adjusted P-value by log-base 10
DESeq2_In_vitro_Rescue3$log10Pval = -log10(DESeq2_In_vitro_Rescue3$AdjP)

  ## Find the mean of the log2(FC) values to set the threshold for downstream analysis
#m = subset(DESeq2_In_vitro_Rescue3, DESeq2_In_vitro_Rescue3$AdjP < 0.01 & DESeq2_In_vitro_Rescue3$log2.FC. > 0)
  ## Copy the value to paste downstream
#mean(m$log2.FC.)

  ## Create the data points of interest for functional display.. ggplot struggles to plot subsets of dataframes, so creating a variable that directly labels all the points accordingly helps.
  ## Create a new column (genes of interest); Using "GeneID" is irrelevant; the values will be replaced in the next step
DESeq2_In_vitro_Rescue3$g.o.i. = DESeq2_In_vitro_Rescue3$GeneID

  ## Fill the points with appropriately indexed data
DESeq2_In_vitro_Rescue3$g.o.i.[which(c(DESeq2_In_vitro_Rescue3$AdjP < 0.01 & DESeq2_In_vitro_Rescue3$log2.FC. > 0))]= "Up DEGs"
DESeq2_In_vitro_Rescue3$g.o.i.[which(c(DESeq2_In_vitro_Rescue3$AdjP < 0.01 & DESeq2_In_vitro_Rescue3$log2.FC. < 0))]= "Down DEGs"
#DESeq2_In_vitro_Rescue3$g.o.i.[which(c(DESeq2_In_vitro_Rescue3$AdjP < 0.01 & (DESeq2_In_vitro_Rescue3$log2.FC. > -1.856636 & DESeq2_In_vitro_Rescue3$log2.FC. < 1.856636)))] = "Stable DEGs"
DESeq2_In_vitro_Rescue3$g.o.i.[which(c(DESeq2_In_vitro_Rescue3$AdjP > 0.01))] = "Non-DEGs"


##### Volcano, MA Plots of In vitro GT DRG Rescue #####
  ## Notes on ggplot commands:
    ## geom-point = alpha is a ggplot add-on function that controls point color transparency
    ## coord_cartesian is a ggplot add-on function that determines axes sizes
    ## labs is a ggplot add-on function that determines the plot labels; x-axis, y-axis, title, legend colors
    ## scale_color_manual is a ggplot add-on function that sets the colors of the plot, corresponding to the groups determined in the the base function ("col" in aes)
    ## Theme is used to manipulate the physical location of the plot title

    ## Pay attention to the exact strings of column names.
volrescue = ggplot(DESeq2_In_vitro_Rescue3) +
  geom_point(data = subset(DESeq2_In_vitro_Rescue3, `g.o.i.` == "Non-DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.2) +
  geom_point(data = subset(DESeq2_In_vitro_Rescue3, `g.o.i.` == "Up DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.3) +
  geom_point(data = subset(DESeq2_In_vitro_Rescue3, `g.o.i.` == "Down DEGs"),
             aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.3) +
  #geom_point(data = subset(DESeq2_In_vitro_Rescue3, `g.o.i.` == "Stable DEGs"),
   #          aes(x = `log2.FC.`, y = `log10Pval`, color = `g.o.i.`), alpha = 0.8) +
  
  coord_cartesian(xlim = c(-10,10), ylim = c(0,200)) +
  labs(x = "log2(FC)", y = "-log10(Adj.P-val)", title = "In vitro e13 Tmem184b GT DRG Rescue Gene Exp. Changes", col = "Gene Data Type") +
  scale_color_manual(values = c("navy", "gray48", "darkgoldenrod4")) +
  #geom_hline(yintercept = 2.0, type = "solid", color = "black") +
  #geom_vline(xintercept = -1.856636, linetype = "dashed", color = "black") +
  #geom_vline(xintercept = 1.856636, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none")

  ## Plot it
volrescue

  ## Transform the mean normalized counts
DESeq2_In_vitro_Rescue3$BaseMean = log(DESeq2_In_vitro_Rescue3$BaseMean)
  ## Desgin the MA plot
MA = ggplot(data = DESeq2_In_vitro_Rescue3) +
  geom_point(data = subset(DESeq2_In_vitro_Rescue3, `g.o.i.` == "Non-DEGs"),
             aes(x = `BaseMean`, y = `log2.FC.`, col = `g.o.i.`), alpha = 0.2) +
  geom_point(data = subset(DESeq2_In_vitro_Rescue3, `g.o.i.` == "Up DEGs"),
             aes(x = `BaseMean`, y = `log2.FC.`, col = `g.o.i.`), alpha = 0.3) +
  geom_point(data = subset(DESeq2_In_vitro_Rescue3, `g.o.i.` == "Down DEGs"),
             aes(x = `BaseMean`, y = `log2.FC.`, col = `g.o.i.`), alpha = 0.3) +
  coord_cartesian(xlim = c(0,13), ylim = c(-10.0,10.0)) + 
  labs(x = "ln(Mean of Normalized Counts)", y = "log2(FC)", title = "e13 in vitro Tmem184b GT/GT DRG Rescue Gene Expression Changes", col = "Gene Data Type") + 
  scale_color_manual(values = c("navy", "gray48", "darkgoldenrod4")) + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

## Plot it
MA


##### Data prep for Rescue bar plots and heatmaps #####
## Perform relevant analysis through Enrichr; in this case, FDR < 0.01
  ## first: GO Bio Processes; this analyzes genes from your dataset in Enrichr, which draws terms from geneontology.org
Up_Rescue_GOs = as.data.frame(
  enrichr(
    c(DESeq2_In_vitro_Rescue3$GeneID[which(DESeq2_In_vitro_Rescue3$g.o.i. == "Up DEGs")]), DBs[1])
)
Down_Rescue_GOs = as.data.frame(
  enrichr(
    c(DESeq2_In_vitro_Rescue3$GeneID[which(DESeq2_In_vitro_Rescue3$g.o.i. == "Down DEGs")]), DBs[1])
)

  ## second: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
Up_Rescue_PATHWAYS = as.data.frame(
  enrichr(
    c(DESeq2_In_vitro_Rescue3$GeneID[which(DESeq2_In_vitro_Rescue3$g.o.i == "Up DEGs")]), DBs[2])
)
Down_Rescue_PATHWAYS = as.data.frame(
  enrichr(
    c(DESeq2_In_vitro_Rescue3$GeneID[which(DESeq2_In_vitro_Rescue3$g.o.i == "Down DEGs")]), DBs[2])
)


  ## Remove unnecessary columns
Up_Rescue_GOs = Up_Rescue_GOs[,-c(5,6,7)]
Down_Rescue_GOs = Down_Rescue_GOs[,-c(5,6,7)]
Up_Rescue_PATHWAYS = Up_Rescue_PATHWAYS[,-c(5,6,7)]
Down_Rescue_PATHWAYS = Down_Rescue_PATHWAYS[,-c(5,6,7)]


  ## Clean up both dataframes and export both dataframes.
  ## Dataframes need additional cleaning in Python;
  ## Following that manipulation, data will be re-imported in other scripts for visualization using ggplot
Up_Rescue_GOs = Up_Rescue_GOs %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Up_Rescue_GOs$Overlap = as.numeric(Up_Rescue_GOs$Overlap)
Up_Rescue_GOs$`Process Size` = as.numeric(Up_Rescue_GOs$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Up_Rescue_GOs = add_column(Up_Rescue_GOs, EnrichrZscore = Up_Rescue_GOs$GO_Biological_Process_2018.Combined.Score/log(Up_Rescue_GOs$GO_Biological_Process_2018.P.value), .before = 6)

  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Up_Rescue_GOs = add_column(Up_Rescue_GOs, Weighted_Overlap_Ratio = Up_Rescue_GOs$Overlap*(Up_Rescue_GOs$Overlap/Up_Rescue_GOs$`Process Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
Up_Rescue_GOs = add_column(Up_Rescue_GOs, Modified_Combined_Score = (abs(Up_Rescue_GOs$EnrichrZscore))*(-log(Up_Rescue_GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)


Down_Rescue_GOs = Down_Rescue_GOs %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Down_Rescue_GOs$Overlap = as.numeric(Down_Rescue_GOs$Overlap)
Down_Rescue_GOs$`Process Size` = as.numeric(Down_Rescue_GOs$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Down_Rescue_GOs = add_column(Down_Rescue_GOs, EnrichrZscore = Down_Rescue_GOs$GO_Biological_Process_2018.Combined.Score/log(Down_Rescue_GOs$GO_Biological_Process_2018.P.value), .before = 6)

  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Down_Rescue_GOs = add_column(Down_Rescue_GOs, Weighted_Overlap_Ratio = Down_Rescue_GOs$Overlap*(Down_Rescue_GOs$Overlap/Down_Rescue_GOs$`Process Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
Down_Rescue_GOs = add_column(Down_Rescue_GOs, Modified_Combined_Score = (abs(Down_Rescue_GOs$EnrichrZscore))*(-log(Down_Rescue_GOs$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)




  ## Repeat the above code for Panther Pathways analyzed from up DEGs
Up_Rescue_PATHWAYS = Up_Rescue_PATHWAYS %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Up_Rescue_PATHWAYS$Overlap = as.numeric(Up_Rescue_PATHWAYS$Overlap)
Up_Rescue_PATHWAYS$`Process Size` = as.numeric(Up_Rescue_PATHWAYS$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Up_Rescue_PATHWAYS = add_column(Up_Rescue_PATHWAYS, EnrichrZscore = Up_Rescue_PATHWAYS$Panther_2016.Combined.Score/log(Up_Rescue_PATHWAYS$Panther_2016.P.value), .before = 6)

  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Up_Rescue_PATHWAYS = add_column(Up_Rescue_PATHWAYS, Weighted_Overlap_Ratio = Up_Rescue_PATHWAYS$Overlap*(Up_Rescue_PATHWAYS$Overlap/Up_Rescue_PATHWAYS$`Process Size`), .before = 4)

  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
Up_Rescue_PATHWAYS = add_column(Up_Rescue_PATHWAYS, Modified_Combined_Score = (abs(Up_Rescue_PATHWAYS$EnrichrZscore))*(-log(Up_Rescue_PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)


  ## Repeat the above code for Panther Pathways analyzed from down DEGs
Down_Rescue_PATHWAYS = Down_Rescue_PATHWAYS %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Down_Rescue_PATHWAYS$Overlap = as.numeric(Down_Rescue_PATHWAYS$Overlap)
Down_Rescue_PATHWAYS$`Process Size` = as.numeric(Down_Rescue_PATHWAYS$`Process Size`)

  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across all terms/pathways when compared to the term/pathway's "pre-rank" based on multiple simulations of the same number of enriched genes.
Down_Rescue_PATHWAYS = add_column(Down_Rescue_PATHWAYS, EnrichrZscore = Down_Rescue_PATHWAYS$Panther_2016.Combined.Score/log(Down_Rescue_PATHWAYS$Panther_2016.P.value), .before = 6)

  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway ("overlap"; numerator) out of the total genes in that term/pathway (denominator) is scaled by the overlap. This weights differentially expressed genes by pathway size.
Down_Rescue_PATHWAYS = add_column(Down_Rescue_PATHWAYS, Weighted_Overlap_Ratio = Down_Rescue_PATHWAYS$Overlap*(Down_Rescue_PATHWAYS$Overlap/Down_Rescue_PATHWAYS$`Process Size`), .before = 4)

  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the "Weighted Overlap Ratio" and Enrichr Z-score. 
Down_Rescue_PATHWAYS = add_column(Down_Rescue_PATHWAYS, Modified_Combined_Score = (abs(Down_Rescue_PATHWAYS$EnrichrZscore))*(-log(Down_Rescue_PATHWAYS$Panther_2016.Adjusted.P.value)), .before = 6)



Up_Rescue_GO_GENES = str_split(Up_Rescue_GOs$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Up_Rescue_GO_GENES)){
  Up_Rescue_GO_GENES[i,] = paste(substr(Up_Rescue_GO_GENES[i,], 1, 1),
                                tolower(substr(Up_Rescue_GO_GENES[i,], 2, 7)), sep = "") 
}
Down_Rescue_GO_GENES = str_split(Down_Rescue_GOs$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Down_Rescue_GO_GENES)){
  Down_Rescue_GO_GENES[i,] = paste(substr(Down_Rescue_GO_GENES[i,], 1, 1),
                                  tolower(substr(Down_Rescue_GO_GENES[i,], 2, 7)), sep = "") 
}



Up_Rescue_PATHWAY_GENES = str_split(Up_Rescue_PATHWAYS$Panther_2016.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Up_Rescue_PATHWAY_GENES)){
  Up_Rescue_PATHWAY_GENES[i,] = paste(substr(Up_Rescue_PATHWAY_GENES[i,], 1, 1),
                                     tolower(substr(Up_Rescue_PATHWAY_GENES[i,], 2, 7)), sep = "") 
}
Down_Rescue_PATHWAY_GENES = str_split(Down_Rescue_PATHWAYS$Panther_2016.Genes, pattern = ";", simplify = TRUE)
for (i in 1:nrow(Down_Rescue_PATHWAY_GENES)){
  Down_Rescue_PATHWAY_GENES[i,] = paste(substr(Down_Rescue_PATHWAY_GENES[i,], 1, 1),
                                       tolower(substr(Down_Rescue_PATHWAY_GENES[i,], 2, 7)), sep = "") 
}


##### In vitro Rescue Bar Plots of GO Terms Using ggplot2 #####
  ## Sort by -log10(Q). Retain the top 10 Processes
Up_Rescue_GOs$neglog10Q = -log10(Up_Rescue_GOs$GO_Biological_Process_2018.Adjusted.P.value)

Up_Rescue_GOs = Up_Rescue_GOs[order(Up_Rescue_GOs$neglog10Q, decreasing = TRUE),]

Up_Rescue_GO_Bar_ADJP = Up_Rescue_GOs[1:10,]
  ## Remove the "GO...." strings in the GO Process column
Up_Rescue_GO_Bar_ADJP$GO_Biological_Process_2018.Term = Up_Rescue_GO_Bar_ADJP$GO_Biological_Process_2018.Term %>% gsub(x = Up_Rescue_GO_Bar_ADJP$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")

  ## Plot the bar plot
  ## Store as a variable
  ## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long GO terms
Up_Rescue_GO_BAR_ADJP = ggplot(Up_Rescue_GO_Bar_ADJP, aes(x = reorder(`GO_Biological_Process_2018.Term`, `neglog10Q`), y = `neglog10Q`)) +
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "2018 GO Biological Process", y = "Enrichr -log10(Adj. P-val)", title = "e13 in vitro Tmem184b GT/GT DRG GO Analysis of Up DEGs") +
  theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + coord_flip()

  ## Plot it.
Up_Rescue_GO_BAR_ADJP

## Repeat for Down DEGs
  ## Sort by -log10(Q). Retain the top 10 Processes
Down_Rescue_GOs$neglog10Q = -log10(Down_Rescue_GOs$GO_Biological_Process_2018.Adjusted.P.value)

Down_Rescue_GOs = Down_Rescue_GOs[order(Down_Rescue_GOs$neglog10Q, decreasing = TRUE),]

Down_Rescue_GO_Bar_ADJP = Down_Rescue_GOs[1:10,]
  ## Remove the "GO...." strings in the GO Process column
Down_Rescue_GO_Bar_ADJP$GO_Biological_Process_2018.Term = Down_Rescue_GO_Bar_ADJP$GO_Biological_Process_2018.Term %>% gsub(x = Down_Rescue_GO_Bar_ADJP$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")

  ## Plot the bar plot
  ## Store as a variable
  ## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long GO terms
Down_Rescue_GO_BAR_ADJP = ggplot(Down_Rescue_GO_Bar_ADJP, aes(x = reorder(`GO_Biological_Process_2018.Term`, `neglog10Q`), y = `neglog10Q`)) +
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "2018 GO Biological Process", y = "Enrichr -log10(Adj. P-val)", title = "e13 in vitro Tmem184b GT/GT DRG GO Analysis of Down DEGs") +
  theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + 
  coord_flip()

  ## Plot it.
Down_Rescue_GO_BAR_ADJP


##### In vitro Rescue Bar Plots of Panther Pathways Using ggplot2 #####
  ## Convert the Q-values to -log10 scale
Up_Rescue_PATHWAYS$neglog10Q = -log10(Up_Rescue_PATHWAYS$Panther_2016.Adjusted.P.value)
  ## Remove all but the top 15 pathways
Up_Rescue_PATHWAYS = Up_Rescue_PATHWAYS[order(Up_Rescue_PATHWAYS$neglog10Q, decreasing = TRUE),]

Up_Rescue_Pathway_Bar_ADJP = Up_Rescue_PATHWAYS[1:10,]
  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
Up_Rescue_Pathway_Bar_ADJP$Panther_2016.Term = Up_Rescue_Pathway_Bar_ADJP$Panther_2016.Term %>% gsub(x = Up_Rescue_Pathway_Bar_ADJP$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")


  ## Plot the bar plot
    ## Store as a variable
    ## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
    ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
    ## scale_fill_continuous() creates the color gradient and legend
    ## coord_flip() switches axes for better visualization of the long pathway names
Up_Rescue_PATHWAY_BAR_ADJP = ggplot(Up_Rescue_Pathway_Bar_ADJP, aes(x = reorder(`Panther_2016.Term`, `neglog10Q`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "Panther 2016 Pathways", y = "Enrichr -log10(Adj. P-val)", title = "e13 in vitro Tmem184b GT/GT DRG Pathway Analysis of Up DEGs") +
  theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + 
  coord_flip()

## Plot it.
Up_Rescue_PATHWAY_BAR_ADJP




  ## Convert the Q-values to -log10 scale
Down_Rescue_PATHWAYS$neglog10Q = -log10(Down_Rescue_PATHWAYS$Panther_2016.Adjusted.P.value)
  ## Remove all but the top 15 pathways
Down_Rescue_PATHWAYS = Down_Rescue_PATHWAYS[order(Down_Rescue_PATHWAYS$neglog10Q, decreasing = TRUE),]

Down_Rescue_Pathway_Bar_ADJP = Down_Rescue_PATHWAYS[1:10,]
  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
Down_Rescue_Pathway_Bar_ADJP$Panther_2016.Term = Down_Rescue_Pathway_Bar_ADJP$Panther_2016.Term %>% gsub(x = Down_Rescue_Pathway_Bar_ADJP$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")


  ## Plot the bar plot
    ## Store as a variable
    ## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
    ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
    ## scale_fill_continuous() creates the color gradient and legend
    ## coord_flip() switches axes for better visualization of the long pathway names
Down_Rescue_PATHWAY_BAR_ADJP = ggplot(Down_Rescue_Pathway_Bar_ADJP, aes(x = reorder(`Panther_2016.Term`, `neglog10Q`), y = `neglog10Q`)) +
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "Panther 2016 Pathways", y = "Enrichr -log10(Adj. P-val)", title = "e13 in vitro Tmem184b GT/GT DRG Pathway Analysis of Down DEGs") +
  theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + 
  coord_flip()

## Plot it.
Down_Rescue_PATHWAY_BAR_ADJP



## Export data as .CSVs to the server
#write.table(Up_Rescue_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 Rescue GT DRG Up GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(Up_Rescue_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 Rescue GT DRG Up Panther Pathways.csv", sep = ",", col.names = T)
#write.table(Down_Rescue_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 Rescue GT DRG Down GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(Down_Rescue_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 Rescue GT DRG Down Panther Pathways.csv", sep = ",", col.names = T)


##### Rescue GO Bio Processes and Panther Pathways Heatmaps with labels using ggplot2 #####
  ## Import the relevant data
Rescue_Up_DEGs_and_Pathways =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr Pathway Clustergram Rescue GT DRG Up.csv")
  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
Rescue_Up_Pathway_Heatmap = melt(Rescue_Up_DEGs_and_Pathways, id = "DEGs")
colnames(Rescue_Up_Pathway_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

  ## Remove excessive pathway name strings
Rescue_Up_Pathway_Heatmap$Pathways = Rescue_Up_Pathway_Heatmap$Pathways %>% gsub(x = Rescue_Up_Pathway_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

  ## Make the digital values factors so it's compatible to graph
Rescue_Up_Pathway_Heatmap$`Presence In Pathway` = as.factor(Rescue_Up_Pathway_Heatmap$`Presence In Pathway`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

Rescue_Up_HEATlabs = ggplot(Rescue_Up_Pathway_Heatmap, aes(`Pathways`, `DEGs`)) + geom_tile(aes(fill = `Presence In Pathway`), color = "black") + scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.01)", title = "D.E.G.s of Affected Pathways in Tmem184b GT/GT e13 Rescue") + 
  guides(fill = guide_legend(title = "Presence in Pathway")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()

  ## Plot the heatmap
Rescue_Up_HEATlabs


  ## Import the relevant data
Rescue_Down_DEGs_and_Pathways =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr Pathway Clustergram Rescue GT DRG Down.csv")
  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
Rescue_Down_Pathway_Heatmap = melt(Rescue_Down_DEGs_and_Pathways, id = "DEGs")
colnames(Rescue_Down_Pathway_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

  ## Remove excessive pathway name strings
Rescue_Down_Pathway_Heatmap$Pathways = Rescue_Down_Pathway_Heatmap$Pathways %>% gsub(x = Rescue_Down_Pathway_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

  ## Make the digital values factors so it's compatible to graph
Rescue_Down_Pathway_Heatmap$`Presence In Pathway` = as.factor(Rescue_Down_Pathway_Heatmap$`Presence In Pathway`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

Rescue_Down_HEATlabs = ggplot(Rescue_Down_Pathway_Heatmap, aes(`Pathways`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Pathway`), color = "black") + scale_fill_manual(values = c("navy", "gold3")) + labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.01)", title = "D.E.G.s of Affected Pathways in Tmem184b GT/GT e13 Rescue") + 
  guides(fill = guide_legend(title = "Presence in Pathway")) + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()

## Plot the heatmap
Rescue_Down_HEATlabs



  ## Import the relevant data
Rescue_Up_DEGs_and_GOs =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr GO Clustergram Rescue GT DRG Up.csv")
  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
Rescue_Up_GO_Heatmap = melt(Rescue_Up_DEGs_and_GOs, id = "DEGs")
colnames(Rescue_Up_GO_Heatmap) = c("DEGs", "GOs", "Presence In Process")

  ## Remove the "GO...." strings in the GO Process column
Rescue_Up_GO_Heatmap$GOs = Rescue_Up_GO_Heatmap$GOs %>% gsub(x = Rescue_Up_GO_Heatmap$GOs, pattern = " \\(.*\\)$", replacement = "")

  ## Make the digital values factors so it's compatible to graph
Rescue_Up_GO_Heatmap$`Presence In Process` = as.factor(Rescue_Up_GO_Heatmap$`Presence In Process`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

Rescue_Up_HEATlabs = ggplot(Rescue_Up_GO_Heatmap, aes(`GOs`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Process`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "2018 GO Biological Process", y = "Differentially Expressed Genes (FDR < 0.05)", title = "D.E.G.s of Affected GO Processes in Tmem184b GT/GT e13 Rescue") + 
  guides(fill = guide_legend(title = "Presence in Process")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


  ## Plot the heatmap
Rescue_Up_HEATlabs



  ## Import the relevant data
Rescue_Down_DEGs_and_GOs =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr GO Clustergram Rescue GT DRG Down.csv")
  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
Rescue_Down_GO_Heatmap = melt(Rescue_Down_DEGs_and_GOs, id = "DEGs")
colnames(Rescue_Down_GO_Heatmap) = c("DEGs", "GOs", "Presence In Process")

  ## Remove the "GO...." strings in the GO Process column
Rescue_Down_GO_Heatmap$GOs = Rescue_Down_GO_Heatmap$GOs %>% gsub(x = Rescue_Down_GO_Heatmap$GOs, pattern = " \\(.*\\)$", replacement = "")

  ## Make the digital values factors so it's compatible to graph
Rescue_Down_GO_Heatmap$`Presence In Process` = as.factor(Rescue_Down_GO_Heatmap$`Presence In Process`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

Rescue_Down_HEATlabs = ggplot(Rescue_Down_GO_Heatmap, aes(`GOs`, `DEGs`)) + geom_tile(aes(fill = `Presence In Process`), color = "black") +
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "2018 GO Biological Process", y = "Differentially Expressed Genes (FDR < 0.05)", title = "D.E.G.s of Affected GO Processes in Tmem184b GT/GT e13 Rescue") + 
  guides(fill = guide_legend(title = "Presence in Process")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()


  ## Plot the heatmap
Rescue_Down_HEATlabs

  ## Export data as .CSVs to the server
#write.table(Up_Rescue_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 Rescue GT DRG Up GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(Up_Rescue_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 Rescue GT DRG Up Panther Pathways.csv", sep = ",", col.names = T)
#write.table(Down_Rescue_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 Rescue GT DRG Down GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(Down_Rescue_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\e13 Rescue GT DRG Down Panther Pathways.csv", sep = ",", col.names = T)


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



##### PANTHER IX Bar plots #####
  ## Import the intersection files generated 
Panther_GO_Bio_Comp_Bar = read_csv("M:/PAPER ASSEMBLY/Itch paper/Molecular Mechanism/eDRGs/Panther GO Biological Processes Complete Overrepresented from genes directly affected by Tmem Trimmed for graph.csv")

  ## Re-order the terms by FDR (Q value aka Adjusted P-value)
Panther_GO_Bio_Comp_Bar = Panther_GO_Bio_Comp_Bar[order(Panther_GO_Bio_Comp_Bar$FDR, decreasing = FALSE),]
  ## Remove unnecessary rows, columns
Panther_GO_Bio_Comp_Bar = Panther_GO_Bio_Comp_Bar[-15,-c(2,3,4,6)]
rownames(Panther_GO_Bio_Comp_Bar) = NULL

  ## Transform the Q
Panther_GO_Bio_Comp_Bar$neglog10 = -log10(Panther_GO_Bio_Comp_Bar$FDR)

  ## Remove the ")" strings in the GO Process column
Panther_GO_Bio_Comp_Bar$`GO biological process complete` = Panther_GO_Bio_Comp_Bar$`GO biological process complete` %>% gsub(x = Panther_GO_Bio_Comp_Bar$`GO biological process complete`, pattern = "^.*\\)", replacement = "")

  ## Take the top 10 of the data
Panther_GO_Bio_Comp_Bar = Panther_GO_Bio_Comp_Bar[1:10,]

  ## Create the bar plot
Panther_GO_Bio_Comp_BAR = ggplot(Panther_GO_Bio_Comp_Bar, aes(x = reorder(`GO biological process complete`, `neglog10`), y = `neglog10`)) +
  geom_col(stat = "identity" , aes(fill = `neglog10`)) +
  scale_fill_continuous(name = "-log10(Adj. P-val)") +
  labs(x = "Panther 2016 GO Biological Process Complete", y = "-log10(Adj. P-val)", title = "Pathways Containing In vitro DEGs Both Down in Mutant and Up upon Rescue") +
  theme_classic() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), legend.position = "none") +
  coord_flip()

  ## Plot it
Panther_GO_Bio_Comp_BAR




##### Enrichr IX supp Bar plots #####
  ## GOs
  ## Sort by -log10(Q). Retain the top 10 Processes 
IX_supp_GOs$neglog10Q = -log10(IX_supp_GOs$GO_Biological_Process_2018.Adjusted.P.value)
IX_supp_GOs = IX_supp_GOs[order(IX_supp_GOs$neglog10Q, decreasing = TRUE),]
IX_supp_GOs_Bar_ADJP = IX_supp_GOs[1:10,]
  ## Remove the "GO...." strings in the GO Process column
IX_supp_GOs_Bar_ADJP$GO_Biological_Process_2018.Term = IX_supp_GOs_Bar_ADJP$GO_Biological_Process_2018.Term %>% gsub(x = IX_supp_GOs_Bar_ADJP$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")

  ## Re-order the indeces
rownames(IX_supp_GOs_Bar_ADJP) = NULL

  ## Shorten the term strings
IX_supp_GOs_Bar_ADJP$GO_Biological_Process_2018.Term[2] = "secondary alcohol\nbiosynthetic process"
IX_supp_GOs_Bar_ADJP$GO_Biological_Process_2018.Term[5] = "fatty acid\ntransmembrane transport"
IX_supp_GOs_Bar_ADJP$GO_Biological_Process_2018.Term[7] = "regulation of alcohol\nbiosynthetic process"
IX_supp_GOs_Bar_ADJP$GO_Biological_Process_2018.Term[8] = "protein localization\nto cell periphery"
IX_supp_GOs_Bar_ADJP$GO_Biological_Process_2018.Term[10] = "protein\nhomooligomerization"

  ## Plot the bar plot
  ## Store as a variable
  ## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long GO terms
IX_supp_GOs_BAR_ADJP = ggplot(IX_supp_GOs_Bar_ADJP, aes(x = reorder(`GO_Biological_Process_2018.Term`, `neglog10Q`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "2018 GO Biological Process", y = "Enrichr -log10(Adj. P-val)", title = "GO Analysis of DEGs Affected by Tmem184b in vitro") +
  theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_flip()

  ## Plot it.
IX_supp_GOs_BAR_ADJP


  ## Pathways
  ## Sort by -log10(Q). Retain the top 10 Processes 
IX_supp_PATHWAYS$neglog10Q = -log10(IX_supp_PATHWAYS$Panther_2016.Adjusted.P.value)
IX_supp_PATHWAYS = IX_supp_PATHWAYS[order(IX_supp_PATHWAYS$neglog10Q, decreasing = TRUE),]
IX_supp_Pathway_Bar_ADJP = IX_supp_PATHWAYS[1:10,]
  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
IX_supp_Pathway_Bar_ADJP$Panther_2016.Term = IX_supp_Pathway_Bar_ADJP$Panther_2016.Term %>% gsub(x = IX_supp_Pathway_Bar_ADJP$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")

  ## re-order the indeces
rownames(IX_supp_Pathway_Bar_ADJP) = NULL

  ## Shorten the term strings
IX_supp_Pathway_Bar_ADJP$Panther_2016.Term[2] = "Alpha adrenergic\nreceptor signaling pathway"
IX_supp_Pathway_Bar_ADJP$Panther_2016.Term[3] = "5HT2 type receptor\nmediated signaling pathway"
IX_supp_Pathway_Bar_ADJP$Panther_2016.Term[4] = "Fructose\ngalactose metabolism"
IX_supp_Pathway_Bar_ADJP$Panther_2016.Term[6] = "Histamine H1 receptor\nmediated signaling pathway"
IX_supp_Pathway_Bar_ADJP$Panther_2016.Term[7] = "Nicotine\npharmacodynamics pathway"
IX_supp_Pathway_Bar_ADJP$Panther_2016.Term[8] = "Heterotrimeric G-protein\nsignaling pathway-\nGq alpha and\nGo alpha mediated pathway"
IX_supp_Pathway_Bar_ADJP$Panther_2016.Term[9] = "Nicotinic acetylcholine\nreceptor signaling pathway"
IX_supp_Pathway_Bar_ADJP$Panther_2016.Term[10] = "2-arachidonoylglycerol\nbiosynthesis"

  ## Plot the bar plot
  ## Store as a variable
  ## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
  ## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
  ## scale_fill_continuous() creates the color gradient and legend
  ## coord_flip() switches axes for better visualization of the long pathway names
IX_supp_PATHWAY_BAR_ADJP = ggplot(IX_supp_Pathway_Bar_ADJP, aes(x = reorder(`Panther_2016.Term`, `neglog10Q`), y = `neglog10Q`)) + 
  geom_col(stat = "identity" , aes(fill = `neglog10Q`)) + 
  scale_fill_continuous(name = "Enrichr \n-log10(Adj. P-val)") + 
  labs(x = "Panther 2016 Pathways", y = "Enrichr -log10(Adj. P-val)", title = "Pathway Analysis of DEGs Affected by Tmem184b in vitro") +
  theme_bw() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_flip()

  ## Plot it.
IX_supp_PATHWAY_BAR_ADJP


  ## Export data as .CSVs to the server
#write.table(IX_supp_GOs, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\IX supp DIV14 GT DRG GO Bio Processes 0 LFC.csv", sep = ",", col.names = T)
#write.table(IX_supp_PATHWAYS, "M:\\Erik\\Data\\Omics\\RNAseq\\Consensus Raw\\IX supp DIV14 GT DRG Panther Pathways 0 LFC.csv", sep = ",", col.names = T)


##### Intersection Rescue (opposite) GO Bio Processes and Panther Pathways Heatmaps with labels using ggplot2 #####
  ## Import the relevant data
IX_supp_Pathways =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr Pathway Clustergram IX supp.csv")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
IX_supp_Pathway_Heatmap = melt(IX_supp_Pathways, id = "DEGs")
colnames(IX_supp_Pathway_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

# # Remove excessive pathway name strings
IX_supp_Pathway_Heatmap$Pathways = IX_supp_Pathway_Heatmap$Pathways %>% gsub(x = IX_supp_Pathway_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

  ## Make the digital (0,1) values factors so it's compatible to graph
IX_supp_Pathway_Heatmap$`Presence In Pathway` = as.factor(IX_supp_Pathway_Heatmap$`Presence In Pathway`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

IX_supp_HEATlabs = ggplot(IX_supp_Pathway_Heatmap, aes(`Pathways`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Pathway`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.01)", title = "Pathway Analysis of DEGs Affected by Tmem184b in vitro") + 
  guides(fill = guide_legend(title = "Presence in Pathway")) + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()

  ## Plot the heatmap
IX_supp_HEATlabs


  ## Import the relevant data
IX_supp_GOs =  read_csv("M:/Erik/Data/Omics/RNAseq/Custom Python Enrichr GO Clustergram IX supp.csv")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
IX_supp_GO_Heatmap = melt(IX_supp_GOs, id = "DEGs")
colnames(IX_supp_GO_Heatmap) = c("DEGs", "GOs", "Presence In Process")

  ## Remove the "GO...." strings in the GO Process column
IX_supp_GO_Heatmap$GOs = IX_supp_GO_Heatmap$GOs %>% gsub(x = IX_supp_GO_Heatmap$GOs, pattern = " \\(.*\\)$", replacement = "")

  ## Make the digital (0,1) values factors so it's compatible to graph
IX_supp_GO_Heatmap$`Presence In Process` = as.factor(IX_supp_GO_Heatmap$`Presence In Process`)

  ## Create and store the heatmap core as a variable
  ## Use geom_tile to map the binary values
  ## Customize by removing the filler surrounding the graph and the tickmarks
  ## Center the Title

IX_supp_HEATlabs = ggplot(IX_supp_GO_Heatmap, aes(`GOs`, `DEGs`)) + 
  geom_tile(aes(fill = `Presence In Process`), color = "black") + 
  scale_fill_manual(values = c("navy", "gold3")) + 
  labs(x = "2018 GO Biological Process", y = "Differentially Expressed Genes (FDR < 0.05)", title = "GO Analysis of DEGs Affected by Tmem184b in vitro") + 
  guides(fill = guide_legend(title = "Presence in Process")) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0)) + 
  coord_flip()

  ## Plot the heatmap
IX_supp_HEATlabs



##### PANTHER IX Bar plots (opposite) #####
  ## Import the intersection files generated from PANTHERDB.ORG
Panther_supp_GO_Bio_Comp_Bar = read_csv("M:/PAPER ASSEMBLY/Itch paper/Molecular Mechanism/eDRGs/Panther GO Biological Processes Complete Overrepresented from genes indirectly affected by Tmem Trimmed for graph.csv")

  ## Re-order the terms by FDR (Q value aka Adjusted P-value)
Panther_supp_GO_Bio_Comp_Bar = Panther_supp_GO_Bio_Comp_Bar[order(Panther_supp_GO_Bio_Comp_Bar$FDR, decreasing = FALSE),]
  ## Remove unnecessary rows, columns
Panther_supp_GO_Bio_Comp_Bar = Panther_supp_GO_Bio_Comp_Bar[-15,-c(2,3,4,6)]
rownames(Panther_supp_GO_Bio_Comp_Bar) = NULL

  ## Transform the Q
Panther_supp_GO_Bio_Comp_Bar$neglog10 = -log10(Panther_supp_GO_Bio_Comp_Bar$FDR)

  ## Remove the ")" strings in the GO Process column
Panther_supp_GO_Bio_Comp_Bar$`GO biological process complete` = Panther_supp_GO_Bio_Comp_Bar$`GO biological process complete` %>% gsub(x = Panther_supp_GO_Bio_Comp_Bar$`GO biological process complete`, pattern = "^.*\\)", replacement = "")

  ## Take the top 10 of the data
Panther_supp_GO_Bio_Comp_Bar = Panther_supp_GO_Bio_Comp_Bar[1:10,]

  ## Condense GO strings
Panther_supp_GO_Bio_Comp_Bar$`GO biological process complete`[9] = "adenylate cyclase-inhibiting\nG protein-coupled\nreceptor signaling pathway"

  ## Create the bar plot
Panther_supp_GO_Bio_Comp_BAR = ggplot(Panther_supp_GO_Bio_Comp_Bar, aes(x = reorder(`GO biological process complete`, `neglog10`), y = `neglog10`)) +
  geom_col(stat = "identity" , aes(fill = `neglog10`)) +
  scale_fill_continuous(name = "-log10(Adj. P-val)") +
  labs(x = "Panther 2016 GO Biological Process Complete", y = "-log10(Adj. P-val)", title = "Pathways Containing In vitro DEGs Both Up in Mutant and Down upon Rescue") +
  theme_classic() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5), legend.position = "none") +
  coord_flip()

  ## Plot it
Panther_supp_GO_Bio_Comp_BAR


