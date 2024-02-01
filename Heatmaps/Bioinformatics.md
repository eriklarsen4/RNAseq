Bioinformatics
================
Erik Larsen
7/17/2021

# Overview

-   I developed this Markdown to document some bioinformatics analyses
    of data published by **M. Bhattacharya Lab** in the publication,
    [Transmembrane protein TMEM184B is necessary for
    interleukin-31–induced
    itch](https://journals.lww.com/pain/Abstract/9000/Transmembrane_protein_TMEM184B_is_necessary_for.97918.aspx).

-   These analyses were performed on a number of different datasets and
    include:

    -   MA plots
    -   volcano plots
    -   heatmaps
    -   bar graphs

1.  `Illumina FASTQ files` were pseudoaligned with the `Salmon`
    algorithm

2.  run through `Usegalaxy.org`’s `DESeq2` wrapper algorithm to obtain
    `DGEA files`

3.  `Galaxy` also returned TPM data for use in transcriptional profiling
    (heatmaps)

-   To obtain pathway and gene ontology info, the
    [Enrichr](https://maayanlab.cloud/Enrichr/) database was queried
    from the console using the
    [enrichR](https://github.com/guokai8/Enrichr#:~:text=Description%20EnrichR%20is%20a%20package%20can%20be%20used,species%20pubished%20by%20ENSEMBL%20and%20included%20with%20Bioconductor.)
    package.

-   The source script of this file is the [Bioinformatics R
    Script](C:/Users/Erik/Desktop/BoxCopy/Programming%20Scripts%20and%20Data/Bio/Scripts/R/Broad/Bioinformatics%20Script.R).

-   Standalone analyses or tutorials, including volcano plotting, and
    `gene ontology` and `pathway analysis` can be found in other
    `R/Github Markdown` files and scripts, [Gene Ontology Analysis
    Snippet](https://github.com/eriklarsen4/RNAseq/blob/master/GO%20Analysis/GO-Analysis.md)
    and [Volcano plot tutorial R
    Markdown](https://github.com/eriklarsen4/ggplot-scripts/blob/master/Bioinformatics/Volcano-plot-tutorial.md)

# Environment Prep

Note that code can be sectioned and condensed with the `Alt + O`
command.

Packages for this script:
[tidyverse](https://cran.r-project.org/package=tidyverse),
[stringr](https://cran.r-project.org/package=stringr),
[readr](https://cran.r-project.org/package=readr),
[biomaRt](https://bioconductor.org/packages/biomaRt/),
[GO.db](http://bioconductor.riken.jp/packages/3.0/data/annotation/html/GO.db.html),
[PANTHER.db](https://bioconductor.org/packages/PANTHER.db/),
[BiocGenerics](https://bioconductor.org/packages/BiocGenerics),
[pheatmap](https://cran.r-project.org/package=pheatmap)

## Install and Load Packages

Install biology-based packages with `BiocManager`.

Load the [BiocGenerics](https://bioconductor.org/packages/BiocGenerics)
package for Bioconductor-relevant functionality (installing packages
from [Bioconductor](https://www.bioconductor.org/))

``` r
library(BiocGenerics) ## Needed to install and/or load Bioconductor packages
```

Load the remaining packages

``` r
library(enrichR) ## Taps the Enrichr database; much better utility than the website

  ## For data wrangling (slicing/adding/removing/melting/rearranging dataframes and their columns and rows):
library(tidyverse)
# library(plyr)
# library(dplyr)
# library(reshape2)

library(readr) ## For importing data and files

library(stringr) ## Awesome for manipulating strings

  ## Databases for downstream analyses of lists via gene ontology:
library(biomaRt)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GO.db)

library(pheatmap) ## For creating awesome heatmaps

library(httpuv) ## For including plots in Markdown file output
```

## Upload Enrichr DBs

Obtain databases from `enrichr` by creating a variable and viewing it.
Subset from the variable a few databases of interest.

``` r
    ## Store it as a variable for better visualization and eventual subsetting
DBs = listEnrichrDbs()

  ## View the variable in the editor to find the relevant databases and their indeces for subsetting
  ## (click on the dataframe in the Global Environment and manually peruse)

diff_DBs = DBs %>%
  dplyr::filter(grepl(libraryName, pattern = "Transcription_Factor_PPIs|Reactome_2022|Most_Popular_Genes|ENCODE_TF|GPCR|Enrichr")) %>%
  dplyr::select(libraryName) %>%
  as.vector() %>%
  unlist() %>%
  as.character()

DBs = DBs %>%
  dplyr::filter(grepl(libraryName, pattern = "GO_.+?(2023)$|Panther_2016|PPI_Hub|Jensen_DISEASES")) %>%
  dplyr::select(libraryName) %>%
  as.vector() %>%
  unlist() %>%
  as.character()
```

Create dataframes to house bioinformatics analyses returned from
enrichr.

``` r
  ## Create dataframes to store returned analyses as variables
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
```

## Source Enrichr Function

Create a data frame to collect data returned by `enrichr`; then create
the function that will tap enrichr’s database and populate all the
dataframes.

``` r
source("~/GitHub/Proteomics/Functions/Enrichr Analysis Function.R")
```

## Data import

Import the adult `DRG DESeq2 file`

``` r
aDRG = read.csv("~/GitHub/ggplot-scripts/Bioinformatics/RNAseq Data Files/DESeq2 Expression Results.csv")

  ## Filter (subset) genes that went undetected or were outliers in terms of counts;
  ## new dataframe should not contain any NAs in p-value columns
aDRG3 = subset(aDRG, (!is.na(aDRG[,"AdjP"])))
 ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs, along with pseudogenes, etc.
aDRG9 = aDRG3 %>% filter(!grepl(GeneID,
                                pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))
#mt.+.?$|  <-- string identifier for mitochondrial tRNAs

  ## Create a column in the DESeq2 dataframe that scales the Adjusted P-value by log-base 10
aDRG9$log10ADJP = -log10(aDRG9$AdjP)

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
aDRG9$g.o.i. = aDRG9$GeneID
  ## Replace DEGs with the string, "DEGs"
aDRG9$g.o.i.[which(aDRG9$AdjP <= 0.05)]= "DEGs"
  ## Replace the remaining genes with "Non-DEGs"
aDRG9$g.o.i.[which(aDRG9$AdjP >= 0.05)]= "Non-DEGs"

aDRG_DEG_list = c(aDRG9$GeneID[which(aDRG9$g.o.i. == "DEGs")])
```

# Query the Enrichr db

Search the `Enrichr` database using the function created above with any
of the above gene lists.

``` r
Enrichr_Analysis(Genes = aDRG_DEG_list)
```

    ## Uploading data to Enrichr... Done.
    ##   Querying GO_Biological_Process_2023... Done.
    ## Parsing results... Done.
    ## Uploading data to Enrichr... Done.
    ##   Querying GO_Cellular_Component_2023... Done.
    ## Parsing results... Done.
    ## Uploading data to Enrichr... Done.
    ##   Querying GO_Molecular_Function_2023... Done.
    ## Parsing results... Done.
    ## Uploading data to Enrichr... Done.
    ##   Querying Panther_2016... Done.
    ## Parsing results... Done.
    ## Uploading data to Enrichr... Done.
    ##   Querying PPI_Hub_Proteins... Done.
    ## Parsing results... Done.
    ## Uploading data to Enrichr... Done.
    ##   Querying Jensen_DISEASES... Done.
    ## Parsing results... Done.
    ## Uploading data to Enrichr... Done.
    ##   Querying ENCODE_TF_ChIP-seq_2015... Done.
    ## Parsing results... Done.
    ## Uploading data to Enrichr... Done.
    ##   Querying Enrichr_Submissions_TF-Gene_Coocurrence... Done.
    ## Parsing results... Done.
    ## Uploading data to Enrichr... Done.
    ##   Querying L1000_Kinase_and_GPCR_Perturbations_down... Done.
    ## Parsing results... Done.
    ## Uploading data to Enrichr... Done.
    ##   Querying L1000_Kinase_and_GPCR_Perturbations_up... Done.
    ## Parsing results... Done.

    ##  [1] "GO BIOLOGICAL PROCESSES"                                                          
    ##  [2] "Lipid Phosphorylation"                                                            
    ##  [3] "Negative Regulation Of Regulated Secretory Pathway"                               
    ##  [4] "Calcium Ion-Regulated Exocytosis Of Neurotransmitter"                             
    ##  [5] "Adenylate Cyclase-Activating Dopamine Receptor Signaling Pathway"                 
    ##  [6] "Vocal Learning"                                                                   
    ##  [7] "Imitative Learning"                                                               
    ##  [8] "Negative Regulation Of Monocyte Chemotaxis"                                       
    ##  [9] "Peripheral Nervous System Neuron Development"                                     
    ## [10] "Regulation Of Regulated Secretory Pathway"                                        
    ## [11] "Adenylate Cyclase-Inhibiting Serotonin Receptor Signaling Pathway"                
    ## [12] "GO CELL COMPONENTS"                                                               
    ## [13] "Heterotrimeric G-protein Complex"                                                 
    ## [14] "Extrinsic Component Of Cytoplasmic Side Of Plasma Membrane"                       
    ## [15] "Outer Dynein Arm"                                                                 
    ## [16] "Potassium Channel Complex"                                                        
    ## [17] "Exocytic Vesicle Membrane"                                                        
    ## [18] "Neuron Projection"                                                                
    ## [19] "Synaptic Vesicle Membrane"                                                        
    ## [20] "Axon"                                                                             
    ## [21] "Neuronal Dense Core Vesicle"                                                      
    ## [22] "Nuclear Envelope Lumen"                                                           
    ## [23] "GO MOLECULAR FUNCTIONS"                                                           
    ## [24] "Diacylglycerol Kinase Activity"                                                   
    ## [25] "Ciliary Neurotrophic Factor Receptor Binding"                                     
    ## [26] "Calcium-Dependent Phospholipid Binding"                                           
    ## [27] "Inositol 1,4,5 Trisphosphate Binding"                                             
    ## [28] "Lipid Kinase Activity"                                                            
    ## [29] "Protein Kinase C Activity"                                                        
    ## [30] "Potassium Ion Leak Channel Activity"                                              
    ## [31] "Calcium-Dependent Protein Serine/Threonine Kinase Activity"                       
    ## [32] "Store-Operated Calcium Channel Activity"                                          
    ## [33] "Leak Channel Activity"                                                            
    ## [34] "PANTHER PATHWAYS"                                                                 
    ## [35] "Histamine H1 receptor mediated signaling pathway"                                 
    ## [36] "Heterotrimeric G-protein signaling pathway-Gq alpha and Go alpha mediated pathway"
    ## [37] "Thyrotropin-releasing hormone receptor signaling pathway"                         
    ## [38] "Oxytocin receptor mediated signaling pathway"                                     
    ## [39] "Endothelin signaling pathway"                                                     
    ## [40] "Endogenous cannabinoid signaling"                                                 
    ## [41] "Synaptic vesicle trafficking"                                                     
    ## [42] "Muscarinic acetylcholine receptor 1 and 3 signaling pathway"                      
    ## [43] "Angiotensin II-stimulated signaling through G proteins and beta-arrestin"         
    ## [44] "5HT2 type receptor mediated signaling pathway"                                    
    ## [45] "PPIs"                                                                             
    ## [46] "ITGB1"                                                                            
    ## [47] "MAP3K7"                                                                           
    ## [48] "PRKCA"                                                                            
    ## [49] "DLG4"                                                                             
    ## [50] "PRKACA"                                                                           
    ## [51] "GRIN2B"                                                                           
    ## [52] "PRKG1"                                                                            
    ## [53] "LCK"                                                                              
    ## [54] "GRIN1"                                                                            
    ## [55] "PIK3CA"                                                                           
    ## [56] "JENSEN DISEASE PHENOTYPES"                                                        
    ## [57] "Adrenal gland cancer"                                                             
    ## [58] "Drug psychosis"                                                                   
    ## [59] "Endometriosis of ovary"                                                           
    ## [60] "Autosomal recessive non-syndromic intellectual disability"                        
    ## [61] "Frontal lobe epilepsy"                                                            
    ## [62] "Kawasaki disease"                                                                 
    ## [63] "Asthenopia"                                                                       
    ## [64] "Esophageal varix"                                                                 
    ## [65] "Thyroid cancer"                                                                   
    ## [66] "Multiple endocrine neoplasia type 1"

## Prep data for heatmaps

Import the adult transcriptional profile (normalized TPM).

``` r
  ## Import the adult normalized counts file.
    ## These are expression estimates for each gene, for each sample/replicate,
    ## where each gene's value is normalized to its sample's effect size
aTPM = read_csv("~/GitHub/ggplot-scripts/Bioinformatics/RNAseq Data Files/RNASeqRepResults.csv", col_names = TRUE)
  ## Rename columns
colnames(aTPM) = c("GeneID", "WT1", "WT2", "WT3", "WT4", "Mut1", "Mut2", "Mut3", "Mut4")

  ## Subset by only genes in the filtered aDRG DESeq2 file
ExpProfile = aTPM[aTPM$GeneID %in% c(aDRG9$GeneID),]
```

Create a function that makes a dataframe housing the Z-scored
transcriptional profile to then graph using `pheatmap`.

Run the Z-score function.

``` r
Find_Row_Z(Expression_Profile = ExpProfile)
```

Arrange the data for generating a full transcriptional profile.

``` r
  ## Don't remove non-DEGs
Z = Z %>%
  filter(GeneID %in% aDRG9$GeneID[])

  ## Create a list of gene names ordered by euclidean distance by which to re-order the dataframe
Euclid_dist_order = hclust(dist(Z[,c(2:9)], method = "euclidean"))$order
  ## The names (not the numbers)
Euclid_dist_ord_Genes = c(Z$GeneID[Euclid_dist_order])

  ## Transform again to order by clusters (Euclidean distances)
Za = Z %>%
  mutate(GeneID =  factor(GeneID, levels = Euclid_dist_ord_Genes)) %>%
  arrange(GeneID)

  ## Find where the "Itch-related DEGs" are in the clustered matrix subsetted to DEGs for the heatmap
itch_DEGs = c("Il31ra", "Cysltr2", "Npy2r", "Sst",
              "Htr1a", "P2rx3", "Lpar3", "Lpar5",
              "Scn11a", "Scn10a", "Mrgprd", "Trpc6",
              "Trpc3", "F2rl2", "Htr1f", "Osmr",
              "Fxyd2", "Htr4", "Mrgprx1", "Ptgdr",
              "Trpa1", "Trpm6", "Hrh1", "Mrgpra3", "Tmem184b")

temp = vector()
for( i in 1:length(itch_DEGs)){
  temp[i] = which(Euclid_dist_ord_Genes == itch_DEGs[i])
}
itch_index = temp
#itch_index

  ## Confirm by indexing; copy and paste to customize the clustered heatmap
#Euclid_dist_ord_Genes[c(itch_index)]
  ## Confirm the location/order of the transformed matrix/dataframe
  ## Should be "Il31ra"
Za[itch_index[1],]
```

    ## # A tibble: 1 × 9
    ##   GeneID   WT1   WT2   WT3   WT4   Mut1   Mut2   Mut3   Mut4
    ##   <fct>  <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
    ## 1 Il31ra  1.19  1.18 0.441 0.830 -0.943 -0.859 -0.927 -0.909

# Plot the heatmaps

## Adult Full

View the full transcriptional profile of **Adult Tmem184b-mutant DRG
neurons**.

``` r
 ## Visualize the Z-scored, Euclidean-clustered and ordered gene TPMs across replicates with "pheatmap"
pheatmap(mat = Za[,2:9],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         cutree_rows = 2,
         treeheight_row = 35,
         treeheight_col = 7,
         show_rownames = F)
```

![](Bioinformatics_files/figure-gfm/Full%20Z-scored,%20Euclidean-clustered%20and%20ordered%20pheatmap-1.png)<!-- -->

## Adult Zoom

**View the cluster near Tmem184b**

``` r
  ## Visualize around the Tmem184b cluster
pheatmap(mat = Za[1597:1637,2:9],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         treeheight_row = 35,
         treeheight_col = 9,
         labels_row = Euclid_dist_ord_Genes[1597:1637])
```

![](Bioinformatics_files/figure-gfm/Tmem184b%20cluster%20Zoom%20on%20Z-scored%20Euclidean-clustered%20and%20ordered%20pheatmap-1.png)<!-- -->

## Adult Itch

**View the cluster including select itch transcripts**.

``` r
  ## Visualize the "Itch-related DEGs"
pheatmap(mat = Za[c(itch_index),2:9],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         treeheight_row = 35,
         treeheight_col = 9,
         cutree_rows = 4,
         labels_row = Euclid_dist_ord_Genes[c(itch_index)])
```

![](Bioinformatics_files/figure-gfm/Itch%20trx%20cluster%20Zoom%20on%20Z-scored%20Euclidean-clustered%20and%20ordered%20pheatmap-1.png)<!-- -->

## Adult DEGs Only

Arrange the data to include only DEGs.

``` r
  ## Transform the profile to include only DEGs
Zb = Z %>%
  filter(GeneID %in% aDRG9$GeneID[c(1:376)])
```

Prep it for a heatmap as before (not shown).

**Plot the results of the profile of only DEGs in Tmem184b-mutant DRG
neurons.**

![](Bioinformatics_files/figure-gfm/Z-scored%20Euclidean-clustered%20and%20ordered%20pheatmap%20of%20DEGS-1.png)<!-- -->

## Adult DEGs Only Zoom

**Zoom to Tmem184b**

![](Bioinformatics_files/figure-gfm/Tmem184b%20Cluster%20Zoom%20Z-scored%20Euclidean-clustered%20ordered%20pheatmap%20of%20DEGs-1.png)<!-- -->

## Adult DEGs Only Itch Zoom

**Zoom to the area comprising many itch genes**

![](Bioinformatics_files/figure-gfm/Itch%20DEG%20Cluster%20Zoom%20Z-scored%20Euclidean-clustered%20ordered%20pheatmap%20of%20DEGs-1.png)<!-- -->

**View select itch genes**.

-   The difference between the first three heatmaps and the last four
    heatmaps is hierarchically clustering only on DEGs in the last four

-   Clustering was performed on all genes in the first three

![](Bioinformatics_files/figure-gfm/Itch%20DEG%20only%20of%20Z-scored%20Euclidean-clustered%20ordered%20pheatmap%20of%20DEGs-1.png)<!-- -->
