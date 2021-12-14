WGCNA
================
Erik Larsen
11/22/2021

The following code was adapted from a compilation of sources, primarily
consisting of the tutorial from the
[WGCNA](https://cran.r-project.org/package=WGCNA) package developed by
Steve Horvath, Peter Langfelder, and Bin Zhang. This was originally a
walkthrough workflow for performing Weighted Gene Co-Expression Analysis
(gene network analysis). There is also a standalone R script and
downstream bioinformatics analysis scripts (GO and pathway analysis,
network analysis) based on the gene modules derived below.

## Environment Prep

Note that code can be sectioned and condensed with the `Alt+O` command.

List of packages for this script:
[BiocGenerics](https://bioconductor.org/packages/BiocGenerics),
[tidyverse](https://cran.r-project.org/package=tidyverse),
[plyr](https://cran.r-project.org/package=plyr),
[dplyr](https://cran.r-project.org/package=dplyr),
[reshape2](https://cran.r-project.org/package=reshape2),
[stringr](https://cran.r-project.org/package=stringr),
[readr](https://cran.r-project.org/package=readr),
[stats4](https://cran.r-project.org/package=stats4),
[matrixStats](https://cran.r-project.org/package=matrixStats),
[dendextend](https://cran.r-project.org/package=dendextend),
[ggdendro](https://cran.r-project.org/package=ggdendro),
[pheatmap](https://cran.r-project.org/package=pheatmap),
[fastcluster](https://cran.r-project.org/package=fastcluster),
[dynamicTreeCut](https://cran.r-project.org/package=dynamicTreeCut),
[grid](https://cran.r-project.org/package=grid),
[igraph](https://cran.r-project.org/package=igraph),
[WGCNA](https://cran.r-project.org/package=WGCNA),
[clValid](https://cran.r-project.org/package=clValid),
[biomaRt](https://bioconductor.org/packages/biomaRt/)

Install biology-based packages with `BiocManager`. Load the
[BiocGenerics](https://bioconductor.org/packages/BiocGenerics) package
for Bioconductor-relevant functionality (installing packages from
[Bioconductor](https://www.bioconductor.org/))

``` r
library(BiocGenerics) ## Needed to install and/or load Bioconductor packages
```

Load some of the packages.

``` r
## For data wrangling (slicing/adding/removing/melting/rearranging dataframes and their columns and rows):
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)

library(stringr) ## Awesome for manipulating strings

library(readr) ## For importing data and files

library(stats4)
library(matrixStats) ## Necessary for matrix multiplication, manipulation
```

    ## Warning: package 'matrixStats' was built under R version 4.1.1

``` r
library(dendextend) ## Additional flexibility in creating dendrograms (dashed line extensions and additional stuff)
library(ggdendro) ## Enables ggplot-quality dendrograms to be attached to heatmaps
library(pheatmap) ## Awesome for creating custom heatmaps
library(fastcluster)
library(dynamicTreeCut) ## Enables flexibly customizing and automatically cutting dendrograms

library(grid) ## grid has been removed from Cran and is incorporated into base R
library(igraph) ## igraph enables network visualization and customization
```

    ## Warning: package 'igraph' was built under R version 4.1.1

``` r
library(WGCNA) ## Creates topological overlap matrices and enables visualizing hierarchically clustered dendrograms and networks
library(clValid) ## Used for validating clustering (Dunn Index, Silhouette)

library(httpuv) ## For including plots in Markdown file output
```

## Data Prep

Import the differential expression dataset of interest and clean the
data.

``` r
## Tmem184b-GT/GT aDRG DGE
DESeq2_Adults = read.csv("M:/Erik/Data/Omics/RNAseq/Adult DRG/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results.csv")

## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
DESeq2_Adults3 = subset(DESeq2_Adults, (!is.na(DESeq2_Adults[,"AdjP"])))

## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
DESeq2_Adults9 = DESeq2_Adults3 %>% filter(!grepl(DESeq2_Adults3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
  ## First create the column and use Gene IDs as place-holders
DESeq2_Adults9$labs = DESeq2_Adults9$GeneID

  ## Replace DEGs with the string, "DEGs"
DESeq2_Adults9$labs[
  which(DESeq2_Adults9$AdjP <= 0.05)]= "DEGs"
    
  ## Replace the remaining genes with "Non-DEGs"
DESeq2_Adults9$labs[
  which(DESeq2_Adults9$AdjP >= 0.05)]= "Non-DEGs"

  ## Create a list of the gene names of differentially expressed genes
DEG_list = c(DESeq2_Adults9$GeneID[
  which(DESeq2_Adults9$labs == "DEGs")])
```

Import the TPM file (normalized counts file).

``` r
  ## Import the TPM file. These are expression estimates for each gene, for each sample/replicate, where each gene's value is normalized to its sample's effect size
Adult_Normalized_Counts = read_csv("M:/Erik/Data/Omics/RNAseq/Adult DRG/Processed Galaxy Output/Counts Files to Upload/RNASeqRepResults.csv", col_names = TRUE)
  ## Rename columns; should know this ahead of time
colnames(Adult_Normalized_Counts) = c("GeneID", "WT1", "WT2", "WT3", "WT4", "Mut1", "Mut2", "Mut3", "Mut4")

  ## Subset the genes most affected by Tmem mutation **if using a laptop** (laptops don't have enough RAM to process >5000)
DESeq2_Short = DESeq2_Adults9[c(1:5000),]
  
  ## Subset the TPM file by the 5000 most-affected genes
Short_Profile = Adult_Normalized_Counts[Adult_Normalized_Counts$GeneID %in% c(DESeq2_Short$GeneID),]
```

## Generate heatmaps from hierarchically clustered data

Create a function that will find row Z scores across the dataset.

Run the function to determine Row Z-score (gene-wise Z)

``` r
Find_Row_Z(Expression_Profile = Short_Profile)
```

Perform hierarchical clustering using Euclidean distance on the Z values
obtained above. Find Tmem184b’s row.

``` r
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
```

    ## [1] 2904

Use that index (Tmem’s row) to center a heatmap around Tmem184b. This
groups the genes most affected by the Tmem gene trap together.

``` r
  ## Create a heatmap around Tmem
pheatmap(mat = Corr_Z[2884:2924,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, labels_row = Euclid_dist_ord_Z_Genes[2884:2924])
```

![](WGCNA3_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

For a comprehensive transcriptional heatmap profile, use all genes

``` r
  ## Create a heatmap of all 5000 genes
pheatmap(mat = Corr_Z[,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, show_rownames = F)
```

![](WGCNA3_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

For non-centered data. This preserves natural expression variation
across replicates. Therefore grouping genes by expression quantification
isn’t a direct effect of the mutation and may just be random noise.
Maybe not as informative.

``` r
  ## Perform hierarchical clustering on the un-centered transcriptional profiles; use Euclidean distance
Euclid_dist_order = hclust(dist(Short_Profile[,2:9], method = "euclidean"))$order

  ## Find gene names
Euclid_dist_ord_Genes = c(Short_Profile$GeneID[Euclid_dist_order])

  ## Re-arrange the clustered profiles into a format compatible with heatmap construction
Corr = Short_Profile %>%
  mutate(GeneID = factor(GeneID, levels = Euclid_dist_ord_Genes)) %>%
  arrange(GeneID)

  ## Find where Tmem is
which(Euclid_dist_ord_Genes == "Tmem184b") # 3213
```

    ## [1] 3213

``` r
  ## Create a heatmap around Tmem
pheatmap(mat = Corr[3193:3233,2:9], color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", angle_col = 0, treeheight_row = 35, treeheight_col = 7, labels_row = Euclid_dist_ord_Genes[3193:3233])
```

![](WGCNA3_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Gene Co-expression Network Prep

First, create a function that will remove NaNs that affect the
construction of the network and prevents mathematical operations.

-   This depends on how many conditions and replicates there are in the
    dataset. In this example, 2 conditions (WT/Mut), 4 replicates each.

This function will return a cleaned dataframe to the global environment,
titled `Expression_Profile_cleaned`. It excludes genes that generate
NaNs because the genes have 0s across replicates

Remove NaN-generating genes by running that function.

``` r
Remove_NAN_Generators(Expression_Profile = Short_Profile, Samples_to_Compare = "All")
```

Use the “cleaned” profile to find what power to raise a correlation
(also called an adjacency) matrix.

-   By increasing the exponent power, it can increase the “connectivity”
    of more related genes from those less related.

Use the WGCNA function, `pickSoftThreshold`, to compute model fits of a
range of exponent powers and compute the “connectivity” metrics of the
data based on these power values. The following code will produce plots.

``` r
powers = c(seq(4,10,by=1), seq(12,20, by=2))

powerTable = list(data = pickSoftThreshold(data = t(Expression_Profile_cleaned[,2:9]), powerVector = powers, verbose = 2))
```

    ## pickSoftThreshold: will use block size 5000.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 5000 of 5000
    ##    Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      4    0.460 -1.09          0.942   460.0     419.0   1020
    ## 2      5    0.543 -1.20          0.953   326.0     291.0    788
    ## 3      6    0.611 -1.30          0.960   240.0     211.0    626
    ## 4      7    0.676 -1.40          0.970   183.0     158.0    509
    ## 5      8    0.717 -1.48          0.974   142.0     121.0    423
    ## 6      9    0.746 -1.56          0.975   113.0      95.3    357
    ## 7     10    0.783 -1.62          0.983    91.9      76.2    306
    ## 8     12    0.817 -1.76          0.985    63.0      51.1    233
    ## 9     14    0.849 -1.86          0.984    45.1      35.8    183
    ## 10    16    0.867 -1.94          0.984    33.5      26.1    148
    ## 11    18    0.879 -2.00          0.981    25.6      19.7    123
    ## 12    20    0.889 -2.04          0.979    20.1      15.3    103

``` r
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
```

    ## [1] 16

``` r
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
```

The “scale-free model fit” graphed as a function of soft threshold (aka
power) should have a y-axis ranging from 0 - 1. Anything larger than
this indicates a noisy dataset, suggesting an overfit or noisy model.

In general, this graph should also resemble a log curve. A power
selection should be made before a horizontal asymptote.

-   This should also be reflected in the maximum connectivity graphed as
    a function of power (power selection that retains the most
    connectivity before it “bottoms out” asymptotically).

The command `powerTable$data$powerEstimate` provides a suggested power
selection for the correlation/adjacency matrix. Usually the power is
4-7, anyway.

## Create the Adjacency and Topological Overlap Matrices (TOM)

The adjacency is a square matrix where the Pearson correlation of every
gene pair is calculated. Raise it to the power as determined above (6)
to amplify genes that correlate better as opposed to genes that don’t
correlate well.

``` r
  ## Create a weighted adjacency (no absolute value of the correlations raised to a power)
  ## Remember to transpose the expression matrix
ADJ = cor(t(Expression_Profile_cleaned[,2:9]), method = "pearson")^6
rownames(ADJ) = Expression_Profile_cleaned$GeneID ## Append gene names for the rows
colnames(ADJ) = Expression_Profile_cleaned$GeneID ## Append gene names for the columns
```

Create the signed TOM. If power were odd-numbered, then the TOM would
contain both positive and negative values.

-   Think of a TOM as an adjacency with context: the correlation of gene
    pairs in the context of their neighbors, extrapolated to all genes.
    -   Example: if looking at a 4-gene module, genes **A**, **B**,
        **C**, **D**, where **Gene A** and **Gene B** directly relate
        (high correlation), **Gene C** and **Gene D** directly relate,
        and **Gene B** and **Gene C** are related, the matrix scores the
        modules’ similarity and weights the direct links strongest

``` r
  ## Create the topological overlap matrix using the WGCNA function "TOMsimilarity": converts an adjacency into a TOM
TOM = TOMsimilarity(adjMat = ADJ, TOMType = "signed")
```

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
rownames(TOM) = rownames(ADJ)
colnames(TOM) = colnames(ADJ)
```

Create the signed dissimilarity TOM; some find this better segregates
modules

``` r
  ## Create the topological overlap dissimilarity matrix
disTOM = 1-TOMsimilarity(adjMat = ADJ, TOMType = "signed")
```

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
rownames(disTOM) = rownames(ADJ)
colnames(disTOM) = rownames(ADJ)
```

## Visualize the WGCNAs

Visualize the results to evaluate the constructed network(s)

``` r
  ## Create a hierarchically clustered dendrogram based on the DISsimilarity (of the TOM); distinguishes clusters better than similarity
Tree = hclust(as.dist(disTOM), method = "average")
  ## Determine clusters to plot based on a "dynamic dendrogram cut"
unmergedLabs = cutreeDynamic(dendro = Tree, distM = disTOM,
                             deepSplit = 1,
                             cutHeight = 0.995,
                             pamRespectsDendro = F)
```

    ##  ..done.

``` r
  ## Convert the clusters from index labels to colors
unmergedCols = labels2colors(unmergedLabs)
  ## Merge clusters that are similar enough based on the subjective/arbitrary "cutHeight" parameter provided within the WGCNA function, "mergeCloseModules"
merge = mergeCloseModules(t(Expression_Profile_cleaned[,2:9]), unmergedLabs, cutHeight = 0.10, verbose = 2) ## cutHeight is the parameter that determines closeness of modules; a smaller cutHeight is more selective- closer to 1 merges more unrelated modules
```

    ##  mergeCloseModules: Merging modules whose distance is less than 0.1
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...

``` r
  ## Convert color indeces into a string
moduleLabs = merge$colors
  ## Convert the string of labels into colors
moduleCols = labels2colors(moduleLabs)
  ## Find the new, merged eigenvectors (modules)
MEs = merge$newMEs

#t(Expression_Profile_cleaned[,2:9])

  ## Find Tmem index
which(Expression_Profile_cleaned$GeneID == "Tmem184b")
```

    ## [1] 4395

``` r
  ## Plot the dendrogram and associated clusters/modules
plotDendroAndColors(Tree, cbind(unmergedCols, moduleCols), 
                    c("Unmerged Modules", "Merged Modules"), 
                    dendroLabels = F, 
                    hang = 0.001,
                    addGuide = T,
                    guideHang = 0.001,
                    main = "Tmem GT DRG Dendrogram")
```

![](WGCNA3_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Investigate the Tmem modules:

``` r
#which(Expression_Profile_cleaned$GeneID == "Il31ra") ## For finding the modules of genes of interest
  ## Determine Tmem's modules
moduleCols[4395]
```

    ## [1] "brown"

``` r
unmergedCols[4395]
```

    ## [1] "yellow"

``` r
  ## Find the indeces from the clustered data of genes within Tmem's clusters
TmemUnmergedModule = c(which(unmergedCols == "yellow"))
TmemModule = c(which(moduleCols == "brown"))
  ## Find the genes of those indeces
TmemUnmergedModule = c(Expression_Profile_cleaned$GeneID[TmemUnmergedModule])
TmemModule = c(Expression_Profile_cleaned$GeneID[TmemModule])
```

``` r
  ## The adjacency values of the genes in the unmerged Tmem module. High values mean higher similarity
ADJ[4395,c(which(unmergedCols == "yellow"))]
```

    ##    A3galt2      Adprh       Alg5     Anapc5      Apaf1    Arhgdia     Arid3a 
    ## 0.83272713 0.30514206 0.31840805 0.31639124 0.59732857 0.90796586 0.32865611 
    ##       Art3      Cadm1      Camk4      Capn1    Carhsp1     Carns1     Casp12 
    ## 0.43434645 0.56625628 0.79885239 0.70121926 0.93525839 0.51457464 0.35326659 
    ##      Cd24a       Cd44       Cd55      Cdh15       Clgn     Cobll1      Cpne3 
    ## 0.84319630 0.60751007 0.51522464 0.42345850 0.78248260 0.56546871 0.55306077 
    ##    Ctdspl2       Ctsa      Ctxn3       Dgka       Dgkh       Dgkz       Dhdh 
    ## 0.80513966 0.69152085 0.64000263 0.74506053 0.59800821 0.55055616 0.22371840 
    ##     Dusp26      Efcc1      Efna1      Esyt1    Fam167a     Fam83d      Fbln2 
    ## 0.68015705 0.65977657 0.51228422 0.64870108 0.78145452 0.75440883 0.64051555 
    ##       Fez2       Fmr1      Gna14       Gng2   Hoxd3os1       Hras       Hrh2 
    ## 0.10335941 0.50217801 0.66469423 0.77713683 0.36547858 0.77418284 0.44235195 
    ##      Htr1a      Htr1f     Il31ra       Isl2      Itpr3      Kcng3     Kcnip4 
    ## 0.94962634 0.64322434 0.84605783 0.69133146 0.44401301 0.58372674 0.52860515 
    ##       Klf5       Ldb2      Lima1      Lpar3       Ly86       Mal2     Mcm3ap 
    ## 0.83239233 0.65028898 0.66459698 0.82106374 0.69674954 0.65454495 0.33846685 
    ##      Mgat3      Moxd1      Myh10    Ngfrap1       Nkap      Nr1d2       Nsg2 
    ## 0.47897222 0.65893652 0.59276737 0.66086299 0.08635315 0.64008205 0.52635370 
    ##     Osbpl3     Osbpl8      Paqr6      Pcbp4      Plcb3     Plcxd3      Plpp6 
    ## 0.56131169 0.73636031 0.54040069 0.74725137 0.81867550 0.32849889 0.83108645 
    ##     Plxnc1      Prdm8        Ret       Rhov       Rmrp      Rn7sk       Scg2 
    ## 0.58936426 0.54743368 0.74187132 0.91523389 0.79152543 0.75418012 0.59056699 
    ##       Scg3     Scn10a     Scn11a    Serinc2      Sf3b2    Slc35f3    Slc37a1 
    ## 0.72463544 0.69709629 0.66026639 0.67780969 0.72714263 0.52747490 0.67853483 
    ##     Stard9      Stat3      Sulf2      Synpr     Tapbpl      Tgfb3   Tmem151b 
    ## 0.48914807 0.71091233 0.68789022 0.80511872 0.44041144 0.35498270 0.77907001 
    ##    Tmem158   Tmem184b    Tmem233        Tnr   Traf3ip2     Trim36      Trnp1 
    ## 0.65724379 1.00000000 0.88488841 0.74490313 0.68986636 0.87000015 0.77640772 
    ##      Trpc6      Tshz2     Tubb2b       Txn1      Vps18     Zbtb41    Zcchc18 
    ## 0.58998866 0.77177574 0.84081448 0.58925922 0.72043357 0.14781876 0.61673222 
    ##    Zdhhc18     Zfp827 
    ## 0.60505017 0.69094927

``` r
  ## The TOM values of the genes in the unmerged Tmem module. High values mean higher similarity
TOM[4395, c(which(unmergedCols == "yellow"))]
```

    ##   A3galt2     Adprh      Alg5    Anapc5     Apaf1   Arhgdia    Arid3a      Art3 
    ## 0.3144887 0.2248785 0.2399185 0.2433551 0.2812256 0.3185770 0.2340291 0.2543616 
    ##     Cadm1     Camk4     Capn1   Carhsp1    Carns1    Casp12     Cd24a      Cd44 
    ## 0.2768624 0.3042338 0.2975909 0.3478533 0.2718560 0.2536838 0.3413768 0.2896568 
    ##      Cd55     Cdh15      Clgn    Cobll1     Cpne3   Ctdspl2      Ctsa     Ctxn3 
    ## 0.2723622 0.2664625 0.3141462 0.2840016 0.2740114 0.3131297 0.2996596 0.2914593 
    ##      Dgka      Dgkh      Dgkz      Dhdh    Dusp26     Efcc1     Efna1     Esyt1 
    ## 0.2941839 0.3050385 0.2761660 0.2116566 0.2965132 0.2881033 0.2744715 0.2863115 
    ##   Fam167a    Fam83d     Fbln2      Fez2      Fmr1     Gna14      Gng2  Hoxd3os1 
    ## 0.3071810 0.3019307 0.2853216 0.1857873 0.2677486 0.3094008 0.3085760 0.2435092 
    ##      Hras      Hrh2     Htr1a     Htr1f    Il31ra      Isl2     Itpr3     Kcng3 
    ## 0.3075845 0.2716093 0.3301151 0.2898748 0.3126815 0.3036309 0.2623828 0.2997978 
    ##    Kcnip4      Klf5      Ldb2     Lima1     Lpar3      Ly86      Mal2    Mcm3ap 
    ## 0.2718596 0.3207743 0.3103710 0.2994885 0.3384732 0.2962441 0.2921793 0.2337064 
    ##     Mgat3     Moxd1     Myh10   Ngfrap1      Nkap     Nr1d2      Nsg2    Osbpl3 
    ## 0.2767245 0.2875109 0.2789044 0.2965377 0.1660100 0.2817616 0.2715264 0.2853397 
    ##    Osbpl8     Paqr6     Pcbp4     Plcb3    Plcxd3     Plpp6    Plxnc1     Prdm8 
    ## 0.3005923 0.2799283 0.2954304 0.3206961 0.2483783 0.3089705 0.2900199 0.2720399 
    ##       Ret      Rhov      Rmrp     Rn7sk      Scg2      Scg3    Scn10a    Scn11a 
    ## 0.3183848 0.3457090 0.3111150 0.3035245 0.2759661 0.3032903 0.3163879 0.3187633 
    ##   Serinc2     Sf3b2   Slc35f3   Slc37a1    Stard9     Stat3     Sulf2     Synpr 
    ## 0.2915631 0.2922522 0.2766179 0.2898838 0.2718132 0.3021945 0.2863415 0.3245482 
    ##    Tapbpl     Tgfb3  Tmem151b   Tmem158  Tmem184b   Tmem233       Tnr  Traf3ip2 
    ## 0.2609314 0.2409881 0.3234383 0.3008143 1.0000000 0.3183205 0.3020523 0.2969370 
    ##    Trim36     Trnp1     Trpc6     Tshz2    Tubb2b      Txn1     Vps18    Zbtb41 
    ## 0.3193978 0.3046786 0.2804652 0.3106010 0.3143706 0.2861631 0.2940486 0.1900470 
    ##   Zcchc18   Zdhhc18    Zfp827 
    ## 0.2823532 0.2864232 0.3182126

``` r
  ## The disTOM values of the genes in the unmerged Tmem module. Low values mean higher similarity
disTOM[4395, c(which(unmergedCols == "yellow"))]
```

    ##   A3galt2     Adprh      Alg5    Anapc5     Apaf1   Arhgdia    Arid3a      Art3 
    ## 0.6855113 0.7751215 0.7600815 0.7566449 0.7187744 0.6814230 0.7659709 0.7456384 
    ##     Cadm1     Camk4     Capn1   Carhsp1    Carns1    Casp12     Cd24a      Cd44 
    ## 0.7231376 0.6957662 0.7024091 0.6521467 0.7281440 0.7463162 0.6586232 0.7103432 
    ##      Cd55     Cdh15      Clgn    Cobll1     Cpne3   Ctdspl2      Ctsa     Ctxn3 
    ## 0.7276378 0.7335375 0.6858538 0.7159984 0.7259886 0.6868703 0.7003404 0.7085407 
    ##      Dgka      Dgkh      Dgkz      Dhdh    Dusp26     Efcc1     Efna1     Esyt1 
    ## 0.7058161 0.6949615 0.7238340 0.7883434 0.7034868 0.7118967 0.7255285 0.7136885 
    ##   Fam167a    Fam83d     Fbln2      Fez2      Fmr1     Gna14      Gng2  Hoxd3os1 
    ## 0.6928190 0.6980693 0.7146784 0.8142127 0.7322514 0.6905992 0.6914240 0.7564908 
    ##      Hras      Hrh2     Htr1a     Htr1f    Il31ra      Isl2     Itpr3     Kcng3 
    ## 0.6924155 0.7283907 0.6698849 0.7101252 0.6873185 0.6963691 0.7376172 0.7002022 
    ##    Kcnip4      Klf5      Ldb2     Lima1     Lpar3      Ly86      Mal2    Mcm3ap 
    ## 0.7281404 0.6792257 0.6896290 0.7005115 0.6615268 0.7037559 0.7078207 0.7662936 
    ##     Mgat3     Moxd1     Myh10   Ngfrap1      Nkap     Nr1d2      Nsg2    Osbpl3 
    ## 0.7232755 0.7124891 0.7210956 0.7034623 0.8339900 0.7182384 0.7284736 0.7146603 
    ##    Osbpl8     Paqr6     Pcbp4     Plcb3    Plcxd3     Plpp6    Plxnc1     Prdm8 
    ## 0.6994077 0.7200717 0.7045696 0.6793039 0.7516217 0.6910295 0.7099801 0.7279601 
    ##       Ret      Rhov      Rmrp     Rn7sk      Scg2      Scg3    Scn10a    Scn11a 
    ## 0.6816152 0.6542910 0.6888850 0.6964755 0.7240339 0.6967097 0.6836121 0.6812367 
    ##   Serinc2     Sf3b2   Slc35f3   Slc37a1    Stard9     Stat3     Sulf2     Synpr 
    ## 0.7084369 0.7077478 0.7233821 0.7101162 0.7281868 0.6978055 0.7136585 0.6754518 
    ##    Tapbpl     Tgfb3  Tmem151b   Tmem158  Tmem184b   Tmem233       Tnr  Traf3ip2 
    ## 0.7390686 0.7590119 0.6765617 0.6991857 0.0000000 0.6816795 0.6979477 0.7030630 
    ##    Trim36     Trnp1     Trpc6     Tshz2    Tubb2b      Txn1     Vps18    Zbtb41 
    ## 0.6806022 0.6953214 0.7195348 0.6893990 0.6856294 0.7138369 0.7059514 0.8099530 
    ##   Zcchc18   Zdhhc18    Zfp827 
    ## 0.7176468 0.7135768 0.6817874

Plot the TOM:

``` r
  ## Create a dendrogram compatible with TOM-TOM plot creation
colorDynamicTOM = labels2colors(cutreeDynamic(Tree, method ="tree"))
  ## Make the diagonal of the TOM NA to re-scale the plot
diag(disTOM) = 0
  ## Create the TOM-TOM plot
TOMplot((disTOM^3), Tree, as.character(moduleCols), main = "Dissim. TOM Plot Top5000 Genes")
```

![](WGCNA3_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Analyze the genes in the merged module using enrichr with the
Enrichr\_Analysis script and function
