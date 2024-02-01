Volcano plots with ggplot
================
Erik Larsen
7/16/2021

I developed the following code as a tutorial for graphical analysis of long-read, `Illumina` transcriptomic data in `Dr. Martha Bhattacharya`'s lab at UA. We processed the reads on a 3rd party server,
[Usegalaxy.org](https://usegalaxy.org/), using a bioinformatics standard
differential expression algorithm, known as
[DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).

The plots, generated with the [ggplot2 R
package](https://cran.r-project.org/package=ggplot2) package are standard for differential gene expression analysis, but I've also added a bit more customization. This `Github R Markdown` is a generic tutorial for those less-experienced with `R` and/or plotting and deriving insights from RNA-seq data.

## Environment Prep

Upload the packages that include `ggplot2` ("Grammar of Graphics")

  + `ggplot2` contains a set of functions that produce images, by mapping data
to aesthetic properties (e.g. colors to gene symbols)

``` r
  ## Install the tidyverse and reshape2 packages if not already done
    ## For example
#install.packages("tidyverse")
#install.packages("reshape2")

  # Load them from your R package library
## For data wrangling
  ## (manipulating dataframes):
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)

library(ggplot2) ## For awesome, publication-quality graphs
library(ggrepel) ## For labeling mapped data on those plots
library(stats4)
library(readr) ## For importing data and files
```

Our paired-end RNA-seq transcript reads were first aligned using `Salmon`, a pseudo alignment algorithm, in a virtual computing environment

  + basically, `Salmon` assigns ("maps") all the read fragments (`fastq` files) from an experiment according to a provided reference sequence
  +   in this case, the reference sequence is the mouse transcriptome
  + having aligned all the read fragments, `Salmon` quantifies each gene's total number of fragments:
      + how many `Gene_X` read fragments aligned to the transcriptome reference `Gene_X` sequence
  + in this analysis, the reads do *not* account for alternative splicing, though `Salmon` is capable of doing so

Import the `Usegalaxy.org` data frame from the lab server that contains
mean transcript counts of all the *Tmem184b*<sup>GT/GT</sup> mice.

``` r
aDRG_TPM = read_csv("https://github.com/eriklarsen4/ggplot-scripts/blob/master/Bioinformatics/RNAseq%20Data%20Files/RNASeqRepResults.csv", col_names = c("GeneID", "WT", "Mut"))
```

Import the differential expression analysis data if not already done
(`Import Dataset` in the `Environment Window`); it is directly called
via the script below. These are the `DESeq2` results from
`Usegalaxy.org`, computing a `Wald Test` on `Salmon` pseudoaligned read counts
files to determine which genes are differentially expressed.

-   Before importing into `R`, **save the Galaxy DESeq2 Output**
    (dowloaded from `Galaxy.org` as a `.txt` file) **as a CSV** with headings to
    import into the `R Global Environment` to then manipulate for
    plotting purposes.

-   Also, in `Microsoft Excel`, add a column that converts `log2 FC` to `% WT` Expression if not already done (copy and paste the formula
    “`POWER(2,"_value_from_log2FC_cell")`” for the entire column). Save
    and import.

``` r
aDRG = read.csv("https://github.com/eriklarsen4/ggplot-scripts/blob/master/Bioinformatics/RNAseq%20Data%20Files/DESeq2%20Expression%20Results.csv")
  
  ## Filter out (subset) genes that went undetected
  ## or were outliers in terms of counts;
    ## new dataframe should not contain any NAs in p-value columns
aDRG3 = subset(aDRG, (!is.na(aDRG[,"AdjP"])))

 ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs.
aDRG9 = aDRG3 %>%
  filter(!grepl(aDRG3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))
#mt.+.?$|  <-- string identifier for mitochondrial tRNAs

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
aDRG9$labs = aDRG9$GeneID
  ## Replace DEGs with the string, "DEGs"
aDRG9$labs[which(aDRG9$AdjP <= 0.05)]= "DEGs"
  ## Replace the remaining genes with "Non-DEGs"
aDRG9$labs[which(aDRG9$AdjP >= 0.05)]= "Non-DEGs"

aDRG_DEG_list = c(aDRG9$GeneID[which(aDRG9$labs == "DEGs")])
```

## Data processing

How many genes are significantly different from WT expression?

``` r
length(
  which(aDRG9$AdjP <= 0.10) ) ## For FDR <= 0.1
```

    ## [1] 535

``` r
length(
  which(aDRG9$AdjP <= 0.05) ) ## For FDR <= 0.05
```

    ## [1] 376

``` r
length( 
  which(aDRG9$AdjP <= 0.01) ) ## For FDR <= 0.01
```

    ## [1] 205

Subset the dataframe to select differentially expressed genes with
Adjusted P-values &lt; 0.05 for use in pathway analysis tools (see [GO Analysis.md](https://github.com/eriklarsen4/Proteomics/blob/main/GO%20Analysis/GO-Analysis.md) to find potential mechanistic
pathways or other genes/interactions of interest).

``` r
HITS.01 = subset(aDRG9, aDRG9$AdjP <= 0.01)
HITS.05 = subset(aDRG9, aDRG9$AdjP <= 0.05)
HITS.10 = subset(aDRG9, aDRG9$AdjP <= 0.10)

  ## Create a new data frame of the average reads from both genotypes of genes with an Adjusted P-value < 0.05
Candidate_Screen = aDRG_TPM %>%
  filter(GeneID %in% HITS.05$GeneID)
```

## Volcano plot prep

Volcano plots display the `fold change` of the manipulation relative to
control (in this case, log<sub>2</sub>(FC) of mutant transcripts relative to
wild type) along the x-axis, and `Adjusted P-values` along the y-axis

Transform the `Adjusted P-values` so that they are integers and are
somewhat shrunk in a manner that is fit for viewing

``` r
aDRG9$log10ADJP = -log10(aDRG9$AdjP)
```

Create a new column that we will use for added visualization insight. Totally arbitrary, but we’ll call them
“genes of interest” or “`g.o.i.`” for short.

``` r
  ## Make the new column a column of characters.
    ## Easiest way is to just use ' "" ' or another column already populated by strings (characters)
aDRG9$g.o.i. = aDRG9$GeneID

    ## Alternatively, you can use dplyr like so (uncomment to execute the command):
# aDRG9 = aDRG9 %>%
#  dplyr::mutate(g.o.i. = GeneID)

  ## Fill the points with appropriately indexed data;
    ## "D.E.G.s", aka "Differentially Expressed Genes",
    ## defined as being "genes of interest" with Adjusted P-Values of <= 0.05
aDRG9$g.o.i.[which(aDRG9$AdjP <= 0.05)] = "D.E.G.s"

aDRG9$g.o.i.[which(aDRG9$AdjP > 0.05)] = "Non-D.E.G.s"

  ## dplyr alternative:
#aDRG9 = aDRG9 %>%
#  dplyr::mutate(g.o.i. = case_when(AdjP <= 0.05 ~ "D.E.G.s",
#                                    TRUE ~ "Non-D.E.G.s"))
```

Create labels on points we’re interested in labeling in the volcano
plots.

``` r
  ## Create the column, "labs"; hide all of the text labels with: ""
aDRG9$labs = ""
```

Label only these selected items; adjust for each experiment; for the
current example and a more global volcano plot, highlight the most
extreme data points on the plot.

``` r
ix_label1 = which(aDRG9$log10ADJP > 50)

  ## Fill in the labels
aDRG9$labs[c(ix_label1)] = aDRG9$GeneID[c(ix_label1)]
```

## First volcano plot

+ Use the `aDRG9` dataframe as the data to plot
+ Set the `x-axis` variable as log-base 2 fold change (`log2FC`)
+ Set the `y-axis` as log-base 10 Adjusted P-value (`log10ADJP`)
+ Set the points to be colored by the `g.o.i.` column previously created

Subset the data by each different group within that g.o.i. column

``` r
  ## Find axis limits
min(aDRG9$log2.FC.)
```

    ## [1] -3.578847

``` r
max(aDRG9$log2.FC.)
```

    ## [1] 1.455933

``` r
max(aDRG9$log10ADJP)
```

    ## [1] 127.9469

Create the plot, and save it as a variable. Add layers:

-   use the `ggplot` function to set the data
-   specify the type of graph, scatter plot points
-   subset to the background/bottom layer
-   map the `x`, `y` variables/axes
-   set `color` to `g.o.i.`
-   set `alpha` (transparency) low (closer to 0 than 1) for points that will be the bottom layer (get plotted over)
-   repeat for the upper layer, using a larger `alpha` than previous layers
-   set the coordinates of the graph with `coord_cartesian`
-   use `geom_hline` to create a horizontal line that will provide a
    clear division between `DEGs` and `Non-DEGs`
-   use `geom_text` to add a label to the line
-   use `labs` to add title, axes, and legend labels
-   use `theme_bw` to turn the graph border and background white with a
    black outline
-   use `theme` to control graph title and legend positions
-   use `scale_color_manual` to map color to the appropriate data points

``` r
vol1 = ggplot(data = aDRG9) + 
  geom_point(data = subset(aDRG9, g.o.i. == "Non-D.E.G.s"),
             aes(x = log2.FC., y = log10ADJP, color = g.o.i.), alpha = 0.1) +
  geom_point(data = subset(aDRG9, g.o.i. == "D.E.G.s"),
             aes(x = log2.FC., y = log10ADJP, color = g.o.i.), alpha = 0.3) +
  coord_cartesian(xlim = c(min(aDRG$log2.FC.),
                           max(aDRG9$log2.FC.)),
                  ylim = c(0, max(aDRG9$log10ADJP))) +
  geom_hline(yintercept = min(aDRG9$log10ADJP[which(aDRG9$g.o.i. == "D.E.G.s")]),
             linetype = "dashed",
             color = "firebrick") +
  geom_text(x = -3, y = 4,
            label = "FDR \U2264 0.05",
            color = "firebrick", size = 3) +
  labs(title = expression(paste(italic("Tmem184b")^"GT/GT", " DRG mRNA Exp. Changes")),
       x = "log2 (FC) rel. to WT",
       y = "-log10 (Adj.P-val)",
       col = "Gene Data Type") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.90,0.87)) + 
  scale_color_manual(values = c("darkgoldenrod4", "gray48"))
```

``` r
  ## Overwrite the variable, adding labels
vol1 = vol1 + geom_text_repel(data = aDRG9,
                              x = aDRG9$log2.FC.,
                              y = aDRG9$log10ADJP,
                              color = "black",
                              aes(label = aDRG9$labs),
                              max.overlaps = Inf) 
```

Plot the first volcano

``` r
  ## Plot the volcano
vol1
```

![](https://github.com/eriklarsen4/ggplot-scripts/blob/master/Bioinformatics/Plots/Broad%20aDRG%20Volcano-1.png)<!-- -->

## Second volcano plot

Adding a little more specific data..

``` r
  ## Re-set the mapping columm
aDRG9$g.o.i. = ""

  ## Fill the points with appropriately indexed data:
    ## Notable pruriceptor genes defined as being "genes of interest", with ADJPs of <= 0.05
aDRG9$g.o.i.[
  c(2,5,8,9,13,14,20,30,33,44,51,69,75,95,112,122,136,192,230,248,327,328,329,331)]= "Itch-related D.E.G.s"

    ## Fill the DEGs not related to itch
aDRG9$g.o.i.[c(1:376)][
  -c(2,5,8,9,13,14,20,30,33,44,51,69,75,95,112,122,136,192,230,248,327,328,329,331)] = "D.E.G.s"

    ## Fill the rest (Non-DEGs)
aDRG9$g.o.i.[which(aDRG9$AdjP > 0.05)] = "Non-D.E.G.s"
```

Add new labels

``` r
  ## Create a new "labs" column, though we could re-use the old one
aDRG9$labs2 = ""

aDRG9$labs2[c(2,5,8,9,13,20,30,51,112,331)] = aDRG9$GeneID[c(2,5,8,9,13,20,30,51,112,331)]
```

+ Store the new graph as a new variable
+ add a new subset
+ zoom in

``` r
vol2 = ggplot(data = aDRG9) + 
  geom_point(data = subset(aDRG9, g.o.i. == "Non-D.E.G.s"),
             aes(x = log2.FC., y = log10ADJP, color = g.o.i.), alpha = 0.1) +
  geom_point(data = subset(aDRG9, g.o.i. == "D.E.G.s"),
             aes(x = log2.FC., y = log10ADJP, color = g.o.i.), alpha = 0.3) +
  geom_point(data = subset(aDRG9, g.o.i. == "Itch-related D.E.G.s"),
             aes(x = log2.FC., y = log10ADJP, color = g.o.i.), alpha = 0.9) +
  coord_cartesian(xlim = c(-3.09,1.6), ylim = c(0,52)) +
  geom_hline(yintercept = min(aDRG9$log10ADJP[which(aDRG9$g.o.i. == "D.E.G.s")]),
             linetype = "dashed",
             color = "firebrick") +
  geom_text(x = -3, y = 2.5,
            label = "FDR \U2264 0.05",
            color = "firebrick", size = 3) +
  labs(title = expression(paste(italic("Tmem184b")^"GT/GT", " DRG mRNA Exp. Changes")),
       x = "log2 (FC) rel. to WT",
       y = "-log10 (Adj.P-val)",
       col = "Gene Data Type") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.85,0.85)) + 
  scale_color_manual(values = c("darkgoldenrod4", "navy", "gray48"))
```

``` r
  ## Overwrite the variable, adding labels
vol2 = vol2 + geom_text_repel(data = aDRG9,
                              x = aDRG9$log2.FC.,
                              y = aDRG9$log10ADJP,
                              color = "black",
                              aes(label = aDRG9$labs2), 
                              max.overlaps = Inf)
```

Plot the second volcano

``` r
  ## Plot the volcano
vol2
```

![](https://github.com/eriklarsen4/ggplot-scripts/blob/master/Bioinformatics/Plots/Itch%20Volcano-1.png)<!-- -->
