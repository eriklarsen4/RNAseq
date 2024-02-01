Itch paper bar plot
================
Erik Larsen
7/17/2021

I developed the following code snippet as a documentation of some of the
graphical analysis of `RNA-seq data` used for our publication,
[Transmembrane protein TMEM184B is necessary for interleukin-31–induced
itch](https://journals.lww.com/pain/Abstract/9000/Transmembrane_protein_TMEM184B_is_necessary_for.97918.aspx).
The data was processed on [Usegalaxy.org](https://usegalaxy.org/) using
a bioinformatics standard algorithm, known as
[DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).  
The bar plots are generated using the
[ggplot2](https://cran.r-project.org/package=ggplot2) package, and are
standard for `gene ontology` and `pathway` analyses.

Additional analyses, including `volcano plotting` and associated data
wrangling in other `R Markdown` files and scripts: [Gene Ontology and
Pathway Analysis R
Markdown](file:///M:/Erik/Data/Omics/Gene-Ontology-and-Pathway-Analysis-for-Lab.html),
companion [GO Analysis R
Script](M:/Erik/Data/Scripts/R/Specific/Searches/GO%20Analysis.R) R
script file, [Volcano plot tutorial R
Markdown](file:///M:/Erik/Data/Omics/Volcano-plot-tutorial.html),
companion [Volcano R
Script](M:/Erik/Data/Scripts/R/Specific/Graphs/Volcano%20Script.R).

## Environment Prep

Note that code can be sectioned and condensed with the `Alt + O`
command.

List of packages for this script:
[tidyverse](https://cran.r-project.org/package=tidyverse),
[plyr](https://cran.r-project.org/package=plyr),
[dplyr](https://cran.r-project.org/package=dplyr),
[reshape2](https://cran.r-project.org/package=reshape2),
[readr](https://cran.r-project.org/package=readr),
[stringr](https://cran.r-project.org/package=stringr),
[ggplot2](https://cran.r-project.org/package=ggplot2),
[ggrepel](https://cran.r-project.org/package=ggrepel),
[ggsignif](https://cran.r-project.org/package=ggsignif)

Load the packages.

``` r
  ## For data wrangling (slicing/adding/removing/melting/rearranging dataframes and their columns and rows):
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)

library(readr) ## For importing data and files

library(stringr) ## Awesome for manipulating strings

library(httpuv) ## For including plots in Markdown file output

library(ggplot2) ## For awesome, publication-quality graphs
library(ggrepel) ## For labeling mapped data on those plots
library(ggsignif) ## For adding significance stars onto plots

library(enrichR) ## Querying the enrichr db, which accesses the GO and PANTHER dbs
```

Import the relevant data.

``` r
e13_DRG = read.csv("https://github.com/eriklarsen4/ggplot-scripts/blob/master/Bioinformatics/RNAseq%20Data%20Files/DESeq2_result_file_on_GT_e13_PARA.csv",
                   header = FALSE)
```

## Data prep

For the `GSEA` data (`DESeq2`):
+ add appropriate column names  
+ remove pseudogenes and rRNAs  
+ subset genes that went undetected or were outliers in terms of
counts/reads

``` r
  ## Add column names
colnames(e13_DRG)[c(1:7)] = c("GeneID", "BaseMean", "log2FC", "stdErr", "WaldStats", "Pval", "AdjP")

  ## Remove pseudogenes or rRNAs
e13_DRG = e13_DRG %>% 
  filter(!grepl(e13_DRG$GeneID,
                pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

  ## Filter (subset) genes that went undetected or were outliers in terms of counts;
    ## Overwritten df should not contain any NAs in p-value columns
e13_DRG = subset(e13_DRG, (!is.na(e13_DRG[,"AdjP"])))
```

Access the `enrichr` server to store the `GO Biological Processes` and
`PANTHER Pathways` dbs as variables.

``` r
DBs = listEnrichrDbs()

  ## In this case, "GO_Biological_Process_2018" (index # 130), and "Panther_2016" (index # 102)
    ## Concatenate them in a list for subsetting, or slice directly
DBs = DBs$libraryName[c(130,102)]
```

Query the `GO Biological Processes` db.

Remove the `GO string identifiers`.  
Filter some of the processes out and **keep the top 19 processes**.  
Convert the `Adjusted P-values` to `-log10` scale.

``` r
e13_GOs = as.data.frame(
  enrichr(
    c(e13_DRG$GeneID[which(e13_DRG$AdjP < 0.01)]), DBs[1])
)
```

    ## Uploading data to Enrichr... Done.
    ##   Querying GO_Biological_Process_2018... Done.
    ## Parsing results... Done.

``` r
  ## Copy into a new object
GO_Bar = e13_GOs

  ## Remove the GO string identifiers
GO_Bar$GO_Biological_Process_2018.Term = GO_Bar$GO_Biological_Process_2018.Term %>%
  gsub(x = GO_Bar$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")
  ## Edit some columns and their data (overlap : genes in your list; 
  ## process size : total genes in the process)
GO_Bar = GO_Bar %>%
  separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
GO_Bar$Overlap = as.numeric(GO_Bar$Overlap)
GO_Bar$`Process Size` = as.numeric(GO_Bar$`Process Size`)

  ## Re-order the dataframe by Q value
GO_Bar = GO_Bar[order(GO_Bar$GO_Biological_Process_2018.Adjusted.P.value, decreasing = FALSE),]
  ## Re-set the index to make sure indeces are correct
row.names(GO_Bar) = NULL
  ## Convert AdjPs (Qvals) to -log10 scale
GO_Bar$GO_Biological_Process_2018.Adjusted.P.value = -log10(GO_Bar$GO_Biological_Process_2018.Adjusted.P.value)

  ## Remove unused columns and filter by the top BPs to graph
GO_Bar = GO_Bar[c(1:19) ,-c(7,8,9)]
```

Create a list of `GO BP Terms` that are more useful for graphing.

``` r
GO_Bar$GO_Biological_Process_2018.Term[c(13,19)] = c("regulation of mRNA splicing\n via spliceosome",
                                                     "RNA splicing, via transesterification reactions with\nbulged adenosine as nucleophile")
```

Repeat the process as above but for `PANTHER Pathways`.

``` r
e13_PATHWAYS = as.data.frame(
  enrichr(
    c(e13_DRG$GeneID[which(e13_DRG$AdjP < 0.01)]), DBs[2])
)
```

    ## Uploading data to Enrichr... Done.
    ##   Querying Panther_2016... Done.
    ## Parsing results... Done.

``` r
Pathway_Bar = e13_PATHWAYS

  ## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
Pathway_Bar$Panther_2016.Term = Pathway_Bar$Panther_2016.Term %>%
  gsub(x = Pathway_Bar$Panther_2016.Term,
       pattern = " Homo sapiens .+.?$",
       replacement = "")
  ## Repeat process for PATHWAYS
Pathway_Bar = Pathway_Bar %>%
  separate(Panther_2016.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
Pathway_Bar$Overlap = as.numeric(Pathway_Bar$Overlap)
Pathway_Bar$`Process Size` = as.numeric(Pathway_Bar$`Process Size`)

Pathway_Bar = add_column(Pathway_Bar,
                         EnrichrZscore = Pathway_Bar$Panther_2016.Combined.Score/log(Pathway_Bar$Panther_2016.P.value), .before = 4)
  ## Re-order the dataframe by Q value
Pathway_Bar = Pathway_Bar[order(Pathway_Bar$Panther_2016.Adjusted.P.value, decreasing = FALSE),]
  ## Re-set the row index to make sure indeces are correct
row.names(Pathway_Bar) = NULL
  ## Remove unused columns and filter by the top hits
Pathway_Bar = Pathway_Bar[c(1:10),-c(7,8,9)]

Pathway_Bar$Panther_2016.Adjusted.P.value = -log10(Pathway_Bar$Panther_2016.Adjusted.P.value)
Pathway_Bar = Pathway_Bar[c(1:10),]
```

## Plot the GO Biological Processes bar plot

Store the `ggplot` as a variable.  
+ **Order the processes by smallest Adjusted P-value to largest**,
simultaneously “mapping” the values of each term to a color (re-order
the x variable, y variable within the ggplot function).  
+ Adding the `geom_col()` function enables the plot to look like a bar
plot, and fills the aesthetics of this function completes the color
mapping.  
+ `scale_fill_continuous()` creates the color gradient and legend.  
+ `coord_flip()` switches axes for better visualization of the long GO
terms

``` r
GO_BAR_ADJP = ggplot(GO_Bar, 
                  aes(x = reorder(`GO_Biological_Process_2018.Term`, `GO_Biological_Process_2018.Adjusted.P.value`),
                      y = `GO_Biological_Process_2018.Adjusted.P.value`)) + 
  geom_col(stat = "identity" , aes(fill = `GO_Biological_Process_2018.Adjusted.P.value`)) + 
  scale_fill_continuous(name = "Enrichr\n-log10(Q)") + 
  labs(x = "GO Biological Process",
       y = "Enrichr (BHM) -log10(Adjusted P-value)",
       title = expression(paste("e13 ", italic("Tmem184b")^"GT/GT", " DRG GO Analysis"))) + 
  theme_light() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())  + 
  geom_signif(stat = "identity",
              aes(x = 19, y = 4.56, xend = 19, yend = 4.56, annotation = "***")) +
  geom_signif(stat = "identity",
              aes(x = 18, y = 4.56, xend = 18, yend = 4.56, annotation = "***")) +
  geom_signif(stat = "identity",
              aes(x = 17, y = 2.39, xend = 17, yend = 2.39, annotation = "**")) +
  geom_signif(stat = "identity",
              aes(x = 16, y = 2.39, xend = 16, yend = 2.39, annotation = "**")) +
  geom_signif(stat = "identity",
              aes(x = 15, y = 2.3, xend = 15, yend = 2.3, annotation = "**")) +
  geom_signif(stat = "identity",
              aes(x = 14, y = 2.3, xend = 14, yend = 2.3, annotation = "**")) +
  geom_signif(stat = "identity",
              aes(x = 13, y = 2.3, xend = 13, yend = 2.3, annotation = "**")) +
  geom_signif(stat = "identity",
              aes(x = 12, y = 2.3, xend = 12, yend = 2.3, annotation = "**")) +
  geom_signif(stat = "identity",
              aes(x = 11, y = 2.2, xend = 11, yend = 2.2, annotation = "**")) +
  geom_signif(stat = "identity",
              aes(x = 10, y = 2.2, xend = 10, yend = 2.2, annotation = "**")) +
  geom_signif(stat = "identity",
              aes(x = 9, y = 2.2, xend = 9, yend = 2.2, annotation = "**")) +
  geom_signif(stat = "identity",
              aes(x = 8, y = 1.8, xend = 8, yend = 1.8, annotation = "*")) +
  geom_signif(stat = "identity",
              aes(x = 7, y = 1.8, xend = 7, yend = 1.8, annotation = "*")) +
  geom_signif(stat = "identity",
              aes(x = 6, y = 1.8, xend = 6, yend = 1.8, annotation = "*")) +
  geom_signif(stat = "identity",
              aes(x = 5, y = 1.8, xend = 5, yend = 1.8, annotation = "*")) +
  geom_signif(stat = "identity",
              aes(x = 4, y = 1.8, xend = 4, yend = 1.8, annotation = "*")) +
  geom_signif(stat = "identity",
              aes(x = 3, y = 1.8, xend = 3, yend = 1.8, annotation = "*")) +
  geom_signif(stat = "identity",
              aes(x = 2, y = 1.8, xend = 2, yend = 1.8, annotation = "*")) +
  geom_signif(stat = "identity",
              aes(x = 1, y = 1.8, xend = 1, yend = 1.8, annotation = "*")) +
  geom_hline(yintercept = 1.30103, linetype = "solid", color = "firebrick") +
  geom_text(x = 2, y = 3, label = "FDR = 0.05", color = "firebrick", size = 4)+
  coord_flip()

  ## Plot the bar graph
GO_BAR_ADJP
```

![](https://github.com/eriklarsen4/ggplot-scripts/blob/master/Bioinformatics/Plots/GO%20Biological%20Processes%20Bar%20Plot-1.png)<!-- -->

## Plot the Panther Pathways bar plot

Repeat as above but for PANTHER Pathways.

``` r
PATHWAY_BAR_ADJP = ggplot(Pathway_Bar,
                          aes(x = reorder(Panther_2016.Term, Panther_2016.Adjusted.P.value),
                              y = Panther_2016.Adjusted.P.value)) +
  geom_col(stat = "identity" , aes(fill = Panther_2016.Adjusted.P.value)) +
  scale_fill_continuous(name = "Enrichr\n-log10(Q)") +
  labs(x = "Panther 2016 Pathways",
       y = "Enrichr (BHM) -log10(Adjusted P-value)",
       title = expression(paste("e13 ", italic("Tmem184b")^"GT/GT", " DRG Pathway Analysis"))) +
  theme_light() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank())  + 
  geom_hline(yintercept = 1.30103, linetype = "solid", color = "firebrick") +
  geom_text(x = 2, y = 1.1, label = "FDR =\n0.05", color = "firebrick", size = 3)+
  coord_flip()

  ## Plot the graph
PATHWAY_BAR_ADJP
```

![](https://github.com/eriklarsen4/ggplot-scripts/blob/master/Bioinformatics/Plots/Pathways%20Bar%20Blot%20by%20Qvalue-1.png)<!-- -->
