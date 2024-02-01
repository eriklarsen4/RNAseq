Binomial Heatmap Script
================
Erik Larsen
7/17/2020

The following code was developed to display which differentially
expressed genes are involved in a given number of `PANTHER` pathways
using `ggplot2`.  
The data was processed using `Salmon`, `DESeq2`, and run through another
script in `Python` to acquire the proper data structure for `ggplot2` to
handle below in a `geom_tile` call.  
Other scripts use the `pheatmap` package for more traditional
transcriptional profiling.

## Environment Prep

Note that code can be sectioned and condensed with the `Alt + O`
command.

List of packages for this script:
[tidyverse](https://cran.r-project.org/package=tidyverse),
[stringr](https://cran.r-project.org/package=stringr),
[reshape2](https://cran.r-project.org/package=reshape2),
[readr](https://cran.r-project.org/package=readr),
[ggplot2](https://cran.r-project.org/package=ggplot2)

``` r
library(readr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringr)
```

``` r
  ## Import the relevant data
e13_DEGs_and_Pathways =  read_csv("https://github.com/eriklarsen4/ggplot-scripts/blob/master/Binomial-Heatmap/Custom%20Python%20Enrichr%20Pathway%20Clustergram%20e13%20GT.csv")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
e13_Heatmap = melt(e13_DEGs_and_Pathways, id = "DEGs")
colnames(e13_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

  ## Remove excessive pathway name strings
Pathway_Names = e13_Heatmap$Pathways %>%
  gsub(x = e13_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

  ## Re-shape the pathways for aesthetic viewing
Term_list = c("Axon guidance mediated by Slit/Robo", "De novo pyrimidine\ndeoxribonucleotide biosynthesis",
              "Dopamine receptor mediated\nsignaling pathway", "FAS signaling pathway", "Heme biosynthesis",
              "Integrin signaling pathway", "Muscarinic acetylcholine\nreceptor 2 and 4 signaling pathway",
              "Oxidative stress response", "Synaptic vesicle trafficking", "Wnt signaling pathway")

  ## Create a for loop to replace the strings
for (i in 1:length(Term_list)){
  Pathway_Names[ which(Pathway_Names == unique(Pathway_Names)[i]) ] = Term_list[i]
}
  ## Put the new strings back into the dataframe
e13_Heatmap[,2] = Pathway_Names

  ## Make the digital values factors so it's compatible to graph
e13_Heatmap$`Presence In Pathway` = as.factor(e13_Heatmap$`Presence In Pathway`)
```

## Create and store the heatmap

(There are 84 genes involved, so for aesthetic readability, they have
been omitted from the plot)

``` r
## Use geom_tile to map the binary values
    ## Customize by removing the filler surrounding the graph and the tickmarks
      ## Center the Title

HEATER = ggplot(e13_Heatmap, aes(e13_Heatmap$Pathways, e13_Heatmap$DEGs)) +
  geom_tile(aes(fill = e13_Heatmap$`Presence In Pathway`), color = "black") +
  scale_fill_manual(values = c("dark blue", "gold")) +
  labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes\n(FDR < 0.01)",
       title = expression(paste("Heatmap of ", italic("Tmem184b")^"GT/GT","-affected Pathways"))) +
  guides(fill = guide_legend(title = "Presence in Pathway")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        panel.grid.major = element_line(colour = "black", size = 0.1),
        panel.grid.minor = element_line(colour = "black", size = 0.1),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0)) +
  coord_flip()

## Plot the heatmap
HEATER
```

![](https://github.com/eriklarsen4/ggplot-scripts/blob/master/Binomial-Heatmap/Heatmap%20of%20DEGs%20Identified%20Across%20PANTHER%20Pathways%20wLabs-1.png)<!-- -->
