

##### Heatmaps ####


library(readr)
library(tidyverse)
library(ggplot)
library(reshape2)
library(dplyr)
library(plyr)

##### e13 Panther GO Bio Processes Heatmap #####

## Import the relevant data
e13_DEGs_and_Pathways =  read_csv("C:/Users/Erik/Desktop/BoxCopy/Lab/Omics/RNAseq/Embryonic DRG/Analysis/Heatmaps/Custom Python Enrichr Pathway Clustergram e13 GT.csv")

## Export the gene list for photoshopping
#write.csv(e13_DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\e13 FDR 01 Heatmap DEGs.csv", sep = ",")


## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
e13_Heatmap = melt(e13_DEGs_and_Pathways, id = "DEGs")
colnames(e13_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

## Remove excessive pathway name strings
Pathway_Names = e13_Heatmap$Pathways %>% gsub(x = e13_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

## Put the new strings back into the dataframe
e13_Heatmap[,2] = Pathway_Names

## Make the digital values factors so it's compatible to graph
e13_Heatmap$`Presence In Pathway` = as.factor(e13_Heatmap$`Presence In Pathway`)


Term_list = c("Axon guidance mediated by Slit/Robo", "De novo pyrimidine\ndeoxribonucleotide biosynthesis", "Dopamine receptor mediated\nsignaling pathway", "FAS signaling pathway", "Heme biosynthesis", "Integrin signaling pathway", "Muscarinic acetylcholine\nreceptor 2 and 4 signaling pathway", "Oxidative stress response", "Synaptic vesicle trafficking", "Wnt signaling pathway")

## Create a for loop to replace the strings
for (i in 1:length(Term_list)){
  Pathway_Names[ which(Pathway_Names == unique(Pathway_Names)[i]) ] = Term_list[i]
}
e13_Heatmap[,2] = Pathway_Names

## Create and store the heatmap core as a variable
## Use geom_tile to map the binary values


## Customize by removing the filler surrounding the graph and the tickmarks
## Center the Title

## Use geom_tile to map the binary values

HEATER = ggplot(e13_Heatmap, aes(e13_Heatmap$Pathways, e13_Heatmap$DEGs)) +
  geom_tile(aes(fill = e13_Heatmap$`Presence In Pathway`), color = "black") +
  scale_fill_manual(values = c("dark blue", "gold")) +
  labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes\n(FDR < 0.01)", title = expression(paste("Heatmap of ", italic("Tmem184b")^"GT/GT","affected Pathways"))) +
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
