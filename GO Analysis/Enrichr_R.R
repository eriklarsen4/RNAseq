#### This script will use the R-version of EnrichR to extract information about gene lists from online tools such as pathways####
#### It is formatted currently to use Panther 2016 - needs to be updated.#####
#### Last modified 4-29-20 by Martha Bhattacharya###########

install.packages("enrichR")
library(enrichR)
library(tidyverse)
library(ggplot2)
library(purrr)
library(reshape2)
library(dplyr)


##### STOP. You will need to load a file containing the Adj P < (your preferred threshold) in the format that DeSeq2 provides. ######
#### Here I am using Adj P < 0.01 as my cutoff, and have already loaded it into R. ########

HITS.01$GeneID <- as.character(HITS.01$GeneID)
str(HITS.01)

#gives a list of functions in the Enrichr databases
dbs <- listEnrichrDbs()

#find the ones that have panther
panther_rows <- grep("Panther", dbs$libraryName)
list(dbs$libraryName[panther_rows])
#this reveals that there is a 2015 and 2016 
 "ftp://ftp.pantherdb.org/panther_library/current_release/PANTHER15.0_hmmscoring.tgz"
 
####### This is an attempt to download the latest Panther version and make EnrichR use it. Not working yet. ########## 

 panther15 <- "ftp://ftp.pantherdb.org/hmm_classifications/current_release/PANTHER15.0_HMM_classifications"
 download.file(panther15,destfile="N:/Martha/Big Data from others/Panther 15.0/panther15")
 untar("panther15.tar.gz",list=TRUE)  ## check contents
 untar("panther15.tar.gz")
 
 if (!requireNamespace("AnnotationHub")) BiocManager::install("AnnotationHub")
 library(AnnotationHub)
 ah <- AnnotationHub()
 query(ah, "PANTHER.db")[[1]]
 
 library(PANTHER.db)
 PANTHER.db
 #restrict to mouse
 pthOrganisms(PANTHER.db) <- "MOUSE"

############## Proceed with Panther 2016 version ################
 
#we'll just use the 2016 version; make a variable for this character string and run
dbsp <- "Panther_2016"
if (websiteLive) 
  pant2016_enriched <- enrichr(HITS.01$GeneID, dbsp)

pant2016_enriched <- as.data.frame(pant2016_enriched)
View (pant2016_enriched)
pant_genelist <- data.frame(cbind(pant2016_enriched$Panther_2016.Term, pant2016_enriched$Panther_2016.Genes))
colnames(pant_genelist) = cbind("Pathway", "Genes In List")
View(pant_genelist)
write.csv(pant2016_enriched, file = "C:\\Users\\marthab1\\Documents\\R\\Volcanos and Heatmaps\\pant2016_enriched.csv", col.names = T)

## Make a data frame of just the pathway, p value, combined score, and genes
PantherDF = subset(pant2016_enriched[c(1,3,8,9)])


## The "str_split" function is part of the "stringr" package which seems to be part of the "tidyverse" package
## The following command will detect semicolons and replace them with "empty" spaces, while retaining the structure and organization of your vector/dataframe
PantherDF$Panther_2016.Genes = str_split(string = PantherDF$Panther_2016.Genes, pattern = ";", simplify = FALSE)

## Clean the dataframe to include the pathway/Panther term and the associated genes
PantherDF = subset(PantherDF[,c(1,2,3,4,8,9)])

class(PantherDF)
str(PantherDF)

genes <- PantherDF$Panther_2016.Genes

genelist <- read.csv("M://Omics//e13_panther_genes.csv")
str(genelist)
genedf <- genelist[3:nrow(genelist),]
genelist$Panther.Term[[2]]


############### Things to try but are not working yet, below ########################333

#okay, so the list of things to tackle include:
  #1 - figure out the rank list of genes appearing in the lists
  #2 - make graph of names vs # appearances (like a histogram)
  #3 (separate) make a 

library(stringi)
#for a given gene, this command gives you whether it exists in the listed row
stri_detect_fixed(PantherDF$Panther_2016.Genes, "SNAP25")
?stri_detect()
#this command gives a count of the number of genes in each list!
stri_count_boundaries(PantherDF$Panther_2016.Genes)

list2 <- stri_extract_all_words(PantherDF$Panther_2016.Genes)

#find unique genes in the list
genesonly <- PantherDF$Panther_2016.Genes
alluniquegenes <- unique(as.data.frame(genelist))

#might be able to use this to find the maximum represented string in the list:
genesonly <- PantherDF$Panther_2016.Genes
topgene<- genesonly.groupBy(identity).maxBy(genesonly.genesonly2.size).genesonly1
tbl %>% flatten_chr %>% unique
