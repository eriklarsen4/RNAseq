Gene Ontology Analysis Snippet
================
Erik Larsen
9/17/2021

# Overview

This markdown serves as documentation of gene ontology analysis of
multiple datasets generated in the lab of **Dr. Martha Bhattacharya**,
two of which have been used for publications
([here](https://journals.lww.com/pain/abstract/2022/05000/transmembrane_protein_tmem184b_is_necessary_for.18.aspx),
[here](https://www.biorxiv.org/content/10.1101/2022.05.27.493793v3#:~:text=Taken%20together%2C%20our%20data%20suggest%20that%20TMEM184B%20is,be%20linked%20to%20neurodevelopmental%20disorders%20than%20to%20dementia.)).

-   `Gene ontology` and `pathway enrichment analyses` provide high-level
    analysis of the cell/mol consequences following perturbations
    (genetic or pharmacological) of interest. This is incredibly
    powerful for studying novel proteins.

    -   As such, the approach has a few **problems**:

        1.  When analyzing a list of genes or proteins quantitatively
            changed by a perturbation, implicated biological processes
            or pathways depend on a database’s curation of data and the
            user’s knowledge of biological processes or pathways.

        2.  Genes and proteins have aliases; these can cause
            discrepancies missed by researchers querying gene ontologies
            when their supplied list does not contain the “correct”
            alias

            -   i.e. `GO Terms` may be associated with some
                genes/proteins but not all aliases

        3.  Some genes/proteins may have been mistakenly omitted from a
            `GO Term`

        4.  It’s probably useful to understand putative connections
            between a given gene/protein & its gene ontologies with
            *other* genes/proteins & *their* gene ontologies— in
            particular, accounting for as many aliases as possible

            -   If a user is interested in a manipulation’s effect on
                cellular processes, this approach may help better
                connect otherwise unconnected genes/proteins to other
                pathways/processes

    -   Thus, I created multiple functions, harnessing the
        `GO database`’s repository, but with more flexibility than just
        plugging and chugging into `geneontology.org`’s or `enrichr`’s
        websites

    -   This snippet demonstrates the creation and use of one particular
        function (another that addresses issues is documented [here]()):

        -   `GO_INFO_fn`, for performing the statistical analyses of the
            genes/proteins in a provided list

            -   More specific questions and their follow-up functions to
                the `GO_INFO_fn` are specified in the
                `Additional Functions` tab of the `Function definitions`
                section below

# Function definitions

## GO_INFO_fn

-   For the GO and pathway over-representation analysis of a list of
    interest (e.g. analysis of differentially expressed genes)

-   The `GO_INFO_fn` takes, as input:

    -   a list of gene or protein IDs
    -   the species of the list

-   it returns to the Global Environment:

    -   a dataframe, `GO_INFO_by_GENE_df`, that houses:

        -   each gene/protein from a provided list (the input list)
            -   all the `# of GO Terms` to which a gene/protein has been
                annotated by the `GO.db`
            -   all the `GO IDs` to which a gene/protein has been
                annotated by the `GO.db`
            -   all the `GO Terms` to which a gene/protein has been
                annotated by the `GO.db`

    -   a list (`Unique_GOs`) of all the unique `GO Terms` associated
        with the entire list

    -   the supplied list of interest, removing genes/proteins that do
        not have terms associated with it ( = `list_of_interest`)

## Additional Functions

The `GO_INFO_by_GENE_df` houses information that can then be used for
more targeted questions, using other functions, such as:

-   Are there any enriched `GO Terms` from the supplied list?

    -   if yes, what are they? ( use
        `Find_Enriched_GOs_from_GO_INFO_df_fn`)

    -   Which genes/proteins from the list of interest are associated
        with a specific `GO Term`? ( use `Find_Genes_With_X_GO_Term_fn`)

    -   Which `GO Terms` are associated with a specific gene/protein
        from the supplied list? ( use `Find_GOs_of_Gene_X_fn`)

    -   What are the `GO Terms` of two genes/proteins of from the
        supplied list? ( use `Find_GOs_of_Two_Genes_fn`)

    -   Which `GO Terms` are shared by all genes/proteins from the
        supplied list? ( use `Find_Genes_Related_By_GO_Term_fn`)

# Example Environment Prep

## Load packages and data

Load packages

``` r
library(stringr)
library(ggplot2)
library(ggrepel)
# library(enrichR)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(GO.db)
library(tidyverse)
# library(tidyr)
library(ggpubr)
library(mgcv)
library(splines)
```

Import data (not shown)

``` r
aDRG = read_csv("~/GitHub/ggplot-scripts/Bioinformatics/RNAseq Data Files/DESeq2 Expression Results.csv")

  ## Filter (subset) genes that went undetected or were outliers in terms of counts;
  ## new dataframe should not contain any NAs in p-value columns
aDRG_sub = subset(aDRG, (!is.na(aDRG[,"AdjP"])))
 ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs, along with pseudogenes, etc.
aDRG_sub = aDRG_sub %>% filter(!grepl(GeneID,
                                pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))
#mt.+.?$|  <-- string identifier for mitochondrial tRNAs

  ## Create a column in the DESeq2 dataframe that scales the Adjusted P-value by log-base 10
aDRG_sub$log10ADJP = -log10(aDRG_sub$AdjP)

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
aDRG_sub = aDRG_sub %>%
  dplyr::mutate(log10ADJP = -log10(AdjP),
                g.o.i. = GeneID,
                g.o.i. = case_when(AdjP <= 0.05 ~ "DEGs",
                                   TRUE ~ "Non-DEGs"),
                labs = case_when(GeneID == "Tmem184b" ~ "Tmem184b",
                                 TRUE ~ ""))

aDRG_DEG_list = aDRG_sub %>%
  dplyr::filter(g.o.i. == "DEGs") %>%
  dplyr::select(GeneID) %>%
  as.vector() %>%
  unlist() %>%
  as.character()
```

Source functions

``` r
source("~/GitHub/RNAseq/Functions/GO Analysis Functions.R")
```

## Visualize Dataset

Visualize the filtered dataset, highlighting the differentially
expressed genes (gold)

**The DEGs are the genes used as the list in this analysis**

![](GO-Analysis_files/figure-gfm/Plot%20the%20Volcano-1.png)<!-- -->

# GO Analysis on List of Interest

Conduct GO Analysis on the example list above using the `GO_INFO_fn`:

-   queries the `GO.db` and extracts `Terms` and `IDs` associated with
    genes/proteins in a provided list

-   output is `GENE_GO_INFO_df` (gene-wise GO information) and
    `GO_INFO_by_TERM_df` (term-wise GO information) in the Global
    Environment

One approach to GO Analysis is by determining `GO Term` or `Pathway`
`Over-representation`, meaning:

-   finding the probability that a particular `GO Term`/`Pathway` is
    more likely than chance to be over- or under-represented in the
    provided list

-   To find this probability, we first need to:

    -   determine the probability of observing the genes/proteins from
        the supplied list in a given `Term` (for all `GO Term`s)

    -   these probabilities (“`Overlap_Prob`s”) are determined using
        **Bayes’ Theorem**, where:

        -   `Overlap_Prob` = n<sub>genes</sub> from `list of interest`
            in the given `GO Term` / total_number<sub>genes</sub> in
            that `GO Term`

    -   “`Weighted_Overlap_Prob`s” are computed in an effort to control
        for list size, heavily annotated genes, and term size, where:

        -   `Weighted_Overlap_Prob` = ( (n<sub>genes</sub> from
            `list of interest` in the given `GO Term`) \*
            (`GO_Term_rate`) \* (`Gene_freq` ) ) / (`list length` /
            `GO Term size`)

        -   `GO_Term_rate` = n<sub>unique_GO_Terms</sub> from
            `list of interest` / n<sub>total</sub> `GO`s

            -   (accounts for the number of `GO Terms` related to the
                provided list)

        -   `Gene_freq` = sum( P<sub>r</sub> (each gene from
            `list of interest` arising in *any* `GO Term`) )

            -   (accounts for heavily annotated genes i.e. **how often
                does gene X get annotated to a different GO Term?**)

        -   `GO Term size` = n<sub>total_genes</sub> in a `GO Term`

-   Finally, a binomial exact test is conducted for each `GO Term` to
    assess how likely the observed probability of overlap is when
    compared to a theoretical probability

## GO_INFO_fn

``` r
GO_INFO_fn(list_of_interest = aDRG_DEG_list,
           species = "mouse")
```

## GO Overlap

Find the `Overlap_Prob`s for the list of interest

-   To determine these probabilities (“`Overlap_Prob`s”), use **Bayes’
    Theorem**

-   Assess their distributions, make corrections

``` r
  ## Store GO term information
for (i in 1:nrow(GO_INFO_by_TERM_df) ) {
  
  x = paste(
    unique(
      unlist(
        str_split(as.vector(GO_INFO_by_TERM_df[i, 4]),
                  pattern = ";")))
    , collapse = ";")
  
  GO_INFO_by_TERM_df[i,3] = length(c(unique(unlist(str_split(as.vector(GO_INFO_by_TERM_df[i, 4]),
                                                             pattern = ";")
                                                   )
                                            )
                                     )
                                   )
  GO_INFO_by_TERM_df[i,4] = x
  
}

  ## Store conditional probabilities
GO_INFO_by_TERM_df = GO_INFO_by_TERM_df %>%
  dplyr::mutate(Num_Genes_in_Term = as.numeric(Num_Genes),
                
                Num_Genes_in_Term_from_List = as.numeric(Num_Genes_from_List),
                
                Num_Genes_not_in_Term = max(Num_Genes_in_Term) - Num_Genes_in_Term,
                
                Num_Genes_not_in_Term_from_List = as.numeric(length(list_of_interest)) - Num_Genes_in_Term_from_List,
                
                Num_Genes_in_Term_not_from_List = Num_Genes_in_Term - Num_Genes_in_Term_from_List,
                
                Num_Genes_not_from_List = max(Num_Genes_in_Term) - as.numeric(length(list_of_interest)),
                
                Num_Genes_not_in_Term_not_from_List = Num_Genes_not_from_List - Num_Genes_in_Term_not_from_List) %>%
  
  dplyr::select(-Num_Genes, -Num_Genes_from_List) %>%
  dplyr::arrange(desc(Num_Genes_in_Term))
```

To find over-representation, determine the probabilities that genes from
the `list of interest` are found in a given `GO Term`:

-   `Overlap_Prob` = n<sub>genes</sub> from `list of interest` in the
    given `GO Term` / tot_n<sub>genes</sub> in that `GO Term`

-   `Weighted_Overlap_Prob` = ( (n<sub>genes</sub> from
    `list of interest` in the given `GO Term`) \* (`GO_Term_rate`) \*
    (`Gene_freq` ) ) / (`list length` / `GO Term size`), where

    -   `GO_Term_rate` = n<sub>unique_GO_Terms</sub> from
        `list of interest` / n<sub>total</sub> `GO`s

    -   `Gene_freq` = sum( P<sub>r</sub> (each gene from
        `list of interest` arising in *any* `GO Term`) )

        -   i.e. **how often does gene X get annotated to a different GO
            Term?**

    -   `GO Term size` = n<sub>total_genes</sub> in a `GO Term`

``` r
  ## Initialize variables for generating the weighted overlap probabilities for each term / gene

  ## how frequently they occur across all GO terms
gene_freq = list()

  ## overlap variables

    ## (aka how many genes in the list of interest were detected in GO Term X? --> Overlap)
    ## for weighted:
      ## account for: 
        ## list size
        ## term size
        ## number of unique terms found from list of interest
        ## how frequently genes from the list are found in any GO terms
GO_INFO_by_TERM_df = GO_INFO_by_TERM_df %>%
  dplyr::mutate(Overlap_Prob = (Num_Genes_in_Term_from_List / max(Num_Genes_in_Term)
                                ),
                Weighted_Overlap_Prob = NA_real_)

  ## Populate
for (i in 1:length(list_of_interest)){
  gene_freq[[i]] = 
    (
      as.numeric(GENE_GO_INFO_df$Num_GO_Terms[which(GENE_GO_INFO_df$GeneID == list_of_interest[i])]) /
        42837
      )
}
names(gene_freq) = c(list_of_interest)

for (i in 1:nrow(GO_INFO_by_TERM_df) ) {
  for (j in 1:length(GO_INFO_by_TERM_df$Num_Genes_in_Term[i]) ) {
    
     GO_INFO_by_TERM_df$Weighted_Overlap_Prob[i] = 
       
       ## overlap rate relative to list size
       (GO_INFO_by_TERM_df$Num_Genes_in_Term_from_List[i]/as.numeric(length(list_of_interest))*
          ## gene annotation frequency
          as.numeric(unlist(gene_freq))[j]*
          ## GO Term frequency
          (as.numeric(length(Unique_GOs))/42837))/
       ## list size relative to the genes in the GO Term
       (as.numeric(length(list_of_interest))/GO_INFO_by_TERM_df$Num_Genes_in_Term[i])
  }
}
```

## Plot Probability Distributions

Plot the `Weighted Overlap Prob` distribution and (unweighted)
`Overlap Prob` distribution

-   distributions of probabilities that given genes/proteins from a
    provided list are all identified within a given `GO Term`
-   `weighted` accounts for list size, term size, frequency of a gene
    occurring in any term

``` r
  ## Inspect the probability distributions
GO_INFO_by_TERM_df %>%
  dplyr::select(GO_Term, Overlap_Prob, Weighted_Overlap_Prob) %>%
  pivot_longer(cols = c(contains("Overlap")),
               names_to = "Overlap Prob Type",
               values_to = "Prob") %>%
  ggplot(aes(x = Prob, fill = `Overlap Prob Type`)) +
  geom_histogram(aes(x = Prob))
```

![](GO-Analysis_files/figure-gfm/Plot%20Overlap%20Probability%20Distributions-1.png)<!-- -->

Weighting broadens the binomial distribution of probabilities away from
zero

# Hypothesis Test

Are the `Overlap Probabilities` significantly different than we would
expect by chance?

-   Perform a statistical test on each term’s `Overlap Probability`:

    -   binomial exact test compares:

        1.  the probability of the number of genes identified from the
            provided list involved in a particular `GO Term`
            (`Overlap Prob`)

        2.  the theoretical (`Overlap Prob`) probability of the number
            of genes involved in a particular `GO Term`

            -   theoretical probability is determined by:
                -   simulating 1000 draws from a binomial distribution
                    for a `GO Term` of the same size with the
                    probability determined from the provided list
                    (above)
                -   averaging those 1000 numbers of simulated
                    “successes”
                -   dividing by the maximum number of genes in all
                    `GO Term`s

    -   test returns a p-value of a two-sided binomial exact test

-   It is [known](https://geneontology.org/docs/faq/) that statistical
    results from this analysis may not be meaningful to the user because
    this approach will return large, vague `GO Terms` OR incredibly
    noisy and specific `GO Terms`, with a heavy dependence on the list
    size

## Binomial Exact Test 1

Perform an exact test across the derived probabilities for each
`GO Term` to determine if that term is over- or under- represented
within the provided list:

-   `Weighted Overlap Probs` vs `Expected`

-   `Unweighted Overlap Probs` vs `Expected`

``` r
  ## Create variables for hypothetical probabilities, drawing from binomial distribution and performing an exact test
GO_INFO_by_TERM_df = GO_INFO_by_TERM_df %>%
  dplyr::mutate(E_Weighted_Overlap = NA_real_,
                E_Weighted_Overlap_Prob = NA_real_,
                Weighted_Overlap_Prob_p = NA_real_,
                E_Overlap = NA_real_,
                E_Overlap_Prob = NA_real_,
                Overlap_Prob_p = NA_real_)

  ## Fill variables
for( i in 1:nrow(GO_INFO_by_TERM_df) ) {
  
  ## what's the expected weighted overlap? probability? -> 1000 draws
    ## take the average of the 1000 draws for the overlap
      ## re-scale by factor of 10
  GO_INFO_by_TERM_df$E_Weighted_Overlap[i] = round(mean(
    rbinom(1000,
           size = GO_INFO_by_TERM_df$Num_Genes_in_Term[i],
           prob = GO_INFO_by_TERM_df$Weighted_Overlap_Prob[i]
           )
    ))
  
  ## same but un-weighted
  GO_INFO_by_TERM_df$E_Overlap[i] = round(mean(
    rbinom(n = 1000,
           size = GO_INFO_by_TERM_df$Num_Genes_in_Term[i],
           prob = GO_INFO_by_TERM_df$Overlap_Prob[i]
           )
    ))
}
  
    ## divide it by the total number of genes
GO_INFO_by_TERM_df$E_Weighted_Overlap_Prob = 
  
  GO_INFO_by_TERM_df$E_Weighted_Overlap / max(GO_INFO_by_TERM_df$Num_Genes_in_Term)


GO_INFO_by_TERM_df$E_Overlap_Prob = 
  
  GO_INFO_by_TERM_df$E_Overlap / max(GO_INFO_by_TERM_df$Num_Genes_in_Term)

for( i in 1:nrow(GO_INFO_by_TERM_df)) {
    ## perform an exact binomial test and extract p-values
  GO_INFO_by_TERM_df$Weighted_Overlap_Prob_p[i] =
    
    binom.test(GO_INFO_by_TERM_df$Num_Genes_in_Term_from_List[i],
               GO_INFO_by_TERM_df$Num_Genes_in_Term[i],
               (GO_INFO_by_TERM_df$E_Weighted_Overlap_Prob[i]),
               alternative = "two.sided")$p.value
  
  
  GO_INFO_by_TERM_df$Overlap_Prob_p[i] = 
    
    binom.test(GO_INFO_by_TERM_df$Num_Genes_in_Term_from_List[i],
                GO_INFO_by_TERM_df$Num_Genes_in_Term[i],
                (GO_INFO_by_TERM_df$E_Overlap_Prob[i]),
                alternative = "two.sided")$p.value
}
 
  ## any significant? remove low count (1 hit) data and add adjusted p-values
GO_INFO_by_TERM_df %>%
  dplyr::arrange(Weighted_Overlap_Prob) %>%
  dplyr::mutate(Weighted_Overlap_Prob_ADJP = (Weighted_Overlap_Prob_p / nrow(GO_INFO_by_TERM_df)) * 0.05,
                 Overlap_Prob_ADJP = (Overlap_Prob_p / nrow(GO_INFO_by_TERM_df)) * 0.05) %>%
  dplyr::filter(Num_Genes_in_Term_from_List > 1) %>%
  dplyr::slice(c(1:15)) %>%
  dplyr::select(c(GO_Term, Num_Genes_in_Term_from_List, Num_Genes_in_Term, contains("Overlap")))
```

    ##                                                                                GO_Term
    ## 1                                                   actin filament bundle distribution
    ## 2                                        heparan sulfate 6-O-sulfotransferase activity
    ## 3           regulation of CDP-diacylglycerol-serine O-phosphatidyltransferase activity
    ## 4  positive regulation of CDP-diacylglycerol-serine O-phosphatidyltransferase activity
    ## 5                        positive regulation of serine C-palmitoyltransferase activity
    ## 6                          intrinsic component of neuronal dense core vesicle membrane
    ## 7                                   intrinsic component of dense core granule membrane
    ## 8                           integral component of neuronal dense core vesicle membrane
    ## 9                                               alpha-1,2-mannosyltransferase activity
    ## 10                                                           GTPase inhibitor activity
    ## 11                                             oncostatin-M-mediated signaling pathway
    ## 12                                                     mu-type opioid receptor binding
    ## 13                                                                   kininogen binding
    ## 14                                                                             C-fiber
    ## 15                                                 modulation by virus of host process
    ##    Num_Genes_in_Term_from_List Num_Genes_in_Term Overlap_Prob
    ## 1                            2                 2 6.910134e-05
    ## 2                            2                 3 6.910134e-05
    ## 3                            2                 3 6.910134e-05
    ## 4                            2                 3 6.910134e-05
    ## 5                            2                 3 6.910134e-05
    ## 6                            2                 3 6.910134e-05
    ## 7                            2                 3 6.910134e-05
    ## 8                            2                 3 6.910134e-05
    ## 9                            2                 4 6.910134e-05
    ## 10                           2                 4 6.910134e-05
    ## 11                           2                 4 6.910134e-05
    ## 12                           2                 4 6.910134e-05
    ## 13                           2                 4 6.910134e-05
    ## 14                           2                 4 6.910134e-05
    ## 15                           2                 4 6.910134e-05
    ##    Weighted_Overlap_Prob E_Weighted_Overlap E_Weighted_Overlap_Prob
    ## 1           1.144145e-09                  0                       0
    ## 2           1.716217e-09                  0                       0
    ## 3           1.716217e-09                  0                       0
    ## 4           1.716217e-09                  0                       0
    ## 5           1.716217e-09                  0                       0
    ## 6           1.716217e-09                  0                       0
    ## 7           1.716217e-09                  0                       0
    ## 8           1.716217e-09                  0                       0
    ## 9           2.288290e-09                  0                       0
    ## 10          2.288290e-09                  0                       0
    ## 11          2.288290e-09                  0                       0
    ## 12          2.288290e-09                  0                       0
    ## 13          2.288290e-09                  0                       0
    ## 14          2.288290e-09                  0                       0
    ## 15          2.288290e-09                  0                       0
    ##    Weighted_Overlap_Prob_p E_Overlap E_Overlap_Prob Overlap_Prob_p
    ## 1                        0         0              0              0
    ## 2                        0         0              0              0
    ## 3                        0         0              0              0
    ## 4                        0         0              0              0
    ## 5                        0         0              0              0
    ## 6                        0         0              0              0
    ## 7                        0         0              0              0
    ## 8                        0         0              0              0
    ## 9                        0         0              0              0
    ## 10                       0         0              0              0
    ## 11                       0         0              0              0
    ## 12                       0         0              0              0
    ## 13                       0         0              0              0
    ## 14                       0         0              0              0
    ## 15                       0         0              0              0
    ##    Weighted_Overlap_Prob_ADJP Overlap_Prob_ADJP
    ## 1                           0                 0
    ## 2                           0                 0
    ## 3                           0                 0
    ## 4                           0                 0
    ## 5                           0                 0
    ## 6                           0                 0
    ## 7                           0                 0
    ## 8                           0                 0
    ## 9                           0                 0
    ## 10                          0                 0
    ## 11                          0                 0
    ## 12                          0                 0
    ## 13                          0                 0
    ## 14                          0                 0
    ## 15                          0                 0

Results are plausible, given what we suspect about this dataset

Sorting by p-value returns a flood of 0 p-values (unhelpful)

## Plot Expected Probs vs Actual Probs

Plot the `Expected Overlap Probabilities` (x axis) and the
`Overlap Probabilities` (y axis)

-   `weighted overlap`

-   `unweighted overlap`

-   include a generalized additive model fit for each

![](GO-Analysis_files/figure-gfm/Prob%20plots%202-1.png)<!-- -->

## Assess GAM Model fits

How does the weighted overlap compare?

``` r
GO_INFO_by_TERM_df %>%
  dplyr::mutate(
                Weighted_Overlap_RMSE = sqrt(mean(weighted_gam_preds - Weighted_Overlap_Prob)^2),
                
                Overlap_RMSE = sqrt(mean(gam_preds - Overlap_Prob)^2),
                
                Weighted_gam_Preds_RMSE = sqrt(mean(gam_fit_on_weighted_gam_preds - Weighted_Overlap_Prob)^2),

                Unweighted_gam_Preds_RMSE = sqrt(mean(gam_fit_on_unweighted_gam_preds - Overlap_Prob)^2),
                
                ) %>%
  dplyr::select(contains("RMSE")) %>%
  dplyr::distinct() %>%
  pivot_longer(cols = c(colnames(.)))
```

    ## # A tibble: 4 × 2
    ##   name                            value
    ##   <chr>                           <dbl>
    ## 1 Weighted_Overlap_RMSE     0.000178   
    ## 2 Overlap_RMSE              0.00163    
    ## 3 Weighted_gam_Preds_RMSE   0.000000426
    ## 4 Unweighted_gam_Preds_RMSE 0.00000520

All fits are comparable

*We’ll use the* `regressed weighted GAM probabilities` in the
`binomial exact test`

## Binomial Exact Test 2

Re-compute the binomial exact test using the
`regressed weighted GAM fit values` as the the
`Weighted Overlap Probabilities`

i.e. Regress the `Weighted Overlap Prob`s and use those regressed values
as the theoretical probs in a `binomial exact test`

``` r
GO_INFO_by_TERM_df = GO_INFO_by_TERM_df %>%
  dplyr::mutate(
    E_Weighted_Overlap_2 = NA_real_,
    E_Weighted_Overlap_Prob_2 = NA_real_,
    Weighted_Overlap_Prob_p_2 = NA_real_,
    Weighted_Overlap_Prob_ADJP_2 = NA_real_,
    E_Overlap_2 = NA_real_,
    E_Overlap_Prob_2 = NA_real_,
    Overlap_Prob_p_2 = NA_real_,
    Overlap_Prob_ADJP_2 = NA_real_
  )

for( i in 1:nrow(GO_INFO_by_TERM_df) ) {
  
  ## what's the expected weighted overlap? probability? -> 1000 draws; USE THE GAM PROB PREDICTS
    ## take the average of the 1000 draws for the overlap
  GO_INFO_by_TERM_df$E_Weighted_Overlap_2[i] = round(mean(
    rbinom(1000,
           size = GO_INFO_by_TERM_df$Num_Genes_in_Term[i],
           prob = GO_INFO_by_TERM_df$gam_fit_on_weighted_gam_preds[i]
           )
    )
  )
  
  GO_INFO_by_TERM_df$E_Overlap_2[i] = round(mean(
    rbinom(1000,
           size = GO_INFO_by_TERM_df$Num_Genes_in_Term[i],
           prob = GO_INFO_by_TERM_df$gam_fit_on_unweighted_gam_preds[i]
           )
  )
  )
}
  
    ## divide predicted overlap counts by total number of genes
GO_INFO_by_TERM_df$E_Weighted_Overlap_Prob_2 = 
  
  GO_INFO_by_TERM_df$E_Weighted_Overlap_2 / max(GO_INFO_by_TERM_df$Num_Genes_in_Term)


GO_INFO_by_TERM_df$E_Overlap_Prob_2 = 
  
  GO_INFO_by_TERM_df$E_Overlap_2 / max(GO_INFO_by_TERM_df$Num_Genes_in_Term)

    ## perform an exact binomial test and extract p-values
for( i in 1:nrow(GO_INFO_by_TERM_df) ) {
  
  GO_INFO_by_TERM_df$Weighted_Overlap_Prob_p_2[i] = 
    
    binom.test(GO_INFO_by_TERM_df$Num_Genes_in_Term_from_List[i],
               GO_INFO_by_TERM_df$Num_Genes_in_Term[i],
               (GO_INFO_by_TERM_df$E_Weighted_Overlap_Prob_2[i]),
               alternative = "two.sided")$p.value
  
  GO_INFO_by_TERM_df$Overlap_Prob_p_2[i] = 
    
    binom.test(GO_INFO_by_TERM_df$Num_Genes_in_Term_from_List[i],
               GO_INFO_by_TERM_df$Num_Genes_in_Term[i],
               (GO_INFO_by_TERM_df$E_Overlap_Prob_2[i]),
               alternative = "two.sided")$p.value
  
}


GO_INFO_by_TERM_df$Weighted_Overlap_Prob_ADJP_2 = 
  
  (GO_INFO_by_TERM_df$Weighted_Overlap_Prob_p_2 / nrow(GO_INFO_by_TERM_df) ) * 0.05


GO_INFO_by_TERM_df$Overlap_Prob_ADJP_2 = 
  
  (GO_INFO_by_TERM_df$Overlap_Prob_p_2 / nrow(GO_INFO_by_TERM_df) ) * 0.05
```

## Do Exact Test Results Confirm?

![](GO-Analysis_files/figure-gfm/GO%20Bar%20Plot-1.png)<!-- -->

    ##                                                                                GO_Term
    ## 1                                                   actin filament bundle distribution
    ## 2                                        heparan sulfate 6-O-sulfotransferase activity
    ## 3           regulation of CDP-diacylglycerol-serine O-phosphatidyltransferase activity
    ## 4  positive regulation of CDP-diacylglycerol-serine O-phosphatidyltransferase activity
    ## 5                        positive regulation of serine C-palmitoyltransferase activity
    ## 6                          intrinsic component of neuronal dense core vesicle membrane
    ## 7                                   intrinsic component of dense core granule membrane
    ## 8                           integral component of neuronal dense core vesicle membrane
    ## 9                                               alpha-1,2-mannosyltransferase activity
    ## 10                                                           GTPase inhibitor activity
    ## 11                                             oncostatin-M-mediated signaling pathway
    ## 12                                                     mu-type opioid receptor binding
    ## 13                                                                   kininogen binding
    ## 14                                                                             C-fiber
    ## 15                                                 modulation by virus of host process
    ##    Num_Genes_in_Term_from_List Num_Genes_in_Term Weighted_Overlap_Prob
    ## 1                            2                 2          1.144145e-09
    ## 2                            2                 3          1.716217e-09
    ## 3                            2                 3          1.716217e-09
    ## 4                            2                 3          1.716217e-09
    ## 5                            2                 3          1.716217e-09
    ## 6                            2                 3          1.716217e-09
    ## 7                            2                 3          1.716217e-09
    ## 8                            2                 3          1.716217e-09
    ## 9                            2                 4          2.288290e-09
    ## 10                           2                 4          2.288290e-09
    ## 11                           2                 4          2.288290e-09
    ## 12                           2                 4          2.288290e-09
    ## 13                           2                 4          2.288290e-09
    ## 14                           2                 4          2.288290e-09
    ## 15                           2                 4          2.288290e-09
    ##    weighted_gam_preds gam_fit_on_weighted_gam_preds E_Weighted_Overlap_Prob_2
    ## 1        0.0001810468                  2.985022e-06                         0
    ## 2        0.0001810468                  2.985022e-06                         0
    ## 3        0.0001810468                  2.985022e-06                         0
    ## 4        0.0001810468                  2.985022e-06                         0
    ## 5        0.0001810468                  2.985022e-06                         0
    ## 6        0.0001810468                  2.985022e-06                         0
    ## 7        0.0001810468                  2.985022e-06                         0
    ## 8        0.0001810468                  2.985022e-06                         0
    ## 9        0.0001810468                  2.985022e-06                         0
    ## 10       0.0001810468                  2.985022e-06                         0
    ## 11       0.0001810468                  2.985022e-06                         0
    ## 12       0.0001810468                  2.985022e-06                         0
    ## 13       0.0001810468                  2.985022e-06                         0
    ## 14       0.0001810468                  2.985022e-06                         0
    ## 15       0.0001810468                  2.985022e-06                         0
    ##    Weighted_Overlap_Prob_ADJP_2
    ## 1                             0
    ## 2                             0
    ## 3                             0
    ## 4                             0
    ## 5                             0
    ## 6                             0
    ## 7                             0
    ## 8                             0
    ## 9                             0
    ## 10                            0
    ## 11                            0
    ## 12                            0
    ## 13                            0
    ## 14                            0
    ## 15                            0

# Conclusions

**No, test results don’t confirm because values are too small to
determine if probabilities or counts are significantly different than by
chance.**

But when removing highly specific `GO Terms` (with only 1 gene) and
arranging by the least probable weighted overlap, the results are
incredibly plausible:

We confirmed `GO Term`s we suspected to be intimately and robustly
releated to the *Tmem184b*, specifically, `C-fibers`, and
`oncostatin-M-mediated signialing pathway`
