
# GO Analysis Functions

# Find_Enriched_GOs_from_GO_INFO_df_fn ----

## Find GOs "Enriched" from a supplied list of genes/proteins (output of GO_INFO_fn)
  ## supply "1:length(list_of_interest)" into the function
Find_Enriched_GOs_from_GO_INFO_df_fn = function(Gene_Indeces){
## Create a variable by which to index the GO terms returned that are annotated to multiple genes in a list
Term_Idx = which(duplicated(c(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
)[which(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
) %in% Unique_GOs)][which(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
)[which(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
) %in% Unique_GOs)] != "")])) == TRUE)

## Find the terms
Enriched_GOs = unique(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
)[which(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
) %in% Unique_GOs)][which(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
)[which(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
) %in% Unique_GOs)] != "")][Term_Idx])


## Find how many times each GO Term is annotated from the list
## Put it in a table
Enriched_GO_counts = table(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
)[which(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
) %in% Unique_GOs)][which(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
)[which(str_split(
  as.vector(GENE_GO_INFO_df[ Gene_Indeces, 4]), pattern = ";", simplify = TRUE
) %in% Unique_GOs)] != "")][Term_Idx])

  ## Transform the table to a dataframe
Enriched_GO_counts_df = Enriched_GO_counts %>%
  as.data.frame() %>%
  arrange(factor(Var1, levels = Enriched_GOs)) %>%
  dplyr::mutate(GO_Term = Var1,
                Num_Genes_in_List = Freq) %>%
  dplyr::select(-Var1, -Freq)

  ## Export to the Global Envir.; return in the console
Enriched_GO_counts_df <<- Enriched_GO_counts_df

return(Enriched_GO_counts_df)
}


# Find_Genes_With_X_GO_Term_fn ----

## Find the genes/proteins that contain specific GO Terms (from GO_INFO_df)
  ## Find which genes/proteins from the list of interest are associated with a specific GO Term
Find_Genes_With_X_GO_Term_fn = function(GO_Term){
  return(GENE_GO_INFO_df$GeneID[which(
    apply(GENE_GO_INFO_df, 1, function(x) any(grepl(GO_Term, x)))
  )])
}

# Find_GOs_of_Gene_X_fn ----

## Find the GO Terms of a gene/protein of interest (from the list_of_interest; indexing the GO_INFO_df)
  ## Find the GO Terms of a particular gene/protein within the list of interest used to generate the GO_INFO_df
Find_GOs_of_Gene_X_fn = function(GENE_ID) {
  return(str_split(
    as.vector(GENE_GO_INFO[ which(GENE_GO_INFO_df$GeneID == GENE_ID) ,   3]), pattern = ";", simplify = TRUE)[
      which(str_split(
        as.vector(GENE_GO_INFO[ which(GENE_GO_INFO_df$GeneID == GENE_ID) ,    3]), pattern = ";", simplify = TRUE) %in% Unique_GOs)
    ]
  )
}

# Find_GOs_of_Two_Genes_fn ----

## Find GO terms of two genes/proteins of interest (from the list_of_interest; indexing the GO_INFO_df)
  ## What are the GO Terms of two genes/proteins of interest? Are they related?
Find_GOs_of_Two_Genes_fn = function(GENE_ID1, GENE_ID2) {
  TermIdx = c(which(  
    str_split(
      as.vector(GENE_GO_INFO_df[ which(GENE_GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
    )[which(str_split(
      as.vector(GENE_GO_INFO_df[ which(GENE_GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
    ) %in% Unique_GOs)] %in%
      
      str_split(
        as.vector(GENE_GO_INFO_df[ which(GENE_GO_INFO_df$GeneID == GENE_ID2), 4]), pattern = ";", simplify = TRUE
      )[which(str_split(
        as.vector(GENE_GO_INFO_df[ which(GENE_GO_INFO_df$GeneID == GENE_ID2), 4]), pattern = ";", simplify = TRUE
      ) %in% Unique_GOs)]
  ))
  ## Find which terms those are
  str_split(
    as.vector(GENE_GO_INFO_df[ which(GENE_GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
  )[which(str_split(
    as.vector(GENE_GO_INFO_df[ which(GENE_GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
  ) %in% Unique_GOs)][TermIdx]
  return(str_split(
    as.vector(GENE_GO_INFO_df[ which(GENE_GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
  )[which(str_split(
    as.vector(GENE_GO_INFO_df[ which(GENE_GO_INFO_df$GeneID == GENE_ID1), 4]), pattern = ";", simplify = TRUE
  ) %in% Unique_GOs)][TermIdx])
}


# Find_Genes_Related_By_GO_Term_fn  ----

## Find which GO Terms are shared by all genes/proteins (from the list_of_interest; indexing the GO_INFO_df)
Find_Genes_Related_By_GO_Term_fn = function(GO_Term1, GO_Term2){
  MatchIdx = c(which(GENE_GO_INFO_df$GeneID[which(
    apply(GENE_GO_INFO_df, 1, function(x) any(grepl(GO_Term1, x)))
  )] %in% GENE_GO_INFO_df$GeneID[which(
    apply(GENE_GO_INFO_df, 1, function(x) any(grepl(GO_Term2, x)))
  )]))
  ## Find which genes/proteins those are
  GENE_GO_INFO_df$GeneID[which(
    apply(GENE_GO_INFO_df, 1, function(x) any(grepl(GO_Term1, x)))
  )][MatchIdx]
  return(GENE_GO_INFO_df$GeneID[which(
    apply(GENE_GO_INFO_df, 1, function(x) any(grepl(GO_Term1, x)))
  )][MatchIdx])
}

# Query_GO_fn ----

 ## Find all genes associated with a regexpr gene/protein term
Query_GO_fn = function(Species_db, GO_db, string_terms){
  
  ## first check that there are GO Terms associated with a provided string
  keys_var = c(keys(GO_db, keytype = "TERM")[
    which(
      grepl(
        keys(GO_db, keytype = "TERM"),
        pattern = string_terms
      ) == TRUE )
  ])
  
  
  
  if (length(keys_var) > 0) { ## if there are GO Terms associated with a provided string, then
    
    ## Concatenate the GO.db terms containing a provided string; will serve as keys; character vector
    GO.Terms = c(keys(GO_db, keytype = "TERM")[
      which(
        grepl(
          keys(GO_db, keytype = "TERM"),
          pattern = string_terms) == TRUE)
    ])
    ## Concatenate the GO.db IDs of the GO.db Terms containing a provided string; character vector
    GO.IDs = c(keys(GO_db, keytype = "GOID")[
      which(
        grepl(
          keys(GO_db, keytype = "TERM"),
          pattern = string_terms) == TRUE)
    ])
    ## Concatenate the organism database terms mapped to the GO.db Terms containing a provided string; list of lists
    GO.list = c(mapIds(Species_db,
                       keys = c(as.character(
                         mapIds(GO_db,
                                keys = c(
                                  keys(GO_db, keytype = "TERM")[
                                    which(
                                      grepl(
                                        keys(GO_db, keytype = "TERM"),
                                        pattern = string_terms)
                                      == TRUE)
                                  ]
                                ),
                                keytype = "TERM",
                                column = "GOID",
                                multiVals = "list")
                       )),
                       keytype = "GOALL",
                       column = "SYMBOL",
                       multiVals = "list"))
    ## Find NAs
    rm_idx = c(which(as.character(is.na(GO.list[][])) == "TRUE"))
    ## Remove NAs
    if (length(rm_idx) >0) {
      GO.list = GO.list[-rm_idx]
      GO.Terms = GO.Terms[-rm_idx]
      GO.IDs = GO.IDs[-rm_idx]
    }
    ## Initialize the matrix that will eventually become a dataframe
    df = matrix(nrow = length(c(GO.Terms)), ncol = 4)
    colnames(df) = c("GO Term ID", "GO Term", "# of Genes", "Genes")
    
    ## Initialize the columns that will contain strings
    df[ , c(1,2,4)] = ""
    ## Initialize the column that will contain numbers/
    df[ , 3] = 0
    ## Fill the columns with GO IDs, GO Terms, the number of genes/proteins within them, and all the genes/proteins within them
    for (i in 1:length(GO.Terms)){
      df[i,1] = GO.IDs[i]
      df[i,2] = GO.Terms[i]
      x = c(as.character(unlist(as.vector(GO.list[i][]))))
      df[i,3] = paste(as.numeric(length(unique(x))))
      df[i,4] =  paste(unique(x), collapse = ";")
    }
    ## Make the matrix a dataframe
    df = as.data.frame(df)
    
    ## Concatenate all the genes/proteins involved in a desired set of GO terms into one list
    all_unique_genes =     
      unique(
        as.character(
          unlist(
            mapIds(Species_db,
                   keys = c(as.character(
                     mapIds(GO_db,
                            keys = c(
                              keys(GO_db, keytype = "TERM")[
                                which(
                                  grepl(
                                    keys(GO_db, keytype = "TERM"),
                                    pattern = string_terms)
                                  == TRUE)
                              ]
                            ),
                            keytype = "TERM",
                            column = "GOID",
                            multiVals = "list")
                   )),
                   keytype = "GOALL",
                   column = "SYMBOL",
                   multiVals = "list")
          )
        )
      )
    
    ## Remove NA terms
    all_unique_genes = all_unique_genes[!is.na(all_unique_genes)]
    
    ## Find aliases of all the unique genes/proteins
    aliases = unique(
      as.character(
        unlist(
          mapIds(Species_db,
                 keys = c(all_unique_genes),
                 keytype = "SYMBOL",
                 column = "ALIAS",
                 multiVals = "list")
          
        )
      )
    )
    
    print(head(all_unique_genes))
    print(c("GO Terms returned as 'GO.Terms' in Global Env.",
            "List of genes in all returned GO Terms returned as 'all_unique_genes' in Global Env.",
            "Dataframe housing all info returned as 'GO_df' in Global Env.",
            "List of list housing raw search results returned as 'GO.list' in Global Env.",
            "Aliases returned as 'aliases'."))
    
    GO.Terms <<- GO.Terms
    GO.IDs <<- GO.IDs
    all_unique_genes = all_unique_genes[!is.na(all_unique_genes)]
    all_unique_genes <<- all_unique_genes
    GO.list <<- GO.list
    GO_df <<- df
    aliases <<- aliases
    
  } else { ## if there are no GO Terms associated with a provided string term,
    
    ## return the same output format, but with "none"
    GO.Terms = "none"
    GO.IDs = "none"
    all_unique_genes = "none"
    GO.list = "none"
    aliases = "none"
    
    if (grepl(string_terms, pattern = "\\|") == TRUE) {
      
      string_terms = c(str_split(as.vector(string_terms), pattern = "\\|", simplify = TRUE))
      
      df = matrix(nrow = length(string_terms), ncol = 4)
      colnames(df) = c("GO Term ID", "GO Term", "# of Genes", "Genes")
      
      for (i in 1:length(string_terms)) {
        
        df[i,2] = string_terms[i]
        
      }
      df[,3] = 0
      df[,c(1,4)] = ""
      df = as.data.frame(df)
       
    } else {
        
      df = matrix(nrow = 1, ncol = 4)
      colnames(df) = c("GO Term ID", "GO Term", "# of Genes", "Genes")
      df = as.data.frame(df)
      df$`GO Term` = string_terms
      df$`GO Term ID` = ""
      df$`# of Genes` = 0
      df$Genes = ""
      
    }
    
    print(c("GO Terms returned as 'GO.Terms' in Global Env.",
            "List of genes in all returned GO Terms returned as 'all_unique_genes' in Global Env.",
            "Dataframe housing all info returned as 'GO_df' in Global Env.",
            "List of list housing raw search results returned as 'GO.list' in Global Env.",
            "Aliases returned as 'aliases'."))
    
    
    GO.Terms <<- GO.Terms
    GO.IDs <<- GO.IDs
    all_unique_genes <<- all_unique_genes
    GO.list <<- GO.list
    GO_df <<- df
    aliases <<- aliases
    
  }
}

# GO_INFO_fn ----

## Create the GO_INFO_df 
## Extract GO information from all genes/proteins in a given (Homo Sapiens or mus musculus)

GO_INFO_fn = function(list_of_interest, species) {
  
  if (species == "human"){
    
    ## Find human Entrez IDs (Ensembl IDs) from the input list
    Ensembl_IDs = as.integer(
      mapIds(org.Hs.eg.db,
             keys = as.character(
               c(list_of_interest[])
             ),
             keytype = "SYMBOL",
             column = "ENTREZID")
    )
    ## Extract GO Terms for all human genes/proteins in the input list
    GO_IDs = c(
      mapIds(org.Hs.eg.db, 
             keys = as.character(
               c(list_of_interest[])
             ),
             keytype = "SYMBOL",
             column = "GOALL",
             multiVals = "list")
    )
    ## Find the human genes/proteins in the input list that don't have Ensembl/Entrez IDs
    indeces_without_Ensembl_IDs = c(
      which(is.na(Ensembl_IDs[]) == TRUE)
    )
    ## Find the human genes/proteins that don't have GO IDs
    indeces_without_GO_IDs = as.integer(c(
      which(is.na(GO_IDs[]) == TRUE)
    ))
    
    ## Remove genes/proteins that don't have IDs (causes errors)
    if (
      
      length(indeces_without_Ensembl_IDs) >= 1 &
      length(indeces_without_GO_IDs) == 0
      
    ) {
      
      list_of_interest = list_of_interest[-indeces_without_Ensembl_IDs]
    
      } else if (
      
        length(indeces_without_Ensembl_IDs) == 0 &
        length(indeces_without_GO_IDs) >= 1
        
    ) {
      
      list_of_interest = list_of_interest[-c(indeces_without_GO_IDs)]
    
      } else if (
        
        length(indeces_without_Ensembl_IDs) >= 1 &
        length(indeces_without_GO_IDs) >= 1
      
    ) {
      
      list_of_interest = list_of_interest[-c(unique(indeces_without_GO_IDs, indeces_without_Ensembl_IDs))]
    
      } else {
      
        list_of_interest
    }
    
    ## Re-query the GO IDs having removed 
    GO_IDs = c(
      mapIds(org.Hs.eg.db, 
             keys = as.character(
               c(list_of_interest[])
             ),
             keytype = "SYMBOL",
             column = "GOALL",
             multiVals = "list")
    )
    
  } else if (species == "mouse") {
    
    ## Find mouse Entrez IDs (Ensembl IDs) from the input list
    Ensembl_IDs = as.integer(
      mapIds(org.Mm.eg.db,
             keys = as.character(
               c(list_of_interest[])
             ),
             keytype = "SYMBOL",
             column = "ENTREZID")
    )
    ## Extract GO BP Terms for all mouse genes/proteins in the input list
    GO_IDs = c(
      mapIds(org.Mm.eg.db, 
             keys = as.character(
               c(list_of_interest[])
             ),
             keytype = "SYMBOL",
             column = "GOALL",
             multiVals = "list")
    )
    ## Find the genes/proteins in the input list that don't have Ensembl/Entrez IDs
    indeces_without_Ensembl_IDs = c(
      which(is.na(Ensembl_IDs[]) == TRUE)
    )
    ## Find the human genes/proteins that don't have GO IDs
    indeces_without_GO_IDs = as.integer(c(
      which(is.na(GO_IDs[]) == TRUE)
    ))
    
    ## Remove genes/proteins that don't have IDs (causes errors)
    if (
      
      length(indeces_without_Ensembl_IDs) >= 1 &
      length(indeces_without_GO_IDs) == 0
      
    ) {
      
      list_of_interest = list_of_interest[-indeces_without_Ensembl_IDs]
      
      } else if (
        
        length(indeces_without_Ensembl_IDs) == 0 &
        length(indeces_without_GO_IDs) >= 1
        
    ) {
        
        list_of_interest = list_of_interest[-c(indeces_without_GO_IDs)]
        
      } else if (
        
        length(indeces_without_Ensembl_IDs) >= 1 &
        length(indeces_without_GO_IDs) >= 1
      
      ) {
        
        list_of_interest = list_of_interest[-c(unique(indeces_without_GO_IDs, indeces_without_Ensembl_IDs))]
      
      } else {
        
        list_of_interest
      
      }
    
    
    ## Re-query the GO IDs having removed genes/proteins without one
    GO_IDs = c(
      mapIds(org.Mm.eg.db, 
             keys = as.character(
               c(list_of_interest[])
             ),
             keytype = "SYMBOL",
             column = "GOALL",
             multiVals = "list")
    )
    
  }
  
  
  
  # GO Term dataframe storage for the genes in the list ----
  ## Create a matrix to house GO information for each protein/gene
  GENE_GO_INFO = matrix(nrow = length(c(list_of_interest)), ncol = 4)
  # rownames(GENE_GO_INFO) = c(list_of_interest)
  colnames(GENE_GO_INFO) = c("GeneID", "Num_GO_Terms", "GO_Term_IDs", "GO_Terms")
  
  for (i in 1:length(list_of_interest)) {
    
    GENE_GO_INFO[i,1] = list_of_interest[i]
    
  }
  
  ## Fill the first column with zeros for genes without GO Terms
  GENE_GO_INFO[ , 2] = 0
  ## Fill the second and third columns with "None"s for genes without GO Terms
  GENE_GO_INFO[ , c(3,4)] = "None"
  ## Fill the second column with the concatenated GO Term IDs associated with each gene, skipping over the genes with 0
  for (i in 1:length(GO_IDs)){
    
    GENE_GO_INFO[i,3] = paste(names(mapIds(GO.db, GO_IDs[[i]][], "TERM", "GOID")), collapse = ";")
    
  }
  ## Extract the GO IDs and remove the redundant ones;
  ## Fill the first column with the numbers of GO Terms associated with each gene, skipping over the genes with 0
  for (i in 1:length(GO_IDs)){
    
    x = str_split(as.vector(GENE_GO_INFO[i,3]), pattern = ";", simplify = TRUE)
    x = c(as.character(x[which(x %in% x == TRUE)]))
    GENE_GO_INFO[i,3] = paste(unique(x), collapse = ";")
    GENE_GO_INFO[i,2] = paste(as.numeric(length(unique(x))))
    
  }
  ## Fill the third column with the concatenated GO Terms associated with each gene, skipping over the genes with 0
  for (i in 1:length(GO_IDs)){
    
    x = str_split(as.vector(GENE_GO_INFO[i,3]), pattern = ";", simplify = TRUE)
    GENE_GO_INFO[i,4] = paste(as.character(mapIds(GO.db, keys = x, keytype = "GOID", "TERM")), collapse = ";")
    
  }
  
  ## Extract all of the terms/IDs into one vector to find unique GO terms
  
  Unique_GOs = c(unique(as.character(str_split(as.vector(GENE_GO_INFO[,4]), pattern = ";", simplify = T))))
  Unique_GO_IDs = c(unique(as.character(str_split(as.vector(GENE_GO_INFO[,3]), pattern = ";", simplify = T))))
  
  if (any(Unique_GOs == "") ) {
    
    Unique_GOs = Unique_GOs[-which(Unique_GOs == "")]
    Unique_GO_IDs = Unique_GO_IDs[-which(Unique_GO_IDs == "")]
    
  }
  
  # GO df storage for all GOs ----
  ## Find all the genes/proteins within the list of terms
  
  GO.list = mapIds(org.Hs.eg.db,
                   keys = Unique_GO_IDs[],
                   keytype = "GOALL",
                   column = "SYMBOL",
                   multiVals = "list")
  
  ## Put the information into a matrix
  GO_INFO = matrix(nrow = length(Unique_GOs), ncol = 6)
  colnames(GO_INFO) = c("GO_Term_ID", "GO_Term", "Num_Genes", "Gene_IDs", "Num_Genes_from_List", "Gene_IDs_from_List")
  ## initialize columns
  GO_INFO[ , c(1,2,4,6)] = ""
  GO_INFO[ , c(3,5)] = 0
  
  ## loop through the terms
  for (i in 1:length(Unique_GOs) ){
    
    ## store the IDs
    GO_INFO[i,1] = Unique_GO_IDs[i]
    ## store the terms
    GO_INFO[i,2] = Unique_GOs[i]
    ## store the number of genes in each term
    GO_INFO[i,3] = as.numeric(length(GO.list[[i]]))
    ## store the genes in each term
    GO_INFO[i,4] = paste(c(GO.list[[i]]), collapse = ";")
    ## store the number of genes in the term also in the provided list
    GO_INFO[i,5] = as.numeric(length(which(c(list_of_interest) %in% GO.list[[i]] == T)))
    ## store the genes in the term also in the provided list
    GO_INFO[i,6] = paste(c(list_of_interest[which(c(list_of_interest) %in% GO.list[[i]] == T)]), collapse = ";")

  }
  
  ## convert to dataframe
  GO_INFO = as.data.frame(GO_INFO)
  
  ## Make the matrix a data frame, and add "GeneIDs" as the first column
  GENE_GO_INFO_df = as.data.frame(GENE_GO_INFO)
  
  # Export ----
  ## Export the GO INFO df, all unique GO Terms, and the filtered list of interest to the Global Envir.
  GENE_GO_INFO_df <<- GENE_GO_INFO_df
  GO_INFO_by_TERM_df <<- GO_INFO
  Unique_GOs <<- Unique_GOs
  Unique_GO_IDs <<- Unique_GO_IDs
  list_of_interest <<- list_of_interest
  
}
