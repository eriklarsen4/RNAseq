# Bioinformatics
Template R code for bioinformatics analysis pipelines.

RNAseq scripts (WGCNA, Hippocampus Bioinformatics, Bioinformatics) analyze gene expression changes in mice with Tmem184b gene trap. These use mRNA TPM and DESeq2 DGEA files.
  WGCNA (As of now, "Manual Net Construction" script) creates TOMplots of TOMs or adjacencies as end-products of WGCNA.
    Hippocampus Bioinformatics creates MA and volcano plots.
    Performs downstream functional analysis of desired modules (through Enrichr metadatabase).
      Creates barplots of these Enrichr results.
        In progress: exporting WGCNA data compatible with cytoscape for network visualization.

MS script analyzes protein interactions from HEKs transfected with mouse TMEM184B. This uses a raw MUDPIT CSV file.
  For AP/IP-MS, creates scatterplots of various data processing methods.
    Performs downstream functional analysis of select proteins as candidate TMEM184B binding partners.
        In progress: SAINT algorithm for determining "true positive" binding partner candidates (against "false positive"s);
        In progress: exporting interaction data compatible with cytoscape for protein-protein interaction network visualization.

Long-term goal: integration of RNAseq gene changes with protein-protein interaction network.
