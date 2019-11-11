
<!-- README.md is generated from README.Rmd. Please edit that file -->
fastaR
======

<!-- <!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) 
[![platform](https://img.shields.io/badge/R-%3E%20v3.5.1-brightgreen)](https://shields.io/category/platform-support) 
<!-- badges: end --> 

Fasta sequence manipulation is required while performing gene set or genome-wide analysis in RNASeq and ChIPSeq. Sequences of group of genes belonging to a pathway or biological term helps to determine the motifs/signature sequences associated with those genes. It also helps to predict regulators for differentially enriched genes. fastaR provides different ways to manipulate and gather information from the set of sequences. For instance, summary of input fasta file `(fa_summary)` or percentage of GC in given sequences `(fa_percent_GC)`. Additionally, gene list specific sequences `(fa_some_records)` or promoters for given genes can also be fetched `(get_promoter_from_feature)`. Together, fastaR provides easy way to analyze and manipulate sequences with only feature file `(bed or gff)`, reference genome file `(.fa or .fasta)` as input.

Install
-------

    devtools::install_github("sethiyap/fastaR")

### faUtils

-   `fa_some_records()`
-   `fa_size()`
-   `fa_summary()`
-   `fa_percent_GC()`

### getUtils

-   `get_fasta_from_bed()`
-   `get_promoter_from_feature()`
