# GOEnrichment

Two complementary approaches for Gene Ontology (GO) enrichment analysis using randomization-based statistical testing.

## Overview

This repository contains R scripts for performing GO term enrichment analysis using randomization tests. 

## Methods Included

### 1. Gene List Randomization (`GoTerm_Randomization.R`)

Tests GO term enrichment in a list of genes by comparing to randomly sampled gene sets.

**Use case**: You have a list of candidate genes (e.g., from differential expression analysis) and want to know which GO terms are overrepresented.

**Input**: Plain text file with gene symbols (HGNC format), one per line

**Method**: For each GO term, randomly samples gene sets of the same size from the genome and compares the observed count to the randomization distribution.

### 2. Genomic Region Randomization (`RandomizationFromGeneRegions.R`)

Tests GO term enrichment in genes overlapping genomic regions by comparing to genes in random genomic locations.

**Use case**: You have genomic regions of interest (e.g., peaks from ChIP-seq, ATAC-seq, or Hi-C) and want to test if genes near/overlapping these regions are enriched for specific functions.

**Input**: BED file with genomic regions (chr, start, end)

**Method**: For each GO term, randomly places regions of the same size on the same chromosomes and compares genes overlapping random regions to your observed regions.


## Requirements

-   R version ≥ 4.0
-   Internet connection (for Ensembl queries)
-   Packages: dplyr, tidyr, purrr, rtracklayer, biomaRt, GenomicRanges

