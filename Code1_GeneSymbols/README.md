# GO Term Enrichment Analysis via Gene Symbols

## Overview

This R script performs Gene Ontology (GO) term enrichment analysis for a list of genes using randomization-based statistical testing. It identifies GO terms that are significantly overrepresented in your gene list compared to random gene sets of the same size. It uses gene symbols. 

## Method

The analysis:

1.  Takes your input list of candidate genes (HGNC symbols)
2.  Retrieves all genes and their GO term annotations from Ensembl
3.  For each GO term present in your candidate genes:
    -   Counts how many candidate genes have that GO term (observed count)
    -   Randomly samples gene sets of the same size from the genome
    -   Counts how many genes in each random sample have that GO term
    -   Compares the observed count to the randomization distribution
4.  Calculates p-values and adjusts for multiple testing (Benjamini-Hochberg)

**Key difference from permutation**: This is a randomization test where genes are sampled _with replacement_ from the full genome, not a permutation where genes are shuffled among fixed positions.
**THIS CODE DOES NOT PERFORM A PERMUTATION**

## Requirements

### R Packages

```r
library(dplyr)
library(tidyr)
library(purrr)
library(rtracklayer)
library(biomaRt)

```

Install missing packages:

```r
install.packages(c("dplyr", "tidyr", "purrr"))
BiocManager::install(c("rtracklayer", "biomaRt"))

```

### Input File

A plain text file with one gene symbol per line (HGNC format):

Example (`HH.txt`):

```
ENSG00000227169
LINC01646MFFP1
TMEM51
ENSG00000239670

```

**Important**: Gene symbols must be in HGNC format (standard human gene symbols).

## Configuration

Edit these parameters in the script:

```r
# Input file
gene_list <- read.table("HH.txt")

# Analysis parameters
n_perm <- 1000  # Number of randomizations (more = more accurate but slower)
seed <- 123     # Random seed for reproducibility
j <- 3          # Minimum genes per GO term to test

# Timeout for Ensembl connection (seconds)
options(timeout = 1200)  # CRITICAL, otherwise the background download would not work

```

## Usage

1.  Place your gene list file in the working directory
2.  Update the file path in the script:
    
    ```r
    setwd("/your/working/directory/")
    
    ```
    
3.  Run the script:
    
    ```r
    source("GoTerm_Randomization.R")
    
    ```
    

## Output

### Console Output

-   Number of background genes retrieved from Ensembl
-   Number of candidate genes found in the database
-   Number of GO terms being tested
-   Progress information for randomization tests
-   Debug information for the first few GO terms

### Files Generated

**`GoTerm_randomized_HH.csv`** - Main results table with columns:

-   `GO`: GO term identifier (e.g., GO:0006355)
-   `observed`: Number of candidate genes with this GO term
-   `pval`: Raw randomization p-value
-   `mean_perm`: Mean count across all randomizations
-   `max_perm`: Maximum count in randomizations
-   `padj`: Adjusted p-value (Benjamini-Hochberg correction)

### Result Interpretation

-   **Low p-values** (< 0.05): GO terms significantly enriched in your gene list
-   **padj < 0.05**: Significant after multiple testing correction (recommended threshold)
-   **observed > mean_perm**: GO term appears more frequently than expected by chance

## Example Output

```
GO             observed  pval    mean_perm  max_perm  padj
GO:0006355     45        0.001   12.3       18        0.015
GO:0045944     38        0.003   10.8       16        0.022
GO:0006357     32        0.012   15.4       22        0.048

```

## Performance Notes

-   **Typical runtime**: 5-30 minutes depending on:
    -   Number of candidate genes
    -   Number of GO terms tested
    -   Number of randomizations (`n_perm`)
   It took 2 min 15 sec for the HH.txt example. 
-   The script displays debug information for the first 3 GO terms tested
-   Progress messages show background gene counts and filtering steps



## Statistical Notes

### P-value Calculation

```
p-value = (# randomizations ≥ observed + 1) / (n_perm + 1)

```

The "+1" in numerator and denominator provides a conservative estimate and ensures p-values are never exactly 0.

### Multiple Testing Correction

The script uses Benjamini-Hochberg (BH) method to control False Discovery Rate (FDR). Use `padj < 0.05` for reporting significant results.

## Gene Length Consideration

The script calculates gene lengths but currently uses **uniform random sampling** (all genes have equal probability). This is appropriate when:

-   Gene length bias is not a concern for your analysis
-   Your candidate genes don't systematically differ in length from background

If length bias is a concern, the sampling could be modified to be length-weighted.

## Note 
Example files belong to supplementary material described in Table S26 by Choin et al., 2021
## References
Choin, J., Mendoza-Revilla, J., Arauna, L.R. _et al._ Genomic insights into population history and biological adaptation in Oceania. _Nature_  **592**, 583–589 (2021). [https://doi.org/10.1038/s41586-021-03236-5](https://doi.org/10.1038/s41586-021-03236-5).

