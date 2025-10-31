# GO Term Enrichment Analysis via Genomic Regions

## Overview

This R script performs Gene Ontology (GO) term enrichment analysis for genomic regions using permutation-based statistical testing. It identifies GO terms that are significantly overrepresented in genes overlapping your input regions compared to random genomic locations.

## Method

The analysis:

1.  Identifies genes that overlap with your input genomic regions
2.  Extracts GO terms associated with these candidate genes
3.  For each GO term, performs permutation tests by:
    -   Randomly sampling genomic regions of the same size on the same chromosomes
    -   Counting how many genes with that GO term overlap the random regions
    -   Comparing observed counts to the permutation distribution
4.  Calculates p-values and adjusts for multiple testing (Benjamini-Hochberg)

## Requirements

### R Packages

```r
library(dplyr)
library(tidyr)
library(purrr)
library(rtracklayer)
library(biomaRt)
library(GenomicRanges)

```

Install missing packages:

```r
install.packages(c("dplyr", "tidyr", "purrr"))
BiocManager::install(c("rtracklayer", "biomaRt", "GenomicRanges"))

```

### Input File

A BED file with genomic regions (tab-separated, no header):

-   Column 1: chromosome (e.g., `chr1` or `1`)
-   Column 2: start position (0-based or 1-based)
-   Column 3: end position

Example (`high_high.bed`):

```
chr1    1000000    1001000
chr2    5000000    5002000
chr3    10000000   10001500

```
It was also executed for other regions described in Table S26 described by Choin et al., 2021 including: 

 - Regions of high frequency in both ancient and present-day (PD)
   samples (high_high.bed)
  - Regions of high frequency in ancient and low in present day (PD)
   samples (high_low.bed)
   - Regions of low frequency in ancient and high in present day (PD)
   samples (low_high.bed)

## Configuration

Edit these parameters in the script:

```r
# Input file
regions <- read.table("high_high.bed", ...)

# Analysis parameters
n_perm <- 1000  # Number of permutations (more = more accurate but slower)
seed <- 123     # Random seed for reproducibility
j <- 3          # Minimum genes per GO term to test

```

## Usage

1.  Place your BED file in the working directory
2.  Update the file path in the script:
    
    ```r
    setwd("/your/working/directory/")
    
    ```
    
3.  Run the script:
    
    ```r
    source("PermutationFromGeneRegions.R")
    
    ```
    It can also be executed using RStudio or the command 
	```
    Rscript RandomizationFromGeneRegions.R
	   ```
    
## Output

### Console Output

-   Progress updates during permutation testing
-   Summary statistics for input regions
-   Number of candidate genes and GO terms tested
-   Estimated runtime

### Files Generated

**`results_permutationFromGeneRegion_HH.csv`** - Main results table with columns:

-   `GO`: GO term identifier
-   `observed`: Number of candidate genes with this GO term
-   `pval`: Raw permutation p-value
-   `mean_perm`: Mean count in permutation distribution
-   `sd_perm`: Standard deviation of permutation distribution
-   `min_perm`: Minimum count in permutations
-   `max_perm`: Maximum count in permutations
-   `padj`: Adjusted p-value (Benjamini-Hochberg correction)

*NOTE:* Remember to change the name of the output file if the input was changed. 

### Visualization

The script generates a histogram showing the distribution of GO term counts in candidate genes, with a red line indicating the minimum threshold.

## Interpretation

-   **Low p-values** (< 0.05): GO terms significantly enriched in your regions
-   **padj < 0.05**: Significant after multiple testing correction (recommended)
-   **observed > mean_perm**: GO term appears more than expected by chance

## Performance Notes

-   Runtime depends on: number of regions, GO terms tested, and permutations
-   Typical runtime: It took 360.55 minutes for example file, but it will range depending on the regions size
-   The script displays progress updates and time estimates
-   Increase `options(timeout = 1200)` if Ensembl download times out

## Additional notes

Although some parts of the code indicate that it performs permutations, it actually performs randomizations. THIS CODE DOES NOT PERFORM PERMUTATIONS. 

## Issues 

 - The code described takes into account all the go terms including those belonging to the genes that only have 1pb of overlap with the region. 
 - Execution time can be long. 
 - This code accounts for genomic region length but they are overlapped over coding regions, so the ideal scenario would be only include genomic regions of coding regions. 

## References 
Choin, J., Mendoza-Revilla, J., Arauna, L.R. _et al._ Genomic insights into population history and biological adaptation in Oceania. _Nature_  **592**, 583–589 (2021). https://doi.org/10.1038/s41586-021-03236-5. 

