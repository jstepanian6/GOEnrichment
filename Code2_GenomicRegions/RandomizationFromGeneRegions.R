library(dplyr)
library(tidyr)
library(purrr)
library(rtracklayer)
library(biomaRt)
library(GenomicRanges)

setwd("/home/jstepanian/Documents/Huber_Lab/FromGeneRegions/")

# Read regions file (expecting columns: chr, start, end)
#regions <- read.table("high_high.bed", header = FALSE, col.names = c("chr", "start", "end"))

regions <- read.table("high_low.bed", header = FALSE, col.names = c("chr", "start", "end"))


# Standardize chromosome naming and coordinate system
regions_df <- regions %>%
  mutate(
    chr = ifelse(grepl("^chr", chr), chr, paste0("chr", chr)), 
    start = start + 1                                            
  )

n_perm <- 1000
seed <- 123
j <- 3 # Minimum number of candidate genes per GO term

options(timeout = 1200)

GO_region_permutation <- function(regions_df, n_perm = 1000, seed = 123, j = 3) {
  set.seed(seed)
  
  if(!all(c("chr", "start", "end") %in% names(regions_df))) {
    stop("Input must have columns: chr, start, end")
  }
  
  regions_df <- regions_df %>%
    mutate(
      region_size = end - start + 1,
      region_id = paste0("region_", row_number())
    )
  
  message("Input: ", nrow(regions_df), " genomic regions")
  message("Region size range: ", min(regions_df$region_size), " - ", max(regions_df$region_size), " bp")
  message("Median region size: ", median(regions_df$region_size), " bp")
  
  # Ensembl gene annotations
  message("Connecting to Ensembl...")
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  message("Fetching all genes with positions and GO terms...")
  all_genes <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", 
                   "start_position", "end_position", "go_id"),
    mart = mart
  )
  
  message("Raw data: ", nrow(all_genes), " gene-GO pairs")
  
  all_genes <- all_genes %>%
    filter(
      chromosome_name %in% c(1:22, "X", "Y"),
      !is.na(go_id),
      go_id != ""
    ) %>%
    mutate(
      chr = paste0("chr", chromosome_name),
      gene_length = abs(end_position - start_position) + 1
    )
  
  message("After filtering: ", nrow(all_genes), " gene-GO pairs on standard chromosomes")
  
  gene_summary <- all_genes %>%
    group_by(ensembl_gene_id) %>%
    summarise(
      hgnc_symbol = dplyr::first(na.omit(hgnc_symbol)),
      chr = dplyr::first(chr),
      start_position = dplyr::first(start_position),
      end_position = dplyr::first(end_position),
      gene_length = dplyr::first(gene_length),
      go_terms = list(unique(go_id)),
      .groups = 'drop'
    ) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  message("Unique genes: ", nrow(gene_summary))
  
  # Function: find overlapping genes
  find_overlapping_genes <- function(chr, start, end, gene_data) {
    gene_data %>%
      filter(
        chr == !!chr,
        end_position >= start,
        start_position <= end
      )
  }
  
  # Check first few regions for overlaps
  message("\nChecking first 5 regions for potential overlaps:")
  for (i in 1:min(5, nrow(regions_df))) {
    region <- regions_df[i, ]
    overlaps <- find_overlapping_genes(region$chr, region$start, region$end, gene_summary)
    message(region$chr, ":", region$start, "-", region$end, " → ", nrow(overlaps), " overlaps")
  }
  message("")
  
  # Find overlapping genes for all regions
  message("Finding genes overlapping input regions...")
  candidate_genes_list <- lapply(1:nrow(regions_df), function(i) {
    region <- regions_df[i, ]
    overlapping <- find_overlapping_genes(region$chr, region$start, region$end, gene_summary)
    if(nrow(overlapping) > 0) {
      overlapping$region_id <- region$region_id
      overlapping$region_size <- region$region_size
    }
    overlapping
  })
  
  candidate_genes <- bind_rows(candidate_genes_list) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  message("Candidate genes found: ", nrow(candidate_genes))
  message("  - From ", length(unique(candidate_genes$region_id)), " regions with overlapping genes")
  
  if(nrow(candidate_genes) == 0) {
    stop("No genes found overlapping input regions! Check chromosome naming or coordinate system.")
  }
  
  # GO term extraction and filtering
  message("Extracting GO terms from candidate genes...")
  all_cand_go_terms <- unlist(candidate_genes$go_terms)
  message("Total GO term annotations in candidates: ", length(all_cand_go_terms))
  
  go_terms_counts <- data.frame(go_id = all_cand_go_terms) %>%
    group_by(go_id) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n >= j) %>%
    arrange(desc(n))
  
  go_terms <- go_terms_counts$go_id
  message(length(go_terms), " GO terms to test (with at least ", j, " genes)")
  
  # ========== DIAGNOSTIC INFORMATION ==========
  message("\n========== DIAGNOSTIC INFORMATION ==========")
  message("Total GO terms to test: ", length(go_terms))
  message("Number of permutations per GO term: ", n_perm)
  message("Total permutations to perform: ", length(go_terms) * n_perm)
  message("Estimated time (assuming 0.01 sec per permutation): ~", 
          round(length(go_terms) * n_perm * 0.01 / 60, 1), " minutes")
  message("===========================================\n")
  
  if(length(go_terms) > 0) {
    message("Top 5 GO terms by count:")
    print(head(go_terms_counts, 5))
  }
  
  # Prepare chromosome-specific gene lists
  genes_by_chr <- split(gene_summary, gene_summary$chr)
  message("Genes available per chromosome:")
  print(sapply(genes_by_chr, nrow))
  
  # Function: sample random region
  sample_random_region <- function(chr, region_size, gene_data_by_chr) {
    chr_genes <- gene_data_by_chr[[chr]]
    if(is.null(chr_genes) || nrow(chr_genes) == 0) return(data.frame())
    
    chr_start <- min(chr_genes$start_position)
    chr_end <- max(chr_genes$end_position)
    chr_length <- chr_end - chr_start
    
    if(region_size > chr_length) region_size <- chr_length
    
    max_start <- chr_end - region_size
    if(max_start <= chr_start) {
      rand_start <- chr_start
    } else {
      rand_start <- sample(chr_start:max_start, 1)
    }
    rand_end <- rand_start + region_size - 1
    
    chr_genes %>%
      filter(
        end_position >= rand_start,
        start_position <= rand_end
      )
  }
  
  # Function: test one GO term (with error handling)
  test_one_GO <- function(go_term, go_index = NULL) {
    tryCatch({
      cand_with_term <- candidate_genes[sapply(candidate_genes$go_terms, 
                                               function(terms) go_term %in% terms), ]
      obs <- nrow(cand_with_term)
      if (obs < j) return(NULL)
      
      perm_counts <- replicate(n_perm, {
        permuted_genes_list <- lapply(1:nrow(regions_df), function(i) {
          region <- regions_df[i, ]
          sample_random_region(region$chr, region$region_size, genes_by_chr)
        })
        
        permuted_genes <- bind_rows(permuted_genes_list) %>%
          distinct(ensembl_gene_id, .keep_all = TRUE)
        
        if(nrow(permuted_genes) == 0) return(0)
        
        has_go_term <- sapply(permuted_genes$go_terms, function(terms) go_term %in% terms)
        sum(has_go_term)
      })
      
      pval <- (sum(perm_counts >= obs) + 1) / (n_perm + 1)
      
      # Show details for first 3 GO terms
      if(!is.null(go_index) && go_index <= 3) {
        message("GO term ", go_term, " (", go_index, "/", length(go_terms), "): obs=", obs,
                ", perm range=[", min(perm_counts), "-", max(perm_counts), "]",
                ", mean=", round(mean(perm_counts), 2))
      }
      
      data.frame(
        GO = go_term,
        observed = obs,
        pval = pval,
        mean_perm = mean(perm_counts),
        sd_perm = sd(perm_counts),
        min_perm = min(perm_counts),
        max_perm = max(perm_counts)
      )
    }, error = function(e) {
      message("ERROR with GO term ", go_term, " (index ", go_index, "): ", e$message)
      return(NULL)
    })
  }
  
  # ========== PERMUTATION TESTS WITH PROGRESS TRACKING ==========
  message("\n========== Starting permutation tests ==========")
  start_time <- Sys.time()
  
  results_list <- vector("list", length(go_terms))
  
  for(i in seq_along(go_terms)) {
    # Progress update every 10 GO terms
    if(i %% 10 == 0 || i == 1) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      if(i > 1) {
        rate <- elapsed / i
        remaining <- rate * (length(go_terms) - i)
        message(sprintf("Progress: %d/%d GO terms (%.1f%%) | Elapsed: %.1f min | Est. remaining: %.1f min",
                        i, length(go_terms), 100*i/length(go_terms), elapsed, remaining))
      } else {
        message(sprintf("Progress: %d/%d GO terms (%.1f%%)",
                        i, length(go_terms), 100*i/length(go_terms)))
      }
    }
    
    results_list[[i]] <- test_one_GO(go_terms[i], go_index = i)
  }
  
  # Combine results
  results <- bind_rows(results_list)
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message("\n========== Permutation tests complete ==========")
  message("Total time: ", round(total_time, 2), " minutes")
  message("GO terms tested: ", nrow(results))
  
  results <- results %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    arrange(pval)
  
  return(list(
    results = results,
    candidate_genes = candidate_genes,
    regions_summary = regions_df
  ))
}

# ==========================================================
# Run the function
# ==========================================================
output <- GO_region_permutation(
  regions_df = regions_df,
  n_perm = n_perm,
  seed = seed,
  j = j
)

# ==========================================================
# View results
# ==========================================================
print("Top enriched GO terms:")
print(head(output$results, 10))

print("\nP-value summary:")
print(summary(output$results$pval))

print("\nNumber of candidate genes:")
print(nrow(output$candidate_genes))

write.csv(output$results, "results_permutationFromGeneRegion_HL.csv", row.names = FALSE)

# ==========================================================
# Plot histogram
# ==========================================================
all_cand_go_terms <- unlist(output$candidate_genes$go_terms)
go_terms_counts <- data.frame(go_id = all_cand_go_terms) %>%
  group_by(go_id) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(desc(n))

hist(go_terms_counts$n,
     breaks = 30,
     main = "Distribution of GO Term Counts in Candidate Genes",
     xlab = "Number of Candidate Genes per GO Term",
     ylab = "Frequency (Number of GO Terms)",
     col = "steelblue",
     border = "white")
abline(v = j, col = "red", lwd = 2, lty = 2)
text(x = j, y = max(hist(go_terms_counts$n, plot = FALSE)$counts) * 0.9,
     labels = paste("Threshold (j =", j, ")"), pos = 4, col = "red")