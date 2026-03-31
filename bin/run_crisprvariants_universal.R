#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Arguments:
# 1: direction (working directory)
# 2: csv_file (sample table)
# 3: amplicon (reference sequence)
# 4: cutSite (target location)
# 5: blat.T.name (chromosome/sequence name)
# 6: blat.T.start (start position)
# 7: blat.T.end (end position)
# 8: restable (results output filename)
# 9: instable (insertions output filename)

if (length(args) < 9) {
  stop("Usage: Rscript script.R <direction> <csv_file> <amplicon> <cutSite> <blat.T.name> <blat.T.start> <blat.T.end> <restable> <instable>")
}

direction <- args[1]
csv_file <- args[2]
amplicon <- args[3]
cutSite <- as.numeric(args[4])
blat.T.name <- args[5]
blat.T.start <- as.numeric(args[6])
blat.T.end <- as.numeric(args[7])
restable <- args[8]
instable <- args[9]

library("CrispRVariants")
library("GenomicRanges")
library(rtracklayer)
library(reshape2)
library(stringr)

# Set working directory
setwd(direction)

# Read sample table
md <- read.table(csv_file, sep = ",", header = TRUE)

# Set up genomic range
blat.strand <- "+"
gd <- GRanges(blat.T.name, blat.strand, 
              ranges = IRanges(blat.T.start, blat.T.end))

# Chimera treatment
treat.chimeras <- "exclude"

# Find all BAM subdirectories (parameter combinations)
bam_base_dir <- file.path(direction, "bam")
param_dirs <- list.dirs(bam_base_dir, full.names = FALSE, recursive = FALSE)

# Filter to only directories that match parameter naming pattern (A*_B*_O*_E*)
param_dirs <- param_dirs[grepl("^A[0-9]+_B[0-9]+_O.*_E", param_dirs)]

if (length(param_dirs) == 0) {
  stop("No parameter combination directories found in ", bam_base_dir)
}

cat("Found", length(param_dirs), "parameter combination(s) to process:\n")
cat(paste(param_dirs, collapse = "\n"), "\n\n")

cat("Found", nrow(md), "sample(s) to process:\n")
cat(paste(md$Run, collapse = "\n"), "\n\n")

# Loop through each parameter combination directory
for (param_name in param_dirs) {
  cat("========================================\n")
  cat("Processing parameter combination:", param_name, "\n")
  cat("========================================\n")
  
  output_dir <- file.path(bam_base_dir, param_name)
  
  # Loop through each sample
  for (i in 1:nrow(md)) {
    sample_name <- md$Run[i]
    
    cat("  Processing sample:", sample_name, "\n")
    
    # Create results directory for this sample and parameter combination
    results_dir <- file.path(direction, "results", param_name, sample_name)
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Build BAM file path for this sample
    bam.file.path <- paste0(output_dir, "/", sample_name, "_mapped_clean.bam")
    
    # Check if BAM file exists
    if (!file.exists(bam.file.path)) {
      cat("    WARNING: BAM file not found:", bam.file.path, "\n")
      cat("    Skipping this sample.\n\n")
      next
    }
    
    cat("    BAM file found. Starting analysis...\n")
    
    ### Analyse ################################################################
    ############################################################################
    
    tryCatch({
      # Run CrispRVariants analysis for single sample
      crispr.set <- readsToTarget(bam.file.path, gd,
                                  reference = amplicon,
                                  target.loc = cutSite,
                                  names = sample_name,
                                  chimeras = treat.chimeras,
                                  collapse.pairs = FALSE,
                                  chimera.to.target = 200,
                                  minoverlap = nchar(amplicon)/2)
      
      # Calculate mutation efficiency
      mut_eff <- mutationEfficiency(crispr.set)
      cat("    Mutation efficiency calculated.\n")
      
      # Get variant counts table
      var_counts <- variantCounts(crispr.set)
      
      # Add indel type classification
      byType <- crispr.set$classifyVariantsByType()
      byType <- unlist(byType)
      var_counts <- cbind.data.frame(var_counts, byType)
      var_counts$byType[var_counts$byType == "SNV"] <- "no variant"
      
      cat("    Variant types:\n")
      print(table(var_counts$byType))
      
      # Output file paths
      results_file <- file.path(results_dir, restable)
      insertions_file <- file.path(results_dir, instable)
      plot_file <- file.path(results_dir, paste0(sample_name, "_", param_name, "_alignments.png"))
      
      # Write results
      write.table(var_counts, results_file, sep = "\t", quote = FALSE, row.names = TRUE)
      cat("    Results written to:", results_file, "\n")
      
      # Get and save insertion sequences
      insseq <- crispr.set$.getInsertions()
      if (nrow(insseq) > 0) {
        insseq$idxs <- vapply(insseq$idxs, paste, collapse = ",", character(1L))
        write.table(insseq, insertions_file, sep = "\t", quote = FALSE, row.names = FALSE)
        cat("    Insertions written to:", insertions_file, "\n")
      } else {
        cat("    No insertions found for this sample.\n")
        # Create empty file to indicate completion
        write.table(data.frame(), insertions_file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
      
      # Create alignment plot
      png(plot_file,
          width = 15, height = 12, bg = "white", units = "in", res = 300)
      plotVariants(crispr.set, col.wdth.ratio = c(3,1),
                   plotFreqHeatmap.args = list(top.n = 50,
                                               type = "proportions",
                                               min.freq = 0.05,
                                               plot.text.size = 2),
                   plotAlignments.args = list(max.insertion.size = 150,
                                              top.n = 50,
                                              legend.cols = 5,
                                              plot.text.size = 2,
                                              target.loc = cutSite,
                                              highlight.guide = FALSE,
                                              highlight.pam = FALSE,
                                              min.freq = 0.05))
      dev.off()
      cat("    Plot saved to:", plot_file, "\n")
      
      cat("    SUCCESS: Analysis complete for", sample_name, "\n\n")
      
    }, error = function(e) {
      cat("\n")
      cat("    ****************************************\n")
      cat("    ERROR in processing sample:", sample_name, "\n")
      cat("    ****************************************\n")
      cat("    Error message:\n")
      cat("    ", conditionMessage(e), "\n")
      cat("    ****************************************\n")
      cat("    Continuing with next sample...\n\n")
      
      # Write error log to file
      error_log <- file.path(results_dir, paste0(sample_name, "_", param_name, "_ERROR.log"))
      writeLines(c(
        paste("Error occurred at:", Sys.time()),
        paste("Parameter combination:", param_name),
        paste("Sample:", sample_name),
        "",
        "Error message:",
        conditionMessage(e),
        "",
        "Full error:",
        paste(capture.output(print(e)), collapse = "\n")
      ), error_log)
      cat("    Error details saved to:", error_log, "\n\n")
    })
  }
  
  cat("Completed all samples for", param_name, "\n\n")
}

cat("========================================\n")
cat("All parameter combinations and samples processed!\n")
cat("========================================\n")