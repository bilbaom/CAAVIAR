#!/usr/bin/env Rscript
library(reshape2)
library(seqinr)

## Analyze insertions and indel statistics
# 0: Open crisprvariants result tables
# 1: Parse crisprvariants insertions table to asing each to a readid
# 2: BLAST insertions to AAV genome
# 3: Match AAV insertions with indels and calculate metrics related to indels (deletion insertion length, ratio...)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: step2_aav_insertions.R <results_csv> <insertions_csv> <sample_name>")
}

resname <- args[1]      # Full path to the CrispRVariants results.csv
insname <- args[2]      # Full path to the CrispRVariants insertions.csv
experiment <- args[3]   # The actual sample or target name

# NOTE: No setwd() is used! Nextflow manages the working directory automatically.

#### FUNCTIONS #################################################################

# Function to classify indels based on CIGAR labels
classify_indels <- function(cigar_labels) {
  classify_single_indel <- function(cigar) {
    operations <- unlist(strsplit(cigar, ","))
    operation_types <- gsub(".*:", "", operations)
    operation_types <- gsub("[0-9]+", "", operation_types)
    insertions <- sum(operation_types == "I")
    deletions <- sum(operation_types == "D")
    
    if (insertions > 1 && deletions == 0) return("multiple insertions")
    else if (insertions >= 1 && deletions >= 1) return("insertion/deletion")
    else if (insertions == 1 && deletions == 0) return("insertion")
    else return("other")
  }
  sapply(cigar_labels, classify_single_indel, USE.NAMES = FALSE)
}

# Function to expand dataframe by splitting comma-separated IDs
expand_dataframe_by_ids <- function(df) {
  expanded_rows <- list()
  for (i in 1:nrow(df)) {
    ids <- trimws(unlist(strsplit(df$idxs[i], ",")))
    n_ids <- length(ids)
    expanded_row <- df[rep(i, n_ids), ]
    expanded_row$idxs <- ids
    expanded_rows[[i]] <- expanded_row
  }
  result <- do.call(rbind, expanded_rows)
  rownames(result) <- NULL
  return(result)
}

# Function to compress df by id and keep all sequences
merge_by_sample_id <- function(df, id_col = "idxs", sample_col = "sample", seq_col = "seq", type_col = "type", cigar_col = "cigar_label") {
  df$sample_id <- paste(df[[sample_col]], df[[id_col]], sep = "_")
  result <- aggregate(
    df[[seq_col]], 
    by = list(sample_id = df$sample_id, sample = df[[sample_col]], idxs = df[[id_col]]), 
    FUN = function(x) paste(unique(x), collapse = ",")
  )
  names(result)[names(result) == "x"] <- seq_col
  
  type_values <- aggregate(df[[type_col]], by = list(sample_id = df$sample_id), FUN = function(x) x[1])
  result[[type_col]] <- type_values$x[match(result$sample_id, type_values$sample_id)]
  
  cigar_values <- aggregate(df[[cigar_col]], by = list(sample_id = df$sample_id), FUN = function(x) x[1])
  result[[cigar_col]] <- cigar_values$x[match(result$sample_id, cigar_values$sample_id)]
  
  return(result)
}

# Function to extract unique sequences and write FASTA file
extract_sequences_to_fasta <- function(df, seq_col = "seq", min_length = 14, output_file = "unique_insertions.fasta", id_prefix = "seqid") {
  all_sequences <- df[[seq_col]]
  split_sequences <- unlist(strsplit(all_sequences, ","))
  unique_sequences <- unique(trimws(split_sequences))
  long_sequences <- unique_sequences[nchar(unique_sequences) > min_length]
  long_sequences <- long_sequences[long_sequences != ""]
  
  if (length(long_sequences) == 0) {
    cat("No sequences longer than", min_length, "characters found.\n")
    return(data.frame(id = character(), sequence = character(), length = integer(), stringsAsFactors = FALSE))
  }
  
  seq_ids <- paste0(id_prefix, sprintf("%02d", seq_along(long_sequences)))
  sequences_df <- data.frame(id = seq_ids, sequence = long_sequences, length = nchar(long_sequences), stringsAsFactors = FALSE)
  
  write.fasta(sequences = as.list(long_sequences), names = seq_ids, file.out = output_file)
  return(sequences_df)
}

#### 0: OPEN VARIANT COUNTS AND INSERTIONS ####
varcounts <- read.table(resname, sep = "\t", header = T)
insseq <- read.table(insname, sep = "\t", header = T)

#### 1: PARSE INSERTIONS ####
last_col <- colnames(varcounts)[ncol(varcounts)]
if ("insertion" %in% varcounts[[last_col]]) {
  colnames(varcounts)[ncol(varcounts)] <- "byType.g1"
} else {
  stop("Error: no 'insertion' value found in the last column.")
}

samplenames <- colnames(varcounts)[-length(colnames(varcounts))]
names(samplenames) <- as.character(1:length(samplenames))
insseq$sample <- samplenames[as.character(insseq$sample)]
insseq$type <- classify_indels(insseq$cigar_label)

insseq_expand <- expand_dataframe_by_ids(insseq)
insseq_expand$count <- NULL
insseq_merged <- merge_by_sample_id(insseq_expand)

#### 2: BLAST INSERTIONS ####
sequences_info <- extract_sequences_to_fasta(
  df = insseq_expand,
  seq_col = "seq",
  min_length = 14,
  output_file = "unique_insertions.fasta",
  id_prefix = "seqid"
)

blast_db_path <- Sys.getenv("BLAST_DB_PATH")
if (blast_db_path == "") {
  blast_db_path <- "/home/bilbaom/mhaav25/ref/blast/blastdb"
}

has_sequences <- nrow(sequences_info) > 0

if (!has_sequences) {
  message("No insertion sequences to BLAST. Setting AAV insertions = 0.")
} else {
  # FIX: Removed the hardcoded 'direction' variable. BLAST writes directly to the current working directory.
  cmd <- paste("blastn -query unique_insertions.fasta ",
               "-db ", blast_db_path, " ",
               "-evalue 0.01 -word_size 12 -out blastResults.csv ",
               '-perc_identity 80 -task blastn -max_target_seqs 1000000 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"',
               sep = '')
  system(cmd)
}

#### 3: MAKE INDEL SUMMARY ####
# FIX: Hardcoded path removed
blast_file <- "blastResults.csv"

if (!has_sequences || !file.exists(blast_file) || file.info(blast_file)$size == 0) {
  message("No BLAST results. Setting AAV insertions = 0.")
  insseq_merged$aavins <- 0
  totalins <- data.frame(Var1 = samplenames, Var2 = TRUE, Freq = 0)
} else {
  blastResults <- read.table(blast_file, sep = "\t", header = FALSE)
  sequences_info_aav <- sequences_info[sequences_info$id %in% unique(blastResults$V1),]
  
  count_matches <- function(s) {
    parts <- unlist(strsplit(s, ","))
    sum(parts %in% sequences_info_aav$sequence)
  }
  
  insseq_merged$aavins <- sapply(insseq_merged$seq, count_matches)
  insseq_merged_aav <- insseq_merged[insseq_merged$aavins > 0,]
  
  totalins <- table(insseq_merged$sample, insseq_merged$aavins > 0)
  totalins <- as.data.frame(totalins)
  totalins <- totalins[totalins$Var2 == TRUE, ]
}

# Sanitize categories
varcounts$byType.g1 <- gsub("[ /]", "_", varcounts$byType.g1)
categories <- unique(varcounts$byType.g1)
total_reads <- colSums(varcounts[, samplenames, drop = FALSE])
cat_reads <- sapply(samplenames, function(s) tapply(varcounts[[s]], varcounts$byType.g1, sum))

complete_cat_reads <- matrix(0, nrow = length(categories), ncol = length(samplenames))
rownames(complete_cat_reads) <- categories
colnames(complete_cat_reads) <- samplenames
existing_categories <- rownames(cat_reads)
complete_cat_reads[existing_categories, ] <- cat_reads
cat_reads <- complete_cat_reads

edited_reads <- total_reads - cat_reads["no_variant", ]
cat_pct <- sweep(cat_reads, 2, total_reads, "/") * 100
edited_reads_pct <- (edited_reads / total_reads) * 100

summary_df <- data.frame(
  sample = samplenames,
  total_reads = total_reads,
  edited_reads = edited_reads,
  edited_reads_pct = edited_reads_pct,
  t(cat_reads),
  t(cat_pct),
  check.names = FALSE
)

ncat <- length(categories)
colnames(summary_df)[(ncol(summary_df)-ncat+1):ncol(summary_df)] <- paste0(categories, "_pct")

summary_df <- merge(summary_df, totalins, by.x = "sample", by.y = "Var1", all.x = TRUE)
summary_df$aavins <- ifelse(is.na(summary_df$Freq), 0, summary_df$Freq)
summary_df$Freq <- NULL
summary_df$Var2 <- NULL
summary_df$aavins_pct <- summary_df$aavins / summary_df$edited_reads * 100

#### 3.1: INS/DEL STATISTICS ####
parse_indel_string <- function(s) {
  events <- unlist(strsplit(s, ","))
  pos  <- as.integer(sub(":.*", "", events))
  size <- as.integer(sub("[ID]$", "", sub(".*:", "", events)))
  type <- sub(".*([ID])$", "\\1", events)
  data.frame(pos=pos, size=size, type=type, stringsAsFactors=FALSE)
}

indel_strings <- rownames(varcounts)
all_events <- list()

for (i in seq_along(indel_strings)) {
  s <- indel_strings[i]
  if (s == "no variant") next
  if (grepl("^SNV", s)) next
  
  events <- parse_indel_string(s)
  
  for (sample in samplenames) {
    n <- varcounts[i, sample]
    if (n > 0) {
      tmp <- events[rep(seq_len(nrow(events)), each=n), , drop=FALSE]
      tmp$sample <- sample
      all_events[[length(all_events)+1]] <- tmp
    }
  }
}

all_events <- do.call(rbind, all_events)

summary_indels <- do.call(rbind, lapply(split(all_events, all_events$sample), function(df) {
  ins <- df[df$type=="I", "size"]
  del <- df[df$type=="D", "size"]
  total_ins <- length(ins)
  total_del <- length(del)
  
  data.frame(
    sample = unique(df$sample),
    total_ins = total_ins,
    ins1 = sum(ins == 1),
    ins14 = sum(ins < 5),
    ins_gt4 = sum(ins > 4),
    ins_gt9 = sum(ins > 9),
    ins_gt14 = sum(ins > 14),
    mean_ins_size = if(total_ins > 0) mean(ins) else NA,
    median_ins_size = if(total_ins > 0) median(ins) else NA,
    total_del = total_del,
    del1 = sum(del == 1),
    del14 = sum(del < 5),
    del_gt4 = sum(del > 4),
    del_gt9 = sum(del > 9),
    del_gt14 = sum(del > 14),
    mean_del_size = if(total_del > 0) mean(del) else NA,
    median_del_size = if(total_del > 0) median(del) else NA,
    ins_gt4_pct  = if(total_ins > 0) sum(ins > 4) / total_ins * 100 else NA,
    ins_gt9_pct  = if(total_ins > 0) sum(ins > 9) / total_ins * 100 else NA,
    ins_gt14_pct = if(total_ins > 0) sum(ins > 14) / total_ins * 100 else NA,
    del_gt4_pct  = if(total_del > 0) sum(del > 4) / total_del * 100 else NA,
    del_gt9_pct  = if(total_del > 0) sum(del > 9) / total_del * 100 else NA,
    del_gt14_pct = if(total_del > 0) sum(del > 14) / total_del * 100 else NA,
    ins_del_ratio = if(total_del > 0) total_ins / total_del else NA,
    nhej = sum(ins == 1) + sum(del == 1),
    nhej_pct = if((total_ins + total_del) > 0) (sum(ins == 1) + sum(del == 1)) / (total_ins + total_del) * 100 else NA,
    stringsAsFactors = FALSE
  )
}))

summary_df <- merge(summary_df, summary_indels, by="sample", all.x=TRUE)
summary_df$target <- experiment

#### SAVE RESULTS ####
# Nextflow automatically captures these because they are written to the current working directory.
write.table(summary_df, file = "summary_df.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(insseq_merged, file = "insseq_merged.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

varcounts_del <- varcounts[varcounts$byType.g1=="deletion",]
indel_strings_del <- rownames(varcounts_del)
all_events_del <- list()

for (i in seq_along(indel_strings_del)) {
  s <- indel_strings_del[i]
  if (s == "no variant") next
  if (grepl("^SNV", s)) next
  
  events <- parse_indel_string(s)
  if (is.null(events) || nrow(events) == 0) next
  
  for (sample in samplenames) {
    n <- varcounts_del[i, sample]
    if (is.na(n)) next
    if (n > 0) {
      tmp <- events[rep(seq_len(nrow(events)), each=n), , drop=FALSE]
      tmp$sample <- sample
      all_events_del[[length(all_events_del)+1]] <- tmp
    }
  }
}

if (length(all_events_del) > 0) {
  all_events_del <- do.call(rbind, all_events_del)
  write.table(all_events_del, file = "all_events_del.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cat("No deletion events found to write.\n")
}