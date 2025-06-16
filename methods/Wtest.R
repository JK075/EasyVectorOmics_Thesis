# === Load libraries ===
library(dplyr)
library(readr)
library(tm)
library(wordcloud)
library(RColorBrewer)
library(purrr)

# === Function to create GO wordclouds for multiple ontologies ===
create_GO_wordclouds <- function(distances_file, 
                                annotations_file, 
                                ontologies = c("BP", "MF", "CC"), 
                                output_dir = "results/wordclouds") {
  
  # Initialize logfile
  logfile_path <- file.path(output_dir, "wordcloud_analysis.log")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  cat("Wordcloud Analysis Log\n", file = logfile_path)
  cat("======================\n\n", file = logfile_path, append = TRUE)
  
  # === Load and filter data ===
  # Read distance data and GO annotations
  distances_df <- read.table(distances_file, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  go_annotations <- read.table(annotations_file, header = TRUE, sep = "\t", 
                             quote = "", stringsAsFactors = FALSE)
  
  # Get extreme outliers from distance data
  extreme_outliers <- distances_df %>%
    filter(outlier_of_outlier == TRUE) %>%
    pull(gene_id) %>%
    unique()
  
  # Get non-outliers (all genes not in extreme_outliers)
  non_outliers <- distances_df %>%
    filter(!gene_id %in% extreme_outliers) %>%
    pull(gene_id) %>%
    unique()
  
  # Log basic information
  cat(paste("Analysis started at:", Sys.time(), "\n"), file = logfile_path, append = TRUE)
  cat(paste("Number of extreme outliers:", length(extreme_outliers), "\n"), 
      file = logfile_path, append = TRUE)
  cat(paste("Number of non-outliers:", length(non_outliers), "\n\n"), 
      file = logfile_path, append = TRUE)
  
  # === Process each ontology ===
  walk(ontologies, function(ontology) {
    
    # Filter annotations for current ontology and extreme outliers
    filtered_go <- go_annotations %>%
      filter(qpid %in% extreme_outliers & ontology == !!ontology) %>%
      distinct(qpid, desc, goid)
    
    # Filter annotations for current ontology and non-outliers (background)
    background_go <- go_annotations %>%
      filter(qpid %in% non_outliers & ontology == !!ontology) %>%
      distinct(qpid, desc, goid)
    
    # Skip if no annotations found for this ontology
    if (nrow(filtered_go) == 0) {
      message(paste("No annotations found for ontology:", ontology))
      cat(paste("\nNo annotations found for ontology:", ontology, "\n"), 
          file = logfile_path, append = TRUE)
      return()
    }
    
    # Log top 15 GO IDs and their frequencies
    top_go_ids <- filtered_go %>%
      count(goid, desc, sort = TRUE) %>%
      head(15)
    
    cat(paste("\n=== Ontology:", ontology, "===\n"), file = logfile_path, append = TRUE)
    cat("\nTop 15 GO IDs and their frequencies:\n", file = logfile_path, append = TRUE)
    capture.output(print(top_go_ids, row.names = FALSE), file = logfile_path, append = TRUE)
    
    # Text processing pipeline
    corpus <- Corpus(VectorSource(filtered_go$desc))
    
    # Apply text transformations without any stopword removal
    corpus <- corpus %>%
      tm_map(content_transformer(tolower)) %>%
      tm_map(removePunctuation) %>%
      tm_map(removeNumbers) %>%
      tm_map(stripWhitespace)
    
    # Create term-document matrix and calculate word frequencies
    tdm <- TermDocumentMatrix(corpus)
    word_freq <- sort(rowSums(as.matrix(tdm)), decreasing = TRUE)
    
    # Log word frequency information
    cat(paste("\nNumber of unique words in wordcloud:", length(word_freq), "\n"), 
        file = logfile_path, append = TRUE)
    cat("\nTop 10 words in wordcloud:\n", file = logfile_path, append = TRUE)
    top_words <- head(data.frame(Word = names(word_freq), Frequency = word_freq), 10)
    capture.output(print(top_words, row.names = FALSE), file = logfile_path, append = TRUE)
    
    # Define output file path
    output_file <- file.path(output_dir, paste0("wordcloud_", ontology, ".svg"))
    
    # Generate and save wordcloud as SVG
    svg(output_file, width = 10, height = 8)
    wordcloud(names(word_freq), word_freq,
             scale = c(4, 0.5),
             max.words = 100,
             random.order = FALSE,
             rot.per = 0.25,
             colors = brewer.pal(8, "Dark2"))
    title(paste("GO:", ontology, "for Extreme Outliers"))
    dev.off()
    
    # Perform Fisher's exact test for word enrichment
    perform_fisher_test(filtered_go, background_go, ontology, output_dir, logfile_path)
    
    # Print status messages
    message(paste("Created wordcloud for", ontology, "with", length(word_freq), "terms"))
    message(paste("Saved to:", output_file))
  })
  
  # Create combined wordcloud if multiple ontologies specified
  if (length(ontologies) > 1) {
    create_combined_wordcloud(go_annotations, extreme_outliers, non_outliers, ontologies, output_dir, logfile_path)
  }
  
  cat(paste("\nAnalysis completed at:", Sys.time(), "\n"), file = logfile_path, append = TRUE)
  message(paste("\nLogfile saved to:", logfile_path))
}

# Function to perform Fisher's exact test for word enrichment
perform_fisher_test <- function(outlier_go, background_go, ontology, output_dir, logfile_path) {
  # Combine all descriptions
  outlier_descriptions <- outlier_go$desc
  background_descriptions <- background_go$desc
  
  # Create corpus for outliers
  outlier_corpus <- Corpus(VectorSource(outlier_descriptions)) %>%
    tm_map(content_transformer(tolower)) %>%
    tm_map(removePunctuation) %>%
    tm_map(removeNumbers) %>%
    tm_map(stripWhitespace)
  
  # Create corpus for background
  background_corpus <- Corpus(VectorSource(background_descriptions)) %>%
    tm_map(content_transformer(tolower)) %>%
    tm_map(removePunctuation) %>%
    tm_map(removeNumbers) %>%
    tm_map(stripWhitespace)
  
  # Create TDM for both
  outlier_tdm <- TermDocumentMatrix(outlier_corpus)
  background_tdm <- TermDocumentMatrix(background_corpus)
  
  # Get word frequencies
  outlier_words <- rowSums(as.matrix(outlier_tdm))
  background_words <- rowSums(as.matrix(background_tdm))
  
  # Get all unique words
  all_words <- unique(c(names(outlier_words), names(background_words)))
  
  # Prepare data for Fisher's test
  fisher_results <- map_dfr(all_words, function(word) {
    # Counts for the word
    a <- ifelse(word %in% names(outlier_words), outlier_words[word], 0)  # word in outliers
    b <- length(outlier_descriptions) - a  # word not in outliers
    c <- ifelse(word %in% names(background_words), background_words[word], 0)  # word in background
    d <- length(background_descriptions) - c  # word not in background
    
    # Skip words with very low counts
    if ((a + c) < 3) return(NULL)
    
    # Fisher's exact test
    fisher_test <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    
    data.frame(
      word = word,
      p_value = fisher_test$p.value,
      odds_ratio = fisher_test$estimate,
      outlier_count = a,
      background_count = c,
      total_outliers = length(outlier_descriptions),
      total_background = length(background_descriptions)
    )
  })
  
  # Adjust for multiple testing
  fisher_results$adjusted_p <- p.adjust(fisher_results$p_value, method = "BH")
  
  # Sort by significance
  fisher_results <- fisher_results %>%
    arrange(adjusted_p, p_value)
  
  # Save results to file
  fisher_file <- file.path(output_dir, paste0("fisher_test_results_", ontology, ".tsv"))
  write_tsv(fisher_results, fisher_file)
  
  # Log top significant results
  cat("\nFisher's exact test results (top 20 significant words):\n", file = logfile_path, append = TRUE)
  top_results <- fisher_results %>%
    filter(adjusted_p < 0.05) %>%
    head(20)
  
  if (nrow(top_results) > 0) {
    capture.output(print(top_results, row.names = FALSE), file = logfile_path, append = TRUE)
  } else {
    cat("No significantly enriched words found at adjusted p-value < 0.05\n", 
        file = logfile_path, append = TRUE)
  }
  
  # Create volcano plot of results
  volcano_file <- file.path(output_dir, paste0("fisher_volcano_", ontology, ".svg"))
  svg(volcano_file, width = 8, height = 6)
  
  plot_data <- fisher_results %>%
    mutate(
      log_odds = log2(odds_ratio),
      significance = -log10(adjusted_p),
      is_significant = adjusted_p < 0.05
    )
  
  plot(plot_data$log_odds, plot_data$significance,
       xlab = "Log2 Odds Ratio", ylab = "-Log10(Adjusted p-value)",
       main = paste("Word Enrichment in Outliers -", ontology),
       col = ifelse(plot_data$is_significant, "red", "gray50"),
       pch = 19)
  abline(h = -log10(0.05), col = "blue", lty = 2)
  abline(v = 0, col = "black", lty = 1)
  grid()
  
  # Label top 5 significant points
  if (nrow(top_results) >= 5) {
    top5 <- head(plot_data[plot_data$is_significant, ], 5)
    text(top5$log_odds, top5$significance, labels = top5$word, pos = 3, cex = 0.7)
  }
  
  dev.off()
  
  message(paste("Fisher's test results saved to:", fisher_file))
  message(paste("Volcano plot saved to:", volcano_file))
}

# Helper function for combined wordcloud across ontologies
create_combined_wordcloud <- function(go_annotations, extreme_outliers, non_outliers, ontologies, output_dir, logfile_path) {
  # Filter annotations for all specified ontologies and outliers
  filtered_go <- go_annotations %>%
    filter(qpid %in% extreme_outliers & ontology %in% ontologies) %>%
    distinct(qpid, desc, goid)
  
  # Filter background for all specified ontologies and non-outliers
  background_go <- go_annotations %>%
    filter(qpid %in% non_outliers & ontology %in% ontologies) %>%
    distinct(qpid, desc, goid)
  
  # Return if no annotations found
  if (nrow(filtered_go) == 0) return()
  
  # Log top 15 GO IDs and their frequencies for combined analysis
  top_go_ids <- filtered_go %>%
    count(goid, desc, sort = TRUE) %>%
    head(15)
  
  cat(paste("\n=== Combined Ontologies:", paste(ontologies, collapse = "+"), "===\n"), 
      file = logfile_path, append = TRUE)
  cat("\nTop 15 GO IDs and their frequencies:\n", file = logfile_path, append = TRUE)
  capture.output(print(top_go_ids, row.names = FALSE), file = logfile_path, append = TRUE)
  
  # Text processing pipeline
  corpus <- Corpus(VectorSource(filtered_go$desc))
  
  corpus <- corpus %>%
    tm_map(content_transformer(tolower)) %>%
    tm_map(removePunctuation) %>%
    tm_map(removeNumbers) %>%
    tm_map(stripWhitespace)
  
  # Create term-document matrix and calculate word frequencies
  tdm <- TermDocumentMatrix(corpus)
  word_freq <- sort(rowSums(as.matrix(tdm)), decreasing = TRUE)
  
  # Log word frequency information for combined analysis
  cat(paste("\nNumber of unique words in combined wordcloud:", length(word_freq), "\n"), 
      file = logfile_path, append = TRUE)
  cat("\nTop 10 words in combined wordcloud:\n", file = logfile_path, append = TRUE)
  top_words <- head(data.frame(Word = names(word_freq), Frequency = word_freq), 10)
  capture.output(print(top_words, row.names = FALSE), file = logfile_path, append = TRUE)
  
  # Define output file path for combined wordcloud
  output_file <- file.path(output_dir, 
                         paste0("wordcloud_combined_", 
                                paste(ontologies, collapse = "_"), ".svg"))
  
  # Generate and save combined wordcloud
  svg(output_file, width = 10, height = 8)
  wordcloud(names(word_freq), word_freq,
           scale = c(4, 0.5),
           max.words = 100,
           random.order = FALSE,
           rot.per = 0.2,
           colors = brewer.pal(8, "Dark2"))
  title(paste("GO:", paste(ontologies, collapse = "+"), "for Extreme Outliers"))
  dev.off()
  
  # Perform Fisher's test for combined ontologies
  perform_fisher_test(filtered_go, background_go, 
                     paste(ontologies, collapse = "_"), output_dir, logfile_path)
  
  # Print status messages
  message(paste("\nCreated COMBINED wordcloud for", paste(ontologies, collapse = "+")))
  message(paste("Saved to:", output_file))
}

# === Example usage ===
create_GO_wordclouds(
  distances_file = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/EVO/results2/majority/ttest/distances_df.tsv",
  annotations_file = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/EVO/results2/pannzer2_results.tsv",
  ontologies = c("BP", "CC", "MF"),
  output_dir = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/EVO/results2/majority/ttest/"
)
