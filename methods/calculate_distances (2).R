library(rstatix)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggpubr)
require(RColorBrewer)

# ====== Global Settings ======
# Define statistical test settings
TEST_METHOD <- "wilcox.test"  # Choose between "t.test" or "wilcox.test"
TEST_ALTERNATIVE <- "greater" # Always use "greater" alternative hypothesis

addAlpha <- function(col, alpha = .25) {
  apply(sapply(col, col2rgb) / 255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha))
}

colors <- brewer.pal(7, "Dark2")
alpha.cols <- sapply(colors, addAlpha)


# Load expression data (gene_id + expression values)
norm_df <- read.table("/storage/EasyVectorOmics/FastQ_GSE125483_JK/kallisto_outputs/combined_output.tsv", header = TRUE, sep = "\t")

# Load centroid data (Orthogroup + centroid vectors)
centroids_df <- read.table("/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/EVO/results/orthogroup_centroids_majorityrule.tsv", header = TRUE, sep = "\t")

# Load orthogroup mapping (Orthogroup ↔ genes)
orthogroups_df <- read.table("/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/Orthogroups_filtered.tsv", header = TRUE, sep = "\t", quote = "", fill = TRUE)

orthologs_df <- read.table("/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_conserved_orthologs_2_formatted_filtered.tsv", header = TRUE, sep = "\t")  # Must have columns: Orthogroup, Gene1, Gene2
paralogs_df <- read.table("/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_special_in_and_out_paralogs_2_formatted_filtered.tsv", header = TRUE, sep = "\t")


colnames(norm_df)[1] <- "gene_id"

# === Convert orthogroups to long format ===
orthogroups_long <- orthogroups_df %>%
  pivot_longer(-Orthogroup, names_to = "species", values_to = "genes") %>%
  filter(!is.na(genes)) %>%
  mutate(genes = strsplit(as.character(genes), ",\\s*")) %>%
  unnest(genes) %>%
  rename(gene_id = genes)

# === Identify expression columns dynamically ===
expr_cols <- setdiff(names(norm_df), "gene_id")

# === Create list of centroid vectors by orthogroup ===
centroid_list <- lapply(split(centroids_df, centroids_df$Orthogroup), function(df) {
  as.numeric(df[, expr_cols])
})

# === Join orthogroup info to expression table ===
norm_with_og <- norm_df %>%
  inner_join(orthogroups_long, by = "gene_id")

# === Compute distances ===
compute_distance <- function(row) {
  gene_vector <- as.numeric(row[expr_cols])
  og <- row[["Orthogroup"]]
  centroid_vector <- as.numeric(centroids_df[centroids_df$Orthogroup == og, expr_cols])
  sqrt(sum((gene_vector - centroid_vector)^2))
}

# Apply function row-wise
norm_with_og$distance_to_centroid <- apply(norm_with_og, 1, compute_distance)

# === Output ===
head(norm_with_og)

# === Add gene_type label (ortholog or paralog) ===
ortholog_gene_map <- orthologs_df %>%
  select(Orthogroup, Gene1, Gene2) %>%
  pivot_longer(cols = c(Gene1, Gene2), names_to = NULL, values_to = "gene_id") %>%
  distinct()  # remove duplicates

ortholog_genes <- ortholog_gene_map$gene_id
paralog_genes <- unique(c(paralogs_df$Gene1, paralogs_df$Gene2))

norm_with_og <- norm_with_og %>%
  mutate(gene_type = case_when(
    gene_id %in% ortholog_genes ~ "ortholog",
    gene_id %in% paralog_genes ~ "paralog",
    TRUE ~ NA_character_
  ))


norm_df_1 <- norm_df %>%
  rename_with(.fn = ~ paste0(.x, "_1"), .cols = all_of(expr_cols))

norm_df_2 <- norm_df %>%
  rename_with(.fn = ~ paste0(.x, "_2"), .cols = all_of(expr_cols))

# --- Join expression data for Gene1 and Gene2 ---
orthologs_with_fc <- orthologs_df %>%
  left_join(norm_df_1, by = c("Gene1" = "gene_id")) %>%
  left_join(norm_df_2, by = c("Gene2" = "gene_id"))

print(head(orthologs_with_fc))

orthologs_with_fc <- orthologs_with_fc %>%
  rowwise() %>%
  mutate(distance = sqrt(sum((c_across(all_of(paste0(expr_cols, "_1"))) -
                              c_across(all_of(paste0(expr_cols, "_2"))))^2, na.rm = TRUE))) %>%
  ungroup()

# --- Get the maximum ortholog distance per orthogroup ---
ortholog_max_pairwise <- orthologs_with_fc %>%
  group_by(Orthogroup) %>%
  summarise(max_ortholog_pairwise_distance = max(distance, na.rm = TRUE), .groups = "drop")

# Preview the result
head(ortholog_max_pairwise)

# Join max pairwise ortholog distance and compute RDI
distances_with_rdi <- norm_with_og %>%
  left_join(ortholog_max_pairwise, by = "Orthogroup") %>%
  mutate(RDI = distance_to_centroid / max_ortholog_pairwise_distance)


# Calculate the 95th percentile of all RDI values (empirical threshold)
p95 <- quantile(distances_with_rdi$RDI[!is.na(distances_with_rdi$RDI) & distances_with_rdi$RDI < Inf], 0.95, na.rm = TRUE)
print(p95)

# Mark outliers
distances_with_rdi <- distances_with_rdi %>%
  mutate(is_outlier = RDI > p95)

get_gene_type <- function(gene_id) {
  if (gene_id %in% ortholog_genes) {
    return("ortholog")
  } else if (gene_id %in% paralog_genes) {
    return("paralog")
  } else {
    return(NA)  
  }
}

distances_with_rdi$gene_type <- sapply(distances_with_rdi$gene_id, get_gene_type)

summary(distances_with_rdi$RDI)

# --- New code for 'outlier_of_outlier' ---
outliers_rdi <- distances_with_rdi %>%
  filter(is_outlier == TRUE)


p95_outlier <- quantile(outliers_rdi$RDI[!is.na(outliers_rdi$RDI) & outliers_rdi$RDI < Inf], 0.95, na.rm = TRUE)
print(p95_outlier)


outliers_rdi <- outliers_rdi %>%
  mutate(outlier_of_outlier = RDI > p95_outlier)


distances_with_rdi <- distances_with_rdi %>%
  left_join(outliers_rdi %>% select(gene_id, outlier_of_outlier), by = "gene_id") %>%
  mutate(outlier_of_outlier = ifelse(is.na(outlier_of_outlier), FALSE, outlier_of_outlier))


# Plot all distances to centroid
all_distances <- ggplot(distances_with_rdi, aes(x = gene_type, y = distance_to_centroid, fill = gene_type)) +
geom_boxplot(outlier.shape = NA, alpha = 0.25, aes(color = gene_type, fill = gene_type)) +  
geom_jitter(width = 0.2, shape = 21, size = 2, aes(color = gene_type, fill = gene_type), alpha = 0.25) +
  scale_fill_manual(values = c(
    "ortholog" = colors[[1]],
    "paralog" = colors[[2]],
    "source-copy" = colors[[3]]
  )) +
  scale_color_manual(values = c(
    "ortholog" = colors[[1]],
    "paralog" = colors[[2]],
    "source-copy" = colors[[3]]
  )) +  labs(
    title = paste("Distance to Centroid by Gene Type Tree"),
    x = "Gene Type",
    y = "Distance to Centroid"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # Zentrierung
  ) +
  stat_compare_means(
    method = TEST_METHOD,
    comparisons = list(c("paralog", "ortholog")),
    alternative = TEST_ALTERNATIVE,  
    label = "p.signif",
    label.y = max(distances_with_rdi$distance_to_centroid, na.rm = TRUE) + 0.1
  )

ggsave("results2/all_distances.svg", plot = all_distances, width = 8, height = 6, device = svg)

effsize_results <- distances_with_rdi %>%
  wilcox_effsize(
    formula = distance_to_centroid ~ gene_type,
    comparisons = list(c("paralog", "ortholog")),
    alternative = TEST_ALTERNATIVE,
    ci = TRUE
  )

# In eine TXT-Datei speichern
write.table(
  effsize_results,
  file = "results2/wilcox_effsize.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Plot distribution of RDI values: distance to centroid
rdi_plot <- ggplot(distances_with_rdi, aes(x = RDI)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  geom_vline(xintercept = p95, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Distribution of Relative Divergence Index (RDI): Distance to centroid",
    x = "RDI",
    y = "Count"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # Zentrierung
  ) +
  annotate("text", x = p95, y = Inf, label = paste0("Threshold: ", round(p95, 3)),
           vjust = -0.5, hjust = 1, color = "red", size = 4)

ggsave("results2/rdi_distribution_plot.svg", plot = rdi_plot, width = 8, height = 6, device = svg)


# Plot outlier distances to centroid
outlier_distances <- ggplot(
  filter(distances_with_rdi, is_outlier == TRUE),
  aes(x = gene_type, y = distance_to_centroid, fill = gene_type)
) +
geom_boxplot(outlier.shape = NA, alpha = 0.25, aes(color = gene_type, fill = gene_type)) +  
geom_jitter(width = 0.2, shape = 21, size = 2, aes(color = gene_type, fill = gene_type), alpha = 0.25) +
  scale_fill_manual(values = c(
    "ortholog" = colors[[1]],
    "paralog" = colors[[2]],
    "source-copy" = colors[[3]]
  )) +
  scale_color_manual(values = c(
    "ortholog" = colors[[1]],
    "paralog" = colors[[2]],
    "source-copy" = colors[[3]]
  )) +  labs(
    title = paste("Distance to Centroid by Gene Type (Outliers Only) Tree (Majority)", TEST_METHOD),
    x = "Gene Type",
    y = "Distance to Centroid"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # Zentrierung
  ) +
  stat_compare_means(
    method = TEST_METHOD,
    comparisons = list(c("paralog", "ortholog")),
    alternative = TEST_ALTERNATIVE,  
    label = "p.signif",
    label.y = max(distances_with_rdi$distance_to_centroid, na.rm = TRUE) + 0.1
  )

ggsave("results2/outlier_distances_t_test.svg", plot = outlier_distances, width = 8, height = 6, device = svg)

write.table(distances_with_rdi, file = "results2/distances_df.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot outlier distances to centroid per family (mean)
outlier_means <- distances_with_rdi %>%
  filter(is_outlier == TRUE) %>%  
  group_by(Orthogroup, gene_type) %>%
  summarise(mean_distance = mean(distance_to_centroid, na.rm = TRUE), .groups = "drop")

outlier_dist_per_family <- ggplot(outlier_means, aes(x = gene_type, y = mean_distance, fill = gene_type)) +
geom_boxplot(outlier.shape = NA, alpha = 0.25, aes(color = gene_type, fill = gene_type)) +  
geom_jitter(width = 0.2, shape = 21, size = 2, aes(color = gene_type, fill = gene_type), alpha = 0.25) +
  scale_fill_manual(values = c(
    "ortholog" = colors[[1]],
    "paralog" = colors[[2]],
    "source-copy" = colors[[3]]
  )) +
  scale_color_manual(values = c(
    "ortholog" = colors[[1]],
    "paralog" = colors[[2]],
    "source-copy" = colors[[3]]
  )) +  labs(
    title = paste("Mean Distance to Centroid per Family (Outliers Only) Tree"),
    x = "Gene Type",
    y = "Distance to Centroid"
  ) +
  theme_classic() +
  theme(legend.position = "none")+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # Zentrierung
  ) +
  stat_compare_means(
    method = TEST_METHOD,
    comparisons = list(c("paralog", "ortholog")),
    alternative = TEST_ALTERNATIVE,  
    label = "p.signif",
    label.y = max(distances_with_rdi$distance_to_centroid, na.rm = TRUE) + 0.1
  )

ggsave("results2/outlier_distances_per_family.svg", plot = outlier_dist_per_family, width = 8, height = 6, device = svg)


effsize_outlier_means <- outlier_means %>%
  wilcox_effsize(
    formula = mean_distance ~ gene_type,
    comparisons = list(c("paralog", "ortholog")),
    alternative = TEST_ALTERNATIVE,
    ci = TRUE
  )

# Ergebnisse an bestehende Datei anhängen
write("\n\n# Effektstärken für outlier_means (Mean Distance per Family)", 
      file = "results2/wilcox_effsize.txt", append = TRUE)

write.table(
  effsize_outlier_means,
  file = "results2/wilcox_effsize.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  append = TRUE
)

# Get source-copy distances

# Prepare norm_df with suffixes for gene1 and gene2
norm_df_1 <- norm_df %>%
  rename_with(.fn = ~ paste0(.x, "_1"), .cols = all_of(expr_cols))

norm_df_2 <- norm_df %>%
  rename_with(.fn = ~ paste0(.x, "_2"), .cols = all_of(expr_cols))

# Join expression vectors for each gene
paralog_distances <- paralogs_df %>%
  select(Orthogroup, Gene1, Gene2) %>%
# Füge Species-Information für Gene1 hinzu (aus orthogroups_long)
  left_join(orthogroups_long %>% select(gene_id, species) %>%
             rename(Gene1_species = species),
           by = c("Gene1" = "gene_id")) %>%
  # Füge Species-Information für Gene2 hinzu (aus orthogroups_long)
  left_join(orthogroups_long %>% select(gene_id, species) %>%
             rename(Gene2_species = species),
           by = c("Gene2" = "gene_id")) %>%
  # Join expression data
  left_join(norm_df_1, by = c("Gene1" = "gene_id")) %>%
  left_join(norm_df_2, by = c("Gene2" = "gene_id")) %>%
  # Filter out missing values in any vector
  filter(if_all(all_of(paste0(expr_cols, "_1")), ~ !is.na(.)) &
         if_all(all_of(paste0(expr_cols, "_2")), ~ !is.na(.))) %>%
  # Calculate Euclidean distance between gene1 and gene2
  rowwise() %>%
  mutate(distance = sqrt(sum((c_across(all_of(paste0(expr_cols, "_1"))) -
                              c_across(all_of(paste0(expr_cols, "_2"))))^2))) %>%
  ungroup() %>%
  mutate(
    gene1_type = case_when(
      Gene1 %in% ortholog_genes ~ "ortholog",
      Gene1 %in% paralog_genes ~ "paralog",
      TRUE ~ NA_character_
    ),
    gene2_type = case_when(
      Gene2 %in% ortholog_genes ~ "ortholog",
      Gene2 %in% paralog_genes ~ "paralog",
      TRUE ~ NA_character_
    )
  )

paralog_distances <- paralog_distances %>%
  left_join(ortholog_max_pairwise, by = "Orthogroup") %>%
  mutate(RDI = distance / max_ortholog_pairwise_distance)

head(paralog_distances)
# Calculate the 95th percentile of all RDI values (empirical threshold)
p95 <- quantile(paralog_distances$RDI[!is.na(paralog_distances$RDI) & paralog_distances$RDI < Inf], 0.95, na.rm = TRUE)

print(p95)

# Mark outliers
paralog_distances <- paralog_distances %>%
  mutate(is_outlier = RDI > p95)

# Jetzt: Outlier-der-Outlier (Top 5% der Outlier)
if (sum(paralog_distances$is_outlier, na.rm = TRUE) > 0) {
  # Filtere nur Outlier und berechne das 95. Perzentil IHRER RDI-Werte
  p95_outlier_of_outlier <- quantile(
    paralog_distances$RDI[paralog_distances$is_outlier & !is.na(paralog_distances$RDI)],
    0.95,
    na.rm = TRUE
  )

  # Markiere die extremsten 5% der Outlier
  paralog_distances <- paralog_distances %>%
    mutate(
      outlier_of_outlier = ifelse(
        is_outlier & RDI > p95_outlier_of_outlier & !is.na(RDI),
        TRUE,
        FALSE
      )
    )
} else {
  # Falls keine Outlier existieren, setze alle outlier_of_outlier auf FALSE
  paralog_distances$outlier_of_outlier <- FALSE
}

paralog_distances <- paralog_distances %>%
  mutate(relationship = case_when(
    Gene1_species == Gene2_species & gene1_type == "ortholog" & gene2_type == "ortholog" ~ "SS",
    Gene1_species == Gene2_species & (gene1_type == "paralog" | gene2_type == "paralog") ~ "SIN",
    Gene1_species != Gene2_species & (gene1_type == "paralog" | gene2_type == "paralog") ~ "SOU",
    Gene1_species != Gene2_species & gene1_type == "ortholog" & gene2_type == "ortholog" ~ "OSS"
  ))


# Save results
write.table(paralog_distances, file = "results2/paralog_pairwise_distances.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot distribution of RDI values: source-copy
rdi_plot <- ggplot(paralog_distances, aes(x = RDI)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  geom_vline(xintercept = p95, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Distribution of Relative Divergence Index (RDI): Source-Copy",
    x = "RDI",
    y = "Count"
  ) +
  annotate("text", x = p95, y = Inf, label = paste0("Threshold: ", round(p95, 3)),
           vjust = -0.5, hjust = 1, color = "red", size = 4) +
  theme_classic()

ggsave("results2/rdi_distribution_source_copy.svg", plot = rdi_plot, width = 8, height = 6, device = svg)


# Source Copy Outliers in Inparalogs und Outparalogs aufteilen
source_copy_outliers <- paralog_distances %>%
  filter(is_outlier == TRUE) %>%
  filter(relationship %in% c("SIN", "SOU")) %>%  # Nur SIN und SOU behalten
  mutate(
    gene_type = case_when(
      relationship == "SIN" ~ "source-copy (inparalog)",
      relationship == "SOU" ~ "source-copy (outparalog)"
    ),
    distance_to_centroid = distance
  ) %>%
  select(gene_type, distance_to_centroid)

# Centroid Outliers (unverändert)
centroid_outliers <- distances_with_rdi %>%
  filter(is_outlier == TRUE) %>%
  select(gene_type, distance_to_centroid)

# Daten kombinieren
combined_outliers <- bind_rows(centroid_outliers, source_copy_outliers)

# Vergleichsgruppen definieren
comparisons <- list(
  c("ortholog", "paralog"),
  c("ortholog", "source-copy (inparalog)"),
  c("ortholog", "source-copy (outparalog)"),
  c("paralog", "source-copy (inparalog)"),
  c("paralog", "source-copy (outparalog)"),
  c("source-copy (inparalog)", "source-copy (outparalog)")
)

# Farben definieren (erweitert für die neuen Gruppen)
plot_colors <- c(
  "ortholog" = colors[[1]],
  "paralog" = colors[[2]],
  "source-copy (inparalog)" = colors[[3]],
  "source-copy (outparalog)" = colors[[4]]  # Annahme: colors hat mindestens 4 Elemente
)

# Plot erstellen
combined_plot <- ggplot(combined_outliers, aes(x = gene_type, y = distance_to_centroid, fill = gene_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, aes(color = gene_type, fill = gene_type)) +
  geom_jitter(width = 0.2, shape = 21, size = 2, aes(color = gene_type, fill = gene_type), alpha = 0.25) +
  scale_fill_manual(values = plot_colors) +
  scale_color_manual(values = plot_colors) +
  stat_compare_means(
    comparisons = comparisons,
    method = TEST_METHOD,
    alternative = TEST_ALTERNATIVE,
    label = "p.signif",
    label.y = max(combined_outliers$distance_to_centroid, na.rm = TRUE) * seq(1.05, 1.6, length.out = length(comparisons))
  ) +
  labs(
    title = paste("Distance by Gene Type (Outliers Only) Tree"),
    x = "Gene Type",
    y = "Distance"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)  # X-Achsenbeschriftung schräg stellen für bessere Lesbarkeit
  )

ggsave("results2/combined_outlier_distances.svg", plot = combined_plot, width = 10, height = 6, device = svg)


effsize_combined_outliers <- combined_outliers %>%
  wilcox_effsize(
    formula = distance_to_centroid ~ gene_type,
    comparisons = comparisons,
    alternative = TEST_ALTERNATIVE,
    ci = TRUE
  )

# Ergebnisse an bestehende Datei anhängen
write("\n\n# Effektstärken für combined_outliers", file = "results2/wilcox_effsize.txt", append = TRUE)
write.table(
  effsize_combined_outliers,
  file = "results2/wilcox_effsize.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  append = TRUE
)

# Define función para ejecutar tests
run_tests <- function(data, group1, group2) {
  subset_data <- data %>%
    filter(gene_type %in% c(group1, group2)) %>%
    mutate(gene_type = factor(gene_type, levels = c(group1, group2))) %>%  # Muy importante
    droplevels()

  t_test <- t.test(distance_to_centroid ~ gene_type, data = subset_data, alternative = "less")
  wilcox_test <- wilcox.test(distance_to_centroid ~ gene_type, data = subset_data, alternative = "less")

  data.frame(
    comparison = paste(group1, "vs", group2),
    method = c("t.test", "wilcox.test"),
    p_value = c(t_test$p.value, wilcox_test$p.value)
  )
}


# Comparaciones
comparisons <- list(
  c("ortholog", "paralog"),
  c("ortholog", "source-copy"),
  c("paralog", "source-copy")
)

# Ejecutar tests
results_list <- lapply(comparisons, function(comp) {
  run_tests(combined_outliers, comp[1], comp[2])
})

# Combinar y ajustar por BH
test_results <- do.call(rbind, results_list)
test_results$BH_adjusted <- p.adjust(test_results$p_value, method = "BH")

# Guardar resultados
write.table(test_results, "results2/p_values_BH_adjusted.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
