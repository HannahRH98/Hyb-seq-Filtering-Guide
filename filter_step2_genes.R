library(tidyverse)
library(ggplot2)
library(ggrepel)
library(cowplot)


# Load data
seq_lengths <- read.table("seq_lengths.tsv", sep = "\t", header = TRUE)


#check number of total genes present. For Angiosperms353 this should be 353.
cat("Number of genes detected:", ncol(seq_lengths) - (start_col - 1), "\n")


start_col <- 2
num_samples <- nrow(seq_lengths) - 1  # exclude header row (assuming first row is header)

# Calculate % Empty, % Present, and counts of present/absent per gene
percent_data <- seq_lengths[-1, start_col:ncol(seq_lengths)] %>%
  summarise(across(everything(), list(
    Empty = ~mean(. == 0) * 100,
    Present = ~mean(. > 0) * 100,
    PresentCount = ~sum(. > 0),
    AbsentCount = ~sum(. == 0)
  ))) %>%
  pivot_longer(everything(),
               names_to = c("Gene", ".value"),
               names_sep = "_") %>%
  arrange(Present)  # sort by Present ascending (low to high)

# Reorder columns so counts follow percentages side-by-side and rename columns
percent_data <- percent_data %>%
  select(Gene, Empty, AbsentCount, Present, PresentCount) %>%
  rename(
    `Empty (%)` = Empty,
    `Absent Count` = AbsentCount,
    `Present (%)` = Present,
    `Present Count` = PresentCount
  )

# Save combined summary file with counts and percentages
write.table(percent_data, file = "gene_missing_present_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Based upon the values generated in the "gene_missing_present_summary.tsv", a threshold can be decided 
# tested, and visualised in the histogram.
# User sets threshold here (minimum percent present required)
threshold_present <- 70  # e.g. keep genes with >= 70% present sequences

# Plot histogram showing distribution of % Present,
# coloring genes below threshold in red and others in blue
ggplot(percent_data, aes(x = `Present (%)`, fill = `Present (%)` < threshold_present)) +
  geom_histogram(bins = 30, color = "black") +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  geom_vline(xintercept = threshold_present, color = "red", linetype = "dashed") +
  labs(title = "Distribution of % Present Data per Gene",
       subtitle = paste("Genes with <", threshold_present, "% present shown in red"),
       x = "% Present per Gene",
       y = "Number of Genes",
       fill = "Below Threshold") +
  theme_minimal()

# Determine genes to remove and keep based on threshold
genes_to_remove <- percent_data %>% filter(`Present (%)` < threshold_present)
genes_to_keep <- percent_data %>% filter(`Present (%)` >= threshold_present)

# Summary output
cat("Filtering threshold (min % present required):", threshold_present, "%\n")
cat("Total genes:", nrow(percent_data), "\n")
cat("Genes to REMOVE:", nrow(genes_to_remove), "\n")
cat("Genes to KEEP:", nrow(genes_to_keep), "\n\n")

cat("Genes to REMOVE:\n")
print(genes_to_remove$Gene)

cat("\nGenes to KEEP:\n")
print(genes_to_keep$Gene)

# Save gene lists for downstream use if desired. These can be used as a genelist file remove undesired 
# genes from the dataset before running trees or can run trees for all of them and remove these gene trees later.
write.table(genes_to_remove$Gene, file = "genes_to_remove.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_to_keep$Gene, file = "genes_to_keep.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
