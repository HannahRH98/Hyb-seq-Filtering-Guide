# Original script developed by N. Meeprom (https://orcid.org/0000-0003-4193-7062) but further adapted by H. Hall
# This script is designed to analyse HybPiper output statistics ("hybpiper_stats.tsv")
# and assist in deciding filtering strategies based on gene/supercontig recovery and length across samples.
# It provides visualisations and numerical summaries showing how many samples
# meet certain gene/supercontig recovery and length recovery thresholds (geneswithseqs, 25%, 50%, 75%, 150% (especially when supercontigs are called) - standard hypiper_stats output).
# Users can explore sample retention across a range of gene/supercontig length recovery cutoffs.
# Can help determine the threshold of data the user takes onto later downstream analyses.
# This should be used as a guide: removal decisions should also consider sample importance,
# not only gene recovery counts. Some users may decide not to filter samples at all.
# REMEMBER --> change gene to supercontig or vice versa when needed.

library(tidyverse)
library(ggplot2)
library(cowplot)

# Define label for plots
dataset <- element_text("loci")

# Load HybPiper stats file output; this file should contain per-sample gene recovery stats
stats <- read.table("duplicates_removed_hybpiper_stat.tsv", sep = "\t", header = TRUE)

# Display first few rows to verify data loaded correctly
head(stats)

## Calculate average number of genes/supercontigs recovered at different length recovery thresholds across all samples
gene_recovery_summary <- summarise(stats,
                                   mean_WithSeqs = mean(GenesWithSeqs, na.rm = TRUE),
                                   mean_25 = mean(GenesAt25pct, na.rm = TRUE),
                                   mean_50 = mean(GenesAt50pct, na.rm = TRUE),
                                   mean_75 = mean(GenesAt75pct, na.rm = TRUE),
                                   mean_150 = if ("GenesAt150pct" %in% colnames(stats)) mean(GenesAt150pct, na.rm = TRUE) else NA_real_
)

# when calling hybpiper stats set to gene, the "GenesAt150pct" column should contain 0s, but when called with supercontigs, this can contain values.


print("Average loci recovered per sample at different sequence recovery thresholds:")
print(gene_recovery_summary)

# Define the range of loci thresholds to examine visually
# IMPORTANT!!!!! This is set to 0-353 for Angiosperms353 bait kit by default; modify if different kit used
genes <- 0:353

# Initialize output dataframe to store counts of samples meeting thresholds
outputs <- data.frame()

# Loop over each gene/supercontig recovery cutoff threshold to calculate how many samples meet that threshold
for (i in 1:length(genes)){
  count_withseqs <- sum(stats$GenesWithSeqs >= genes[i])
  count_25 <- sum(stats$GenesAt25pct >= genes[i])
  count_50 <- sum(stats$GenesAt50pct >= genes[i])
  count_75 <- sum(stats$GenesAt75pct >= genes[i])
  count_150 <- if ("GenesAt150pct" %in% colnames(stats)) sum(stats$GenesAt150pct >= genes[i]) else NA_real_
  
  outputs[i,1] <- genes[i]
  outputs[i,2] <- (genes[i] * 100/length(genes))
  outputs[i,3] <- count_withseqs
  outputs[i,4] <- count_25
  outputs[i,5] <- count_50
  outputs[i,6] <- count_75
  outputs[i,7] <- count_150
}

# Rename output columns with shorter names
colnames(outputs) <- c("Loci", "PercentLoci", "WithSeqs", "At25pct", "At50pct", "At75pct", "At150pct")

# View the resulting summary table
head(outputs)

# Convert from wide to long format for ggplot2 visualization
threshold_cols <- c("WithSeqs", "At25pct", "At50pct", "At75pct")
if (!all(is.na(outputs$At150pct))) {
  threshold_cols <- c(threshold_cols, "At150pct")
}

outputs_longer <- pivot_longer(outputs,
                               cols = all_of(threshold_cols),
                               names_to = "ThresholdType",
                               values_to = "Samples")

# Plot number of samples retained across loci recovery thresholds, separated by threshold type
samples_retained_plot <- ggplot(outputs_longer, aes(x = Loci, y = Samples, color = ThresholdType)) + 
  geom_line(size = 1) + 
  theme_minimal() +
  labs(title = "Retention Across Sequence Recovery Thresholds",
       x = "Number of Loci Recovered",
       y = "Number of Samples Retained",
       color = "Length Recovery Threshold")

print(samples_retained_plot)

# Save plot as PNG
ggsave("samples_retained_supercontig.png", samples_retained_plot, dpi = 300, width = 8, height = 5)

# You may want to filter based on where the number of recovered loci begins to drop off OR drops off significantly.
# OR you may want to look at the maximum number of samples you can retain without losing to many loci
# i.e. where the number of samples you wish to retain hits the threshold lines
# Some users care more about sequence length recovery thresholds than just the total number of loci with sequences.
# Recovery thresholds (e.g., 25%, 50%, 75%) can greatly influence how many loci are counted as successfully recovered per sample.
# However, if recovery is generally good, there may be little difference between thresholds,
# and the user can treat them similarly. The key is to balance the trade-off between retaining samples
# recovering more loci, and obtaining longer sequences.
# REMEMBER, when calling supercontigs, you are MUCH more likely to get longer "genes" as these include not only
# the exon but introns/flanking regions too!
# IT IS ALL ABOUT BALANCE!

## Calculate slopes (rate of change) of sample retention as loci threshold increases
slope_outputs <- data.frame(
  Loci = genes,
  SlopeWithSeqs = NA_real_,
  Slope25 = NA_real_,
  Slope50 = NA_real_,
  Slope75 = NA_real_,
  Slope150 = NA_real_
)

for (i in 2:length(genes)) {
  slope_outputs$SlopeWithSeqs[i] <- outputs$WithSeqs[i] - outputs$WithSeqs[i - 1]
  slope_outputs$Slope25[i] <- outputs$At25pct[i] - outputs$At25pct[i - 1]
  slope_outputs$Slope50[i] <- outputs$At50pct[i] - outputs$At50pct[i - 1]
  slope_outputs$Slope75[i] <- outputs$At75pct[i] - outputs$At75pct[i - 1]
  if (!all(is.na(outputs$At150pct))) {
    slope_outputs$Slope150[i] <- outputs$At150pct[i] - outputs$At150pct[i - 1]
  }
}

# Reshape slope data for plotting
slope_outputs_longer <- pivot_longer(slope_outputs,
                                    cols = c("SlopeWithSeqs", "Slope25", "Slope50", "Slope75", "Slope150"),
                                    names_to = "ThresholdType",
                                    values_to = "Slope")

# Plot slope of sample retention change across loci length thresholds
slope_plot <- ggplot(slope_outputs_longer, aes(x = Loci, y = Slope, color = ThresholdType)) + 
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Slope of Sample Retention Change Across Loci Thresholds",
       x = "Number of Loci Recovered (Gene Threshold)",
       y = "Slope (Samples Lost per Locus Increase)",
       color = "Recovery Threshold")

print(slope_plot)

# Combine plots vertically
combined_plots <- plot_grid(samples_retained_plot, slope_plot, ncol = 1, align = "v")

# Save combined plot
ggsave("samples_loci_slope.png", combined_plots, dpi = 300, width = 8, height = 10)

## Join loci and slope data for export and further analysis
joined_loci_slope <- left_join(outputs, slope_outputs, by = "Loci")

write_csv(joined_loci_slope, "joined_loci_slope.csv")


