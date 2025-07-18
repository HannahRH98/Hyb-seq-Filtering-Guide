# Original script developed by N. Meeprom (https://orcid.org/0000-0003-4193-7062) but further adapted by H. Hall
# This script is designed to analyse HybPiper output statistics ("hybpiper_stats.tsv")
# and assist in deciding filtering strategies based on gene length recovery across samples.
#
# It provides visualisations and numerical summaries showing how many samples
# meet certain gene length recovery thresholds (25%, 50%, 75% - standard hypiper_stats output).
# Users can explore sample retention across a range of gene length recovery cutoffs.
# Can help determine the length threshold of data the user takes onto later downstream analyses.
# This should be used as a guide: removal decisions should also consider sample importance,
# not only gene recovery counts.

library(tidyverse)
library(ggplot2)
library(cowplot)

# Define label for plots
dataset <- element_text("Genes")

# Load HybPiper stats file output; this file should contain per-sample gene recovery stats
stats <- read.table("duplicates_removed_hybpiper_stat.tsv", sep = "\t", header = TRUE)

# Display first few rows to verify data loaded correctly
head(stats)

## Calculate average number of genes recovered at different length recovery thresholds across all samples
gene_recovery_summary <- summarise(stats,
                                   mean_25 = mean(GenesAt25pct),
                                   mean_50 = mean(GenesAt50pct),
                                   mean_75 = mean(GenesAt75pct))

print("Average genes recovered per sample at different sequence recovery thresholds:")
print(gene_recovery_summary)

## Define the range of loci thresholds to examine visually
# This is set to 0-353 for Angiosperms353 bait kit by default; modify if different kit used
genes <- 0:353

# Initialize output dataframe to store counts of samples meeting thresholds
outputs <- data.frame()

# Loop over each gene recovery cutoff threshold to calculate how many samples meet that threshold
for (i in 1:length(genes)){
  count1 <- sum(stats$GenesAt25pct >= genes[i])
  count2 <- sum(stats$GenesAt50pct >= genes[i])
  count3 <- sum(stats$GenesAt75pct >= genes[i])
  outputs[i,1] <- genes[i]
  outputs[i,2] <- (genes[i] * 100/length(genes))
  outputs[i,3] <- count1
  outputs[i,4] <- count2
  outputs[i,5] <- count3
}


colnames(outputs) <- c("Loci", "PercentLoci", "SamplesAt25pct", "SamplesAt50pct", "SamplesAt75pct")

# View the resulting summary table
head(outputs)

# Convert from wide to long format for ggplot2 visualization
outputs_longer <- pivot_longer(outputs,
                               cols = c("SamplesAt25pct", "SamplesAt50pct", "SamplesAt75pct"),
                               names_to = "ThresholdType",
                               values_to = "Samples")

# Plot number of samples retained across loci recovery thresholds, separated by threshold type
samples_retained_plot <- ggplot(outputs_longer, aes(x = Loci, y = Samples, color = ThresholdType)) + 
  geom_line(size = 1) + 
  theme_minimal() +
  labs(title = "Sample Retention Across Sequence Recovery Thresholds",
       subtitle = "Lines show number of samples meeting sequence recovery cutoff",
       x = "Number of Loci Recovered",
       y = "Number of Samples Retained",
       color = "Length Recovery Threshold")

print(samples_retained_plot)

# Save plot as PNG
ggsave("samples_retained_supercontig.png", samples_retained_plot, dpi = 300, width = 8, height = 5)

# You may want to filter based on where the number of recovered loci drops off significantly.
# Different sequence length recovery thresholds can greatly affect how many loci are recovered per sample.
# However, if recovery is generally good, there may be little difference between thresholds,
# and the user can treat them similarly. The key is to balance the trade-off between retaining samples
# recovering more loci, and obtaining longer sequences.

## Calculate slopes (rate of change) of sample retention as loci threshold increases
slope_outputs <- data.frame(Loci = genes, Slope25 = NA, Slope50 = NA, Slope75 = NA)

for (i in 2:length(genes)) {
  slope_outputs$Slope25[i] <- (outputs$SamplesAt25pct[i] - outputs$SamplesAt25pct[i - 1]) / (outputs$Loci[i] - outputs$Loci[i - 1])
  slope_outputs$Slope50[i] <- (outputs$SamplesAt50pct[i] - outputs$SamplesAt50pct[i - 1]) / (outputs$Loci[i] - outputs$Loci[i - 1])
  slope_outputs$Slope75[i] <- (outputs$SamplesAt75pct[i] - outputs$SamplesAt75pct[i - 1]) / (outputs$Loci[i] - outputs$Loci[i - 1])
}

# Reshape for plotting
slope_outputs_longer <- pivot_longer(slope_outputs,
                                    cols = c("Slope25", "Slope50", "Slope75"),
                                    names_to = "ThresholdType",
                                    values_to = "Slope")

# Plot slope of sample retention change across loci length thresholds
slope_plot <- ggplot(slope_outputs_longer, aes(x = Loci, y = Slope, color = ThresholdType)) + 
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Slope of Sample Retention Change Across Loci Thresholds",
       subtitle = "Shows rate at which samples drop out as threshold increases",
       x = "Number of Loci Recovered (Gene Threshold)",
       y = "Slope (Samples Lost per Locus Increase)",
       color = "Recovery Threshold")

print(slope_plot)

# Combine plots vertically
combined_plots <- plot_grid(samples_retained_plot, slope_plot, ncol = 1, align = "v")

# Save combined plot
ggsave("samples_supercontig_slope.png", combined_plots, dpi = 300, width = 8, height = 10)

## Join loci and slope data for export and further analysis
joined_gene_slope <- left_join(outputs, slope_outputs, by = "Loci")

write_csv(joined_gene_slope, "joined_gene_slope.csv")

## Plot histograms of gene recovery counts per sample for additional overview
stats_longer <- pivot_longer(stats, cols = c("GenesWithSeqs", "GenesAt25pct", "GenesAt50pct", "GenesAt75pct"),
                             names_to = "RecoveryThreshold",
                             values_to = "GeneCount")

histogram_plot <- ggplot(stats_longer, aes(x = GeneCount, fill = RecoveryThreshold)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribution of Gene Recovery Across Samples",
       x = "Number of Genes Recovered",
       y = "Number of Samples",
       fill = "Recovery Threshold")

print(histogram_plot)

# Save histogram plot
ggsave("gene_recovery_histogram.png", histogram_plot, dpi = 300, width = 8, height = 5)
