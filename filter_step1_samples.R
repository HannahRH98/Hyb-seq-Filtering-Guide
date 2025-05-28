# Developed by N. Meeprom (https://orcid.org/0000-0003-4193-7062) and H. Hall 
# designed to use standard HybPiper stats output as input

library(tidyverse)
library(ggplot2)
library(cowplot)

# plot a graph to show the samples retained in different thresholds

dataset <- element_text("Genes")

stats <- read.table("hybpiper_stats.tsv", sep = "\t", header = TRUE)
stats

## Count average number of genes recovered

summarise(stats, mean_25 = mean(GenesAt25pct), mean_50 = mean(GenesAt50pct), mean_75 = mean(GenesAt75pct))


## Now, let's use for loop to calculate the number of samples retained from 0-353

genes <- c(0:353)

outputs <- data.frame()

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

colnames(outputs) <- c("Loci", "Percent", "25%","50%", "75%")

outputs

outputs_longer <- pivot_longer(outputs, 
                               cols=c("25%","50%", "75%"), 
                               names_to = "Lengths", values_to = "Samples")
outputs_longer

samples_retained <- ggplot(outputs_longer) + 
  geom_line(aes(x = Loci, y = Samples, col = Lengths)) + 
  theme_minimal() + 
  theme() + 
  labs(title = dataset, 
       x = "Number of Loci", 
       y = "Number of samples",
       col = "Count lengths")
samples_retained

ggsave("samples_retained_supercontig.png",samples_retained, dpi = 300)



## Calculate slope from all observation x and y
## 353 refers to the 353 genes in the angiosperms353 but this can be changed for different target enrichment sets

genes <- c(0:353)

slope_outputs <- data.frame()

for (i in 2:length(genes)){
  slope_25 <- function(x){
    d1 <- outputs[x,1] - outputs[x-1,1]
    d2_25 <- outputs[x,3] - outputs[x-1,3]
    dydx_25 <- d2_25/d1
    return(dydx_25)
  }
  slope_50 <- function(x){
    d1 <- outputs[x,1] - outputs[x-1,1]
    d2_50 <- outputs[x,4] - outputs[x-1,4]
    dydx_50 <- d2_50/d1
    return(dydx_50)
  }
  slope_75 <- function(x){
    d1 <- outputs[x,1] - outputs[x-1,1]
    d2_75 <- outputs[x,5] - outputs[x-1,5]
    dydx_75 <- d2_75/d1
    return(dydx_75)
  }

  slope_outputs[i,1] <- genes[i]
  slope_outputs[i,2] <- slope_25(i)
  slope_outputs[i,3] <- slope_50(i)
  slope_outputs[i,4] <- slope_75(i)
}

colnames(slope_outputs) <- c(dataset, "25%","50%", "75%")

slope_outputs_longer <- pivot_longer(slope_outputs, cols=c("25%","50%", "75%"), 
                                      names_to = "Lengths", values_to = "Slope")
slope_outputs_longer

## Now, plot the graphs from the data above

samples_retained <- ggplot(outputs_longer) + 
  geom_line(aes(x = Supercontigs, y = Samples, col = Lengths)) + 
  theme_minimal() + theme(axis.text.x = element_blank(), 
                          axis.title.x = element_blank(), 
                          legend.position = "none")
samples_retained

slope_sample_retained <- ggplot(slope_outputs_longer) + 
  geom_line(aes(x = Supercontigs, y = Slope, col = Lengths)) + 
  theme_minimal() + theme(legend.position = "bottom")

slope_sample_retained

plot_grid(samples_retained, slope_sample_retained, ncol = 1, align = "v")

ggsave("samples_supercontig_slope.png", dpi = 300)


## Join the data (genes recoverd + slope)

joined_gene_slope <- left_join(outputs, slope_outputs, by = "Genes")
joined_gene_slope

write_csv(joined_gene_slope, "joined_gene_slope.csv")


## histogram

stats_longer <- pivot_longer(stats, cols = c("25%","50%","75%"), names_to = "Count", values_to = "Genes")

ggplot(stats_longer) + geom_histogram(aes(x = Genes, fill = Count), bins = 30, position = "dodge", alpha = 0.5) + theme_minimal()

################################################################################

#Bad samples, using the stat files we imported before

stats_longer <- pivot_longer(stats, cols=c("GenesWithSeqs","GenesAt25pct", "GenesAt50pct", "GenesAt75pct"), 
             names_to = "CountCriteria", values_to = "Count")

ggplot(data = stats_longer) + geom_histogram(aes(x = Count, fill = CountCriteria), bins = 20, position = "dodge")
