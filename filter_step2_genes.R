# Developed by H. Hall & N. Meeprom (https://orcid.org/0000-0003-4193-7062)

library(tidyverse)
library(ggplot2)
library(cowplot)

citation("ape")

seq_lengths <- read.table("seq_lengths_genes_manual.tsv", sep = "\t", header = TRUE)

# Make sure that the first column and row that contain our sequence length values are correct

seq_lengths[2,2]

# Count seqeunces with 0 bp for each gene (empty sequence)

num.rows <- nrow(seq_lengths)
num.cols <- ncol(seq_lengths)

# designate the first column that contain sequence lengths (usually 2, but depending on your table format)
from <- 2

count.empty <- data.frame(colnames(seq_lengths[2:num.cols]))
zero <- data.frame()

for (x in from:num.cols){
  zero[x, 1] <- sum(seq_lengths[2:num.rows, x] == 0)
}

count.empty[,2] <- data.frame(zero[2:num.cols,1])
count.empty[,3] <- zero[2:num.cols,1]*100/(num.rows-1)

colnames(count.empty) <- c("Gene", "Empty", "PercentEmpty")

View(count.empty)

# Plot a histogram showing how many empty sequences are found in each gene

ggplot(data = count.empty) + geom_histogram(aes(x = Empty), bins = 50)

# Or plot a line graph to see it more clearly

ggplot(data = count.empty) + 
  geom_bar(aes(x = reorder(Gene, -Empty), y = PercentEmpty), stat = "identity") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Genes", y = "Percent of Empty Sequences")

## Print genes that contain certain amount/percentage of empty sequences
## Look into the plots above and count.empty file to see if there is any clear-cut percentage of empty genes

threshold.zero <-80

empty_select <- filter(count.empty, PercentEmpty >= threshold.zero)

print(empty_select[1], row.names = FALSE)

sink("gene_empty_75.txt")
print(empty_select[1], row.names = FALSE)
sink()


# Then, open the output text file and manually remove spaces. This file is now ready to be used as a genelist file to remove undesired genes from the dataset before running trees or can run trees for all of them and remove these gene trees later.

## Now, see how many genes would are there if we set a criteria of how much empty sequence we would accept for those genes

sum(count.empty[1:nrow(count.empty), 2] >= threshold.zero)

accept <- data.frame()

for (i in 1:353){
  accept[i ,1] <- i
  accept[i, 2] <- sum(count.empty[1:nrow(count.empty), 3] >= i)
}

colnames(accept) <- c("ThreholdEmptyAccepted", "GenesToRemove")

head(accept)

ggplot(data = accept) +
  geom_line(aes(x = ThreholdEmptyAccepted, y = GenesToRemove))

# histogram for frequency of genes to remove over the thresholds

ggplot(data = accept) + geom_histogram(aes(x = GenesToRemove), bins = 100)


## Count NOT EMPTY

num.rows <- nrow(seq_lengths)
num.cols <- ncol(seq_lengths)

count.present <- data.frame(colnames(seq_lengths[2:num.cols]))
present <- data.frame()

for (x in from:num.cols){
  present[x, 1] <- sum(seq_lengths[2:num.rows, x] > 0)
}

count.present[,2] <- data.frame(present[2:num.cols,1])
count.present[,3] <- present[2:num.cols,1]*100/(num.rows-1)

colnames(count.present) <- c("Gene", "Presnt", "PercentPresent")

count.present
