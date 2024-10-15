## Setup
### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)

## Load in data
data <- read.table("Data/count.txt", header=T, row.names=1)
meta <- read.table("Data/meta.txt", header=T, row.names=1)

### Check classes of the data we just brought in
class(meta)  # Check the class of the metadata
class(data)  # Check the class of the count data

# Create DESeqDataSet from matrix input
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

## Run analysis
dds <- DESeq(dds)

# Uncomment these lines if you need to install these packages
#BiocManager::install("apeglm")
#install.packages('ashr')

#### Define contrasts, extract results table, and shrink the log2 fold changes
contrast_oe <- c("sampletype", "SARS", "MOCK")
res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)
res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken, type="ashr")

# Plot MA (mean-difference) plots
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
plotMA(res_tableOE, ylim=c(-2,2))

# Check class and metadata columns of the results table
class(res_tableOE)
mcols(res_tableOE, use.names=T)

# Convert results table to data frame for viewing
res_tableOE %>% data.frame() %>% View()

# Summarize the results table
summary(res_tableOE)

# Define significance cutoffs
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

# Convert results table to tibble and filter for significant results
res_tableOE_tb <- res_tableOE %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE <- res_tableOE_tb %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sigOE

# Write significant results to DESeq2.csv
write.csv(sigOE, "Output/DESeq2.csv")

# Load additional libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)

# Create a column to indicate significant genes
res_tableOE_tb <- res_tableOE_tb %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(res_tableOE_tb) + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) + 
  ggtitle("Overexpression") + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))

## Create a column to indicate which genes to label
res_tableOE_tb <- res_tableOE_tb %>% 
  arrange(padj) %>% 
  mutate(genelabels = "")

res_tableOE_tb$genelabels[1:10] <- res_tableOE_tb$gene[1:10]
View(res_tableOE_tb)

# Updated volcano plot with gene labels
ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(colour = threshold_OE)) + 
  geom_text_repel(aes(label = genelabels)) + 
  ggtitle("Overexpression") + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))

# Display session information
sessionInfo()
