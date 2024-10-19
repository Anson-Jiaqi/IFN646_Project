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
class(meta)
class(data)

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

dds <- estimateSizeFactors(dds)                              
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
dds
### Plot PCA 
plotPCA(rld, intgroup="sampletype")


# Input is a matrix of log transformed values
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100



# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
#ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sampletype))
ggplot(df, aes(x = PC3, y = PC4, color = sampletype)) +
  geom_point() +
  labs(title = "PCA Plot", 
       x = paste0("PC3: ", round(percentVar[3], 2), "% variance"),
       y = paste0("PC4: ", round(percentVar[4], 2), "% variance")) +
  theme_minimal()

### Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames
### Plot heatmap
pheatmap(rld_cor)


dds <- DESeq(dds)
plotDispEsts(dds)
