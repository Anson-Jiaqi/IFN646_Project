library(edgeR)
library(ggplot2)
library(ggrepel)

# Load the counts data from a file
Counts <- read.table("Data/count.txt", header=TRUE, row.names=1)
dim(Counts)   # Check the dimensions of the data
head(Counts)  # Display the first few rows of the counts data

# Create a DGEList object from the counts data
dgList <- DGEList(counts=Counts, genes=rownames(Counts))
dgList
dgList$samples
head(dgList$counts)  # Display the first few rows of the counts in DGEList
head(dgList$genes)   # Display the first few genes in DGEList

# Calculate counts per million (CPM)
countsPerMillion <- cpm(dgList)
summary(countsPerMillion)  # Summary of CPM values

# Filter out low expression genes
countCheck <- countsPerMillion > 1
head(countCheck)  # Display the first few rows of the countCheck matrix
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
summary(cpm(dgList))  # Summary of CPM values after filtering

# Calculate normalization factors using TMM (trimmed mean of M-values)
dgList <- calcNormFactors(dgList, method="TMM")
plotMDS(dgList)  # Multi-dimensional scaling plot

# Define sample types
sampleType <- rep("N", ncol(dgList))  # N=normal; T=tumor
sampleType[sample(grep("SARS", colnames(dgList)))] <- "S"
sampleType

# Create sample replicate identifiers
sampleReplicate <- paste("S", rep(1:2, each=2), sep="")
sampleReplicate

# Create design matrix for the GLM
designMat <- model.matrix(~sampleReplicate + sampleType)
designMat

# Estimate dispersions for the GLM
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
plotBCV(dgList)  # Biological coefficient of variation plot

# Fit the GLM and perform likelihood ratio test
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=3)
edgeR_result <- topTags(lrt, Inf)
edgeR_result

# Filter significant results
sig_infected <- edgeR_result$table[edgeR_result$table$FDR < 0.05 & abs(edgeR_result$table$logFC) > 0.58, ]
sig_infected

# Save significant results to a CSV file
write.csv(sig_infected, "Output/edgeR.csv")

# Add a significance column for plotting
edgeR_result$table$Significant <- "Not Significant"
edgeR_result$table$Significant[edgeR_result$table$FDR < 0.05 & abs(edgeR_result$table$logFC) > 0.58] <- "Significant"

# Plot the Volcano plot
ggplot(edgeR_result$table, aes(x=logFC, y=-log10(FDR), color=Significant)) + 
  geom_point(alpha=0.4, size=1.75) + 
  scale_color_manual(values=c("grey", "red")) + 
  theme_minimal() + 
  labs(title="Volcano Plot", x="Log Fold Change", y="-Log10 FDR") + 
  theme(legend.position="right") + 
  geom_text_repel(aes(label=ifelse(Significant == "Significant", rownames(edgeR_result$table), "")), 
                  size=2.5, max.overlaps=10)
