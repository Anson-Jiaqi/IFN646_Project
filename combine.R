# Read edgeR and DESeq2 results from CSV files
edgeR_data <- read.csv("Output/edgeR.csv")
DESeq2_data <- read.csv("Output/DESeq2.csv")

# Print the column names of each data frame to ensure they are read correctly
print(names(edgeR_data))
print(names(DESeq2_data))

# Extract the gene columns from each data frame
edgeR_genes <- edgeR_data$genes
DESeq2_genes <- DESeq2_data$gene

# Find the common genes between the two data sets
common_elements <- intersect(edgeR_genes, DESeq2_genes)
print(common_elements)

# Filter the edgeR data to include only the common genes
filtered_genes_ed <- edgeR_data[edgeR_data$genes %in% common_elements, ]

# Filter the DESeq2 data to include only the common genes
filtered_genes_de <- DESeq2_data[DESeq2_data$gene %in% common_elements, ]

# Print the filtered data frames to check the results
print(filtered_genes_ed)
print(filtered_genes_de)

# Write the filtered data frames to new CSV files
write.csv(filtered_genes_ed, "Output/edgeR_combined.csv")
write.csv(filtered_genes_de, "Output/DESeq2_combined.csv")
