library(DESeq2)
library(ggplot2)
library(pheatmap)

# Read count matrix
count_data <- read.csv("count1hs2.csv", row.names = 1)

# Read metadata (sample information)
metadata <- read.csv("metadata1hs1.csv", row.names = 1)


# Ensure row names match column names in count matrix
all(rownames(metadata) %in% colnames(count_data))  # Should return TRUE
all(colnames(count_data) %in% rownames(metadata))  # Should return TRUE

dds <- DESeqDataSetFromMatrix(countData = count_data, 
                              colData = metadata, 
                              design = ~ condition)  # Modify design formula based on your experiment


colnames(count_data) <- gsub("\\.", "-", colnames(count_data))


# Remove lowly expressed genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Normalize counts
dds <- DESeq(dds)

# Regularized log transformation (for visualization)
rld <- rlog(dds, blind = FALSE)

# Define reference level if needed
dds$condition <- relevel(dds$condition, ref = "control")  # Adjust reference level

# Run differential expression
res <- results(dds, name = "condition_infected_vs_control")
# Adjust groups

# Order by p-value
res <- res[order(res$padj), ]

# Save results
write.csv(as.data.frame(res), "DEG_results_vishwas1hs2.csv")

res$logP <- -log10(res$pvalue)
ggplot(res, aes(x = log2FoldChange, y = logP)) +
  geom_point(aes(color = padj < 0.05)) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 p-value")

plotMA(res, main = "DESeq2 MA Plot", ylim = c(-5, 5))

top_genes <- head(order(res$padj, na.last = NA), 50)  # Select top 50 DEGs
pheatmap(assay(rld)[top_genes, ], cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, annotation_col = metadata)

# Select significant genes (FDR < 0.05 & |log2FC| > 1)
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Separate Upregulated Genes (log2FC > 1)
upregulated <- subset(sig_genes, log2FoldChange > 1)

# Separate Downregulated Genes (log2FC < -1)
downregulated <- subset(sig_genes, log2FoldChange < -1)

# Save results as CSV
write.csv(as.data.frame(sig_genes), "Significant_DEGs_vishwas1hs2.csv")         # All DEGs
write.csv(as.data.frame(upregulated), "Upregulated_Genes_vishwas1hs2.csv")      # Upregulated
write.csv(as.data.frame(downregulated), "Downregulated_Genes_vishwas1hs2.csv")  # Downregulated

# Sort DEGs by padj (smallest first)
top_genes <- head(sig_genes[order(sig_genes$padj), ], 10)
print(top_genes)

top_up <- head(upregulated[order(-upregulated$log2FoldChange), ], 10)
print(top_up)

top_down <- head(downregulated[order(downregulated$log2FoldChange), ], 10)
print(top_down)


