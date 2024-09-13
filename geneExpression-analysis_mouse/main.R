# Load the dataset
mouse_data <- readRDS("~/Documents/KutumLab/Summer 2024/BroadSeq Project/geneExpression_mouse/mouse_salmon_E-MTAB-6798.rds")
str(mouse_data)
head(mouse_data)
summary(mouse_data) #Output: "DESeqDataSet object of length 35655 with 6 metadata columns"
rownames(mouse_data)
colnames(mouse_data)
# View the metadata for a specific gene
rowData(mouse_data)["ENSMUSG00000102041", ]

#statistics of the column data
coldata <- colData(mouse_data)
head(coldata)

#statistics of the row data
rowdata <- rowData(mouse_data)
summary(rowdata)
head(rowdata)

#the count data is the matrix which contains raw gene expression counts for each gene across various samples.
#row contains the gene Id while the column contains the sample Id. 
# The value in the matrix represent the raw counts of gene expression. O means the gene was not expressed in that sample
counts <- assay(mouse_data)
head(counts)

# Bar plot of tissue types
ggplot(as.data.frame(coldata), aes(x = factor_organism_part)) +
  geom_bar() +
  theme_minimal() +
  ggtitle("Distribution of Tissue Types")
# Bar plot of developmental stages
ggplot(as.data.frame(coldata), aes(x = factor_dev_stage)) +
  geom_bar() +
  theme_minimal() +
  ggtitle("Distribution of Developmental Stages")

# Check available assays
assayNames(mouse_data)
# Log-transform for normalization : Converts raw counts to log2 counts, which stabilizes variance
logcounts <- log2(assay(mouse_data, "counts") + 1)
assays(mouse_data)$logcounts <- logcounts
# Visualize the distribution of log-transformed counts
hist(logcounts, breaks = 50, main = "Distribution of Log-Transformed Counts", xlab = "Log2(Counts + 1)", col = "green")


# Perform PCA using logcounts
computedPCA_logcounts <- broadSeq::prcompTidy(mouse_data, scaledAssay = "logcounts", ntop = 500)

# Define custom color palette and shape scale
color_palette <- scales::hue_pal()(length(unique(colData(mouse_data)$factor_organism_part)))
shape_scale <- seq(0, length(unique(colData(mouse_data)$factor_dev_stage)) - 1)

# Plot PCA logcounts
pca_logcounts_plot <- plotAnyPC(
  computedPCA = computedPCA_logcounts, 
  x = 1, y = 2, 
  color = "factor_organism_part", shape = "factor_dev_stage", 
  legend = "bottom"
) + 
  scale_color_manual(values = color_palette) + 
  scale_shape_manual(values = shape_scale)

# Display the plot
print(pca_logcounts_plot)


#Plot PCA logcounts with size as factor_dev_stage.
pca_logcounts_plot2 <- plotAnyPC(
  computedPCA = computedPCA_logcounts, 
  x = 1, y = 2, 
  color = "factor_organism_part", size = "factor_dev_stage", 
  legend = "bottom"
) + 
  scale_color_manual(values = color_palette) + 
  scale_shape_manual(values = shape_scale)

# Display the plot
print(pca_logcounts_plot2)

