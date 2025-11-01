
# Perform differential expression analysis between SLE and HC groups
diff <- seurat_diff2(datafilt = CD14Mono,
                     group.by = "disease",
                     group1 = "SLE",
                     group2 = "HC",
                     assay = "RNA",
                     min.pct = 0.1,
                     thres.fc = 0)

# Transpose expression matrix
CD14Mono_exp <- as.data.frame(t(CD14Mono_exp))

# Selected features from ML model
fea_sel

# Prepare input for ECDF plot
input <- data.frame(type = CD14Mono@meta.data$disease,
                    value = CD14Mono_exp$IFIT2)

# Plot ECDF for a specific gene
ecdf_plot(input = input,
          ylim = c(0,1),
          title = "IFIT2")

# Define gene to plot
i = "IFIT2"

# Define custom color palette for FeaturePlot
my_colors <- c("#F0F0F0",'#EDD1D8', '#f4a3a8', '#e38191', '#cc607d', '#ad466c', '#8b3058', '#672044')

# Plot gene expression on UMAP
FeaturePlot(GSE135779_combined, features = i, reduction = "umap_harmony", pt.size = 1, raster = TRUE) + 
  theme_bw() +  # Use a clean white background
  theme(
    panel.grid = element_blank(),            # Remove grid lines
    axis.ticks = element_blank(),            # Remove axis ticks
    axis.text = element_blank(),             # Remove axis text
    axis.title = element_text(colour = "black", size = 15),  # Customize axis title
    plot.title = element_text(size = 17, hjust = 0.5),       # Center-align plot title
    panel.border = element_blank(),          # Remove panel border
    legend.position = "right"                # Keep color legend on the right
  ) + 
  scale_color_gradientn(colors = my_colors) +  # Apply custom color gradient
  labs(x = ' ', y = ' ', title = i)

# Loop through selected features to generate UMAP plots
select <- fea_sel
my_colors <- c("#F0F0F0",'#EDD1D8', '#f4a3a8', '#e38191', '#cc607d', '#ad466c', '#8b3058', '#672044')

for (i in select) {
  # Set working directory for saving figures
  setwd("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/113机器学习方案/10个umap")
  
  # Generate UMAP FeaturePlot for each gene
  fig2 = FeaturePlot(GSE135779_combined, features = i, reduction = "umap_harmony",
                     pt.size = 1, raster = TRUE) + 
    theme_bw() +  # Clean white background
    theme(
      panel.grid = element_blank(),            # Remove grid lines
      axis.ticks = element_blank(),            # Remove axis ticks
      axis.text = element_blank(),             # Remove axis text
      axis.title = element_text(colour = "black", size = 15),  # Customize axis title
      plot.title = element_text(size = 17, hjust = 0.5),       # Center-align title
      panel.border = element_blank(),          # Remove panel border
      legend.position = "right"                # Keep legend on the right
    ) + 
    scale_color_gradientn(colors = my_colors) +  # Apply custom color gradient
    labs(x = ' ', y = ' ', title = i)  # Set gene as plot title
  
  # Construct output file name
  file_name <- paste0(i, "umap_harmony.png")
  
  # Save figure as PNG file
  ggsave(file_name, fig2, height = 3.5, width = 4.7)
}

# Filter differential genes based on selected features
diff <- diff[diff$gene %in% select,]

## Begin further analysis

## Filter for protein-coding mRNA only
human_mrna <- read.table("/home/ug1268/tools/human_sc/mRNA_list.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)

# Retain common protein-coding genes
common <- intersect(rownames(CD14Mono), rownames(human_mrna))
CD14Mono_mRNA <- CD14Mono[common,]

# Set cell identities based on Seurat clusters
Idents(CD14Mono_mRNA) <- CD14Mono_mRNA$seurat_clusters
