### GSE138458

library(Biobase)
library(GEOquery)
library(limma)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(hgu133plus2.db)

## Expression data preprocessing

geo_file <- "GSE61635_series_matrix.txt.gz"

# Read and skip annotation lines starting with "!"
geo_data <- read.delim(geo_file, comment.char = "!", header = TRUE, stringsAsFactors = FALSE)

head(geo_data[, 1:5])
dim(geo_data)

## Probe annotation

gpl_table <- read.delim2("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/数据集/GSE61635/GPL570-55999.txt", row.names=NULL, comment.char="#")

# Extract probe ID and gene symbol
anno <- gpl_table[, c("ID", "Gene.Symbol")]
colnames(anno) <- c("ID_REF", "gene_symbol")

anno$gene_symbol <- sub(" ///.*", "", anno$gene_symbol)

rownames(geo_data) <- geo_data$ID_REF
rownames(anno) <- anno$ID_REF

exp <- cbind(anno, geo_data)

exp <- exp[!duplicated(exp$gene_symbol),]
exp <- na.omit(exp)
exp <- exp[rownames(exp) != "1552829_at", ]
rownames(exp) <- exp$gene_symbol
exp <- exp[,-c(1:3)]
exp <- as.data.frame(t(exp))

## Clinical information

group_list <- factor(c(rep("SLE", 99), rep("HC", 30)))

names(group_list) <- rownames(exp)

exp$type <- group_list

exp <- exp[, c("type", setdiff(colnames(exp), "type"))]

exp$type <- as.character(exp$type)

GSE61635 <- exp

SLE_cohorts_diffexp[["GSE61635"]] <- GSE61635

### GSE49454

library(Biobase)
library(GEOquery)
library(limma)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(hgu133plus2.db)

## Expression data preprocessing

geo_file <- "GSE49454_series_matrix.txt.gz"

# Read and skip annotation lines starting with "!"
geo_data <- read.delim(geo_file, comment.char = "!", header = TRUE, stringsAsFactors = FALSE)

head(geo_data[, 1:5])
dim(geo_data)

## Probe annotation

gpl_table <- read.delim2("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/数据集/GSE49454/GPL10558-50081.txt", row.names=NULL, comment.char="#")

colnames(gpl_table)

# Extract probe ID and gene symbol
anno <- gpl_table[, c("ID", "ILMN_Gene")]
colnames(anno) <- c("ID_REF", "gene_symbol")

anno$gene_symbol <- sub(" ///.*", "", anno$gene_symbol)

anno <- anno[anno$gene_symbol != "" & !is.na(anno$gene_symbol), ]

rownames(geo_data) <- geo_data$ID_REF
rownames(anno) <- anno$ID_REF

exp <- cbind(anno, geo_data)

exp <- exp[!duplicated(exp$gene_symbol),]
exp <- na.omit(exp)
# exp <- exp[rownames(exp) != "1552829_at", ]
rownames(exp) <- exp$gene_symbol
exp <- exp[,-c(1:3)]
exp <- as.data.frame(t(exp))

## Clinical information

# Read metadata lines (skip expression matrix)
meta_lines <- readLines("GSE49454_series_matrix.txt.gz")

# Extract lines containing sample characteristics
meta_info <- meta_lines[grepl("^!Sample_characteristics_ch1", meta_lines)]

# Extract all fields containing characteristics
meta <- grep("^!Sample_characteristics_ch1", meta_lines, value = TRUE)

# Remove prefix
meta_clean <- sub("^!Sample_characteristics_ch1\\t", "", meta)

# Split into columns by tab
meta_df <- as.data.frame(do.call(rbind, strsplit(meta_clean, "\t")))

group_list <- factor(c(rep("SLE", 157), rep("HC", 20)))

names(group_list) <- rownames(exp)

exp$type <- group_list

exp <- exp[, c("type", setdiff(colnames(exp), "type"))]

exp$type <- as.character(exp$type)

GSE49454 <- exp

SLE_cohorts_diffexp[["GSE49454"]] <- GSE49454

### GSE112087

library(Biobase)
library(GEOquery)
library(limma)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(hgu133plus2.db)

## Expression data preprocessing

GSE112087_count <- read.delim("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/数据集/GSE112087/GSE112087_counts-matrix-EnsembIDs-GRCh37.p10.txt")

## Probe annotation

library(org.Hs.eg.db)
library(AnnotationDbi)

# Remove version number from Ensembl IDs (e.g., ENSG00000123456.1 → ENSG00000123456)
ids <- gsub("\\..*", "", GSE112087_count$X)

# Map to gene symbol
gene_symbol <- mapIds(org.Hs.eg.db, keys = ids,
                      column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Merge
GSE112087_count$symbol <- gene_symbol
GSE112087_count <- GSE112087_count[!is.na(GSE112087_count$symbol), ]
GSE112087_count <- GSE112087_count[!duplicated(GSE112087_count$symbol), ]
rownames(GSE112087_count) <- GSE112087_count$symbol
GSE112087_count <- GSE112087_count[, -ncol(GSE112087_count)]
GSE112087_count <- GSE112087_count[,-1]

# Convert count data using DESeq2

library(DESeq2)

# Construct DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = GSE112087_count,
                              colData = data.frame(row.names = colnames(GSE112087_count)),
                              design = ~ 1)  # No grouping, only transformation

# VST transformation
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Extract normalized expression matrix
exp <- as.data.frame(assay(vsd))

exp <- as.data.frame(t(exp))

## Clinical information

group_list <- factor(c(rep("SLE", 62), rep("HC", 58)))

names(group_list) <- rownames(exp)

exp$type <- group_list

exp <- exp[, c("type", setdiff(colnames(exp), "type"))]

exp$type <- as.character(exp$type)

GSE112087 <- exp

SLE_cohorts_diffexp[["GSE112087"]] <- GSE112087

### GSE72509

library(Biobase)
library(GEOquery)
library(limma)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(hgu133plus2.db)

## Expression data preprocessing

GSE72509_SLE_RPKMs <- read.csv("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/数据集/GSE72509/GSE72509_SLE_RPKMs.csv")

GSE72509_SLE_RPKMs <- GSE72509_SLE_RPKMs[!duplicated(GSE72509_SLE_RPKMs$SYMBOL),]

rownames(GSE72509_SLE_RPKMs) <- GSE72509_SLE_RPKMs$SYMBOL

GSE72509_SLE_RPKMs <- GSE72509_SLE_RPKMs[,-1]

exp <- as.data.frame(t(GSE72509_SLE_RPKMs))

## Clinical information

exp$type <- ifelse(grepl("SLE", rownames(exp)), "SLE",
                  ifelse(grepl("control", rownames(exp)), "HC", NA))

table(exp$type)

# Move the type column to the first position
exp <- exp[, c(ncol(exp), 1:(ncol(exp) - 1))]

# Check results
head(exp)[1:5,1:5]

GSE72509 <- exp

SLE_cohorts_diffexp[["GSE72509"]] <- GSE72509

### GSE169080

GSE169080_SLE <- read.delim("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/数据集/GSE169080/GSE169080_all.counts.SLE.txt", row.names=1)

GSE169080_control <- read.delim("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/数据集/GSE169080/GSE169080_control_all_counts.txt")
GSE169080_control <- GSE169080_control[!duplicated(GSE169080_control$GeneID),]
rownames(GSE169080_control) <- GSE169080_control$GeneID
GSE169080_control <- GSE169080_control[,-1]

commonid <- intersect(rownames(GSE169080_SLE),rownames(GSE169080_control))

GSE169080_control <- GSE169080_control[commonid,]
GSE169080_SLE <- GSE169080_SLE[commonid,]

exp <- cbind(GSE169080_SLE, GSE169080_control)

# Convert count data using DESeq2

library(DESeq2)

# Construct DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = data.frame(row.names = colnames(exp)),
                              design = ~ 1)  # No grouping, only transformation

# VST transformation
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Extract normalized expression matrix
exp <- as.data.frame(assay(vsd))

exp <- as.data.frame(t(exp))

## Clinical information

group_list <- factor(c(rep("SLE", 4), rep("HC", 3)))

names(group_list) <- rownames(exp)

exp$type <- group_list

exp <- exp[, c("type", setdiff(colnames(exp), "type"))]

exp$type <- as.character(exp$type)

GSE169080 <- exp

SLE_cohorts_diffexp[["GSE169080"]] <- GSE169080

### GSE144390

library(Biobase)
library(GEOquery)
library(limma)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(hgu133plus2.db)

## Expression data preprocessing

geo_file <- "GSE144390_series_matrix.txt"

# Read and skip annotation lines starting with "!"
geo_data <- read.delim(geo_file, comment.char = "!", header = TRUE, stringsAsFactors = FALSE)

head(geo_data[, 1:5])
dim(geo_data)

## Probe annotation

gpl_table <- read.delim("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/数据集/GSE144390/GPL6244-17930.txt", comment.char="#")

# Extract probe ID and gene symbol
anno <- gpl_table[, c("ID", "gene_assignment")]
colnames(anno) <- c("ID_REF", "gene_symbol")

anno <- anno[anno$gene_symbol != "---", ]

anno$gene_symbol <- sapply(strsplit(anno$gene_symbol, " // "), `[`, 2)

rownames(geo_data) <- geo_data$ID_REF
rownames(anno) <- anno$ID_REF

comonid <- intersect(rownames(geo_data),rownames(anno))

geo_data <- geo_data[comonid,]

exp <- cbind(anno, geo_data)

exp <- exp[!duplicated(exp$gene_symbol),]
exp <- na.omit(exp)
# exp <- exp[rownames(exp) != "1552829_at", ]
rownames(exp) <- exp$gene_symbol
exp <- exp[,-c(1:3)]
exp <- as.data.frame(t(exp))

## Clinical information

group_list <- factor(c(rep("HC", 3), rep("SLE", 3)))

names(group_list) <- rownames(exp)

exp$type <- group_list

exp <- exp[, c("type", setdiff(colnames(exp), "type"))]

exp$type <- as.character(exp$type)

GSE144390 <- exp

SLE_cohorts_diffexp[["GSE144390"]] <- GSE144390

#### Save and check
saveRDS(SLE_cohorts_diffexp, file = "SLE_cohorts_diffexp_10.rds")
