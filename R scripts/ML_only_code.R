# Set working directory
setwd("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/ML_only")

# Load required packages
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)

# Source custom ML functions
source("refer.ML.R")

# Load training dataset
train_data <- read.csv("training_data.csv", row.names = 1, stringsAsFactors = FALSE)

# Separate features and labels
train_features <- train_data[, -ncol(train_data), drop = FALSE]
train_labels <- train_data[, ncol(train_data), drop = F]

# Load testing dataset
test_data <- read.csv("testing_data.csv", row.names = 1, stringsAsFactors = FALSE)
test_features <- test_data[, -ncol(test_data), drop = FALSE]
test_labels <- test_data[, ncol(test_data), drop = F]

# Extract cohort information from sample names
test_labels$Cohort <- gsub("(.+)\\_(.+)\\_(.+)", '\\1', rownames(test_data))

# Ensure training and testing datasets share the same features
common_features <- intersect(colnames(train_features), colnames(test_features))
train_data <- as.matrix(train_data[, common_features])
test_data <- as.matrix(test_data[, common_features])

# Standardize data
train_data = scaleData(train_data, centerFlags = TRUE, scaleFlags = TRUE)
test_data = scaleData(test_data, cohort = test_labels$Cohort, centerFlags = TRUE, scaleFlags = TRUE)

# ==== Diabetes machine learning analysis workflow ====

# Load list of machine learning methods
methods <- read.table("113_ML_methods.txt", header = TRUE, sep = "\t", check.names = FALSE)
methods <- methods$x

# Prepare ML model parameters
classVar = "Group"            # Name of classification variable
min.selected.var = 5          # Threshold for number of selected genes
Variable = colnames(train_features)
preTrain.method = c("Lasso", "glmBoost", "RF", "Stepglm[both]", "Stepglm[backward]")

###################### Build ML models using training data ######################
# Initialize list to store preselected variables
preTrain.var <- list()

# Set seed for reproducibility
set.seed(123)

# Initialize list to store runtime for each method
time.list <- list()

# Loop through each pre-training method to select features
for (method in preTrain.method) {
  time.taken <- system.time({
    preTrain.var[[method]] <- RunML(
      method = method,
      Train_set = train_data,
      Train_label = train_labels,
      mode = "Variable",  # Run in variable selection mode
      classVar = classVar
    )
  })
  
  # Store elapsed time in seconds
  time.list[[method]] <- time.taken[["elapsed"]]
  
  # Print runtime for current method
  cat(sprintf("Method [%s] took %.2f seconds\n", method, time.taken[["elapsed"]]))
}

# Add all features as a simple baseline
preTrain.var[["simple"]] <- colnames(train_data)

# Load UpSetR package for visualizing feature intersections
library(UpSetR)

# Extract selected gene lists from each algorithm
gene_lists <- list(
  Lasso = preTrain.var$Lasso,
  glmBoost = preTrain.var$glmBoost,
  RF = preTrain.var$RF,
  Step_both = preTrain.var$`Stepglm[both]`,
  Step_back = preTrain.var$`Stepglm[backward]`
)

# Create UpSet matrix
upset_data <- fromList(gene_lists)
allcolour = c('#DF0A1F','#1C5BA7','#019E73','#ED621B','#E477C1')

# Plot UpSet diagram
pdf("gene_intersection_upset.pdf", width = 10, height = 6)
upset(upset_data, sets = names(gene_lists),
      order.by = "freq",
      text.scale = 1.2,
      matrix.color = allcolour,
      mainbar.y.label = "Gene number intersected",
      sets.x.label = "Gene number selected")
dev.off()

# Identify core intersecting genes
core_genes <- Reduce(intersect, gene_lists)
write.table(core_genes, "core_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

###################### Build ML models using selected features ######################
model <- list()           # Initialize list to store models
set.seed(123)             # Set seed for reproducibility
Train_set_bk = train_data # Backup training data

# Loop through each ML method to build models
for (method in methods) {
  cat(match(method, methods), ":", method, "\n")
  method_name = method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)
  
  Variable = preTrain.var[[method[1]]]
  train_data = Train_set_bk[, Variable]
  Train_label = train_labels
  
  # Build ML model
  model[[method_name]] <- RunML(
    method = method[2],        # ML algorithm
    Train_set = train_data,    # Training expression data
    Train_label = train_labels,# Training labels
    mode = "Model",            # Run in model building mode
    classVar = classVar
  )
  
  # Remove model if number of selected variables is below threshold
  if (length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}

# Restore original training data
train_data = Train_set_bk; rm(Train_set_bk)

# Save all ML model results
saveRDS(model, "model.MLmodel.rds")

###################### Calculate sample prediction scores ######################
# Load ML model results
model <- readRDS("model.MLmodel.rds")

# Identify valid models based on selected features
methodsValid <- names(model)

# Predict risk scores for each sample using valid models
RS_list <- list()
for (method in methodsValid) {
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], new_data = rbind.data.frame(train_data, test_data))
}

# Combine risk scores into a matrix
riskTab = as.data.frame(t(do.call(rbind, RS_list)))
riskTab = cbind(id = row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Predict class labels for each sample
Class_list <- list()
for (method in methodsValid) {
  Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = rbind.data.frame(train_data, test_data))
}

Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
classTab = cbind(id = row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Extract variables selected by each ML model
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

fea_df <- lapply(model, function(fit) {
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file = "model.genes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###################### Calculate AUC for each model ######################
AUC_list <- list()
for (method in methodsValid) {
  AUC_list[[method]] <- RunEval(
    fit = model[[method]],        # ML model
    Test_set = test_data,         # Testing expression data
    Test_label = test_labels,     # Testing labels
    Train_set = train_data,       # Training expression data
    Train_label = train_labels,   # Training labels
    Train_name = "Train",         # Training set name
    cohortVar = "Cohort",         # Cohort ID
    classVar = classVar           # Classification variable
  )
}

# Combine AUC results into a matrix
AUC_mat <- do.call(rbind, AUC_list)
aucTab = cbind(Method = row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

###################### Plot AUC heatmap ######################
# Load AUC matrix
AUC_mat <- read.table("model.AUCmatrix.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)

# Sort ML models by mean AUC
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = TRUE)
AUC_mat <- AUC_mat[names(avg_AUC),]

# Select features from the best model
fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

# Define colors for cohort annotations
CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
names(CohortCol) <- colnames(AUC_mat)

# Plot heatmap
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(
  Cindex_mat = AUC_mat,       # AUC matrix
  avg_Cindex = avg_AUC,       # Average AUC
  CohortCol = CohortCol,      # Cohort colors
  barCol = "steelblue",       # Color of side bar
  cellwidth = cellwidth, cellheight = cellheight,
  cluster_columns = FALSE, cluster_rows = FALSE
)

pdf(file = "ML_AUC_heatmap.pdf", width = cellwidth * ncol(AUC_mat) + 6, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

###################### ROC analysis per cohort ######################
library(pROC)

rsFile = "model.riskMatrix.txt"      # Risk score matrix file
method = "glmBoost+LDA"              # Selected ML method (based on heatmap)

# Load risk score file
riskRT = read.table(rsFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Extract cohort ID from sample names
CohortID = gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID = gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort = CohortID

# Loop through cohorts to plot ROC curves
for (Cohort in unique(riskRT$Cohort)) {
  rt = riskRT[riskRT$Cohort == Cohort, ]
  
  # Extract sample class (Control vs Case)
  y = gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
  y = ifelse(y == "Control", 0, 1)
  
  # Generate ROC curve
  roc1 = roc(y, as.numeric(rt[, method]))
  ci1 = ci.auc(roc1, method = "bootstrap")
  ciVec = as.numeric(ci1)
  
  # Save ROC plot
  pdf(file = paste0("ROC.", Cohort, ".pdf"), width = 5, height = 4.75)
  plot(roc1, print.auc = TRUE, col = "red", legacy.axes = TRUE, main = Cohort)
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "red")
  dev.off()
}

###################### ROC analysis for individual genes ######################
library(glmnet)
library(pROC)

expFile = "merged_normalized_data.txt"  # Expression data file
geneFile = "model.genes.txt"            # Selected genes file

# Load expression data
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Extract sample class (Control vs Case)
y = gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(rt))
y = ifelse(y == "Control", 0, 1)

# Load gene list
geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
geneRT <- geneRT$V1[geneRT$V2 == 'glmBoost+RF']

# Define color palette for plotting
bioCol = rainbow(length(geneRT), s = 0.9, v = 0.9)

# Loop through genes to plot ROC curves
aucText = c()
k = 0
for (x in as.vector(geneRT)) {
  k = k + 1
  roc1 = roc(y, as.numeric(rt[x, ]))
  
  if (k == 1) {
    pdf(file = "ROC.genes.pdf", width = 9, height = 9)
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "", lwd = 3)
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  } else {
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "", lwd = 3, add = TRUE)
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  }
}

# Add legend showing AUC for each gene
legend("bottomright", aucText, lwd = 3, bty = "n", cex = 0.8, col = bioCol[1:(ncol(rt)-1)], inset = c(0.05, 0))
dev.off()
