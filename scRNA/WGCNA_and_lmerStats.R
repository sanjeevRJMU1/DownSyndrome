### WGCNA ANALYSIS PIPELINE FOR scRNA-seq DATA
library(Seurat)
library(WGCNA)
library(enrichR)
library(lmerTest)
library(ggplot2)
set.seed(123)
options(stringsAsFactors = FALSE)


### STEP 1: DATA PREPARATION
## Load Seurat object (should be already processed/filtered)
seurat_obj <- readRDS("path/to/seurat_object.RDS")

## Extract expression data
## Option A: Use variable features only
var_features <- VariableFeatures(seurat_obj)
# Or use scaled data features
# var_features <- rownames(seurat_obj@assays$SCT@scale.data)

## Option B: Use all genes
# var_features <- rownames(seurat_obj)

## Prepare expression matrix (cells as rows, genes as columns)
expr_matrix <- t(as.matrix(GetAssayData(seurat_obj, slot = "data")))
expr_matrix_subset <- expr_matrix[, var_features]

print(paste("Matrix dimensions:", nrow(expr_matrix_subset), "cells x", 
            ncol(expr_matrix_subset), "genes"))


### STEP 2: NETWORK CONSTRUCTION
## Enable parallel processing
allowWGCNAThreads(nThreads = 4)

## Construct network
net <- blockwiseModules(
  expr_matrix_subset,
  power = 10,                    # Soft threshold (typically 6-12 for scRNA)
  corType = "bicor",            # Robust correlation
  networkType = "signed",       # Signed network preserves direction
  minModuleSize = 10,           # Minimum module size
  reassignThreshold = 0,
  mergeCutHeight = 0.15,        # Module merging threshold
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM",
  verbose = 3
)

## Save network object
saveRDS(net, "WGCNA_network.RDS")


### STEP 3: MODULE VISUALIZATION
## Plot dendrogram with module colors
plotDendroAndColors(net$dendrograms[[1]], 
                    net$colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)

## Module summary
print(table(net$colors))
modules <- unique(net$colors)
print(paste("Total modules:", length(modules)))

## Save gene-to-module mapping
gene_module_df <- data.frame(
  gene = names(net$colors),
  module = net$colors
)
write.csv(gene_module_df, "gene_module_mapping.csv")


### STEP 4: MODULE SCORING IN SEURAT
## Add module scores to Seurat object
for(module_color in modules[modules != "grey"]) {
  
  ## Get genes in module
  module_genes <- names(net$colors)[net$colors == module_color]
  
  ## Add module score
  seurat_obj <- AddModuleScore(
    seurat_obj,
    features = list(module_genes),
    name = module_color
  )
}

## Save updated object
saveRDS(seurat_obj, "seurat_with_module_scores.RDS")


### STEP 5: MODULE VISUALIZATION
## Create output directory
dir.create("WGCNA_results", showWarnings = FALSE)

## Visualize module scores
metadata <- seurat_obj@meta.data
module_scores <- paste0(modules[modules != "grey"], "1")

## Boxplots by condition
for(module_score in module_scores) {
  p <- ggplot(metadata, aes(x = condition, y = .data[[module_score]], 
                            fill = condition)) +
    geom_boxplot() +
    labs(title = paste("Module:", gsub("1$", "", module_score)),
         y = "Module Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0("WGCNA_results/", module_score, "_boxplot.pdf"), 
         p, width = 6, height = 5)
}

## Feature plots on UMAP
for(module_score in module_scores[1:min(6, length(module_scores))]) {
  p <- FeaturePlot(seurat_obj, features = module_score)
  ggsave(paste0("WGCNA_results/", module_score, "_umap.pdf"), 
         p, width = 6, height = 6)
}


### STEP 6: FUNCTIONAL ENRICHMENT
## Run GO enrichment for each module
for(module_color in modules[modules != "grey"]) {
  
  ## Create output directory
  dir.create(paste0("WGCNA_results/", module_color), showWarnings = FALSE)
  
  ## Get module genes
  module_genes <- names(net$colors)[net$colors == module_color]
  
  ## Run enrichR
  enrich_results <- enrichr(module_genes, 
                            databases = c("GO_Biological_Process_2021"))
  
  ## Filter significant results
  go_results <- enrich_results$GO_Biological_Process_2021
  go_sig <- go_results[go_results$Adjusted.P.value < 0.05, ]
  
  ## Save results
  write.csv(go_sig, 
           paste0("WGCNA_results/", module_color, "/GO_enrichment.csv"))
  
  ## Plot top terms
  if(nrow(go_sig) > 0) {
    p <- plotEnrich(go_sig[1:min(20, nrow(go_sig)), ], 
                   showTerms = 20, 
                   numChar = 40, 
                   y = "Count", 
                   orderBy = "P.value")
    ggsave(paste0("WGCNA_results/", module_color, "/GO_plot.pdf"), 
           p, width = 10, height = 8)
  }
  
  ## Save module genes
  write.csv(data.frame(gene = module_genes),
           paste0("WGCNA_results/", module_color, "/genes.csv"))
}


### STEP 7: STATISTICAL TESTING (OPTIONAL)
## Test module differences between conditions using mixed effects model
## This is useful if you have batch effects or paired samples

test_module_differences <- function(metadata, module_score, 
                                   test_var = "condition", 
                                   random_var = "batch") {
  
  ## Prepare data
  metadata[[module_score]] <- as.numeric(metadata[[module_score]])
  
  ## Run mixed effects model if random effect exists
  if(!is.null(random_var) && random_var %in% colnames(metadata)) {
    formula_str <- paste0(module_score, " ~ ", test_var, " + (1|", random_var, ")")
    model <- lmer(as.formula(formula_str), data = metadata)
  } else {
    ## Simple linear model
    formula_str <- paste0(module_score, " ~ ", test_var)
    model <- lm(as.formula(formula_str), data = metadata)
  }
  
  return(summary(model))
}

## Test all modules
if("condition" %in% colnames(metadata)) {
  
  results_list <- list()
  
  for(module_score in module_scores) {
    test_result <- test_module_differences(metadata, module_score)
    results_list[[module_score]] <- test_result
    
    ## Extract p-value (adjust based on your model structure)
    if(class(test_result)[1] == "summary.merMod") {
      p_val <- coef(test_result)[2, "Pr(>|t|)"]
    } else {
      p_val <- coef(test_result)[2, "Pr(>|t|)"]
    }
    
    print(paste(module_score, "p-value:", p_val))
  }
  
  ## Multiple testing correction
  p_values <- sapply(module_scores, function(x) {
    if(class(results_list[[x]])[1] == "summary.merMod") {
      coef(results_list[[x]])[2, "Pr(>|t|)"]
    } else {
      coef(results_list[[x]])[2, "Pr(>|t|)"]
    }
  })
  
  p_adjusted <- p.adjust(p_values, method = "fdr")
  
  ## Save results
  stats_results <- data.frame(
    module = gsub("1$", "", module_scores),
    p_value = p_values,
    p_adjusted = p_adjusted
  )
  write.csv(stats_results, "WGCNA_results/module_statistics.csv")
}


### STEP 8: MODULE EIGENGENE ANALYSIS (OPTIONAL)


## Calculate module eigengenes
ME_data <- moduleEigengenes(expr_matrix_subset, net$colors)
module_eigengenes <- ME_data$eigengenes

## Correlate with metadata variables
if("condition" %in% colnames(metadata)) {
  
  ## Convert condition to numeric if categorical
  if(is.factor(metadata$condition) || is.character(metadata$condition)) {
    cond_numeric <- as.numeric(factor(metadata$condition))
  } else {
    cond_numeric <- metadata$condition
  }
  
  ## Calculate correlations
  correlations <- cor(module_eigengenes, cond_numeric, use = "pairwise.complete.obs")
  
  ## Plot heatmap
  pdf("WGCNA_results/module_condition_correlation.pdf", width = 8, height = 6)
  heatmap(as.matrix(correlations), 
          Colv = NA, 
          scale = "none",
          col = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()
}

print("WGCNA analysis complete")
print(paste("Results saved in:", "WGCNA_results/"))