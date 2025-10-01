### STATISTICAL TESTING FOR MODULE SCORES IN scRNA-seq
library(Seurat)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(dplyr)
set.seed(123)


### OPTION 1: LINEAR MIXED EFFECTS (LMER)
## Use when you have batch effects or paired samples
run_lmer_analysis <- function(seurat_obj, 
                              module_names,
                              test_var = "condition",
                              random_var = "batch",
                              control_level = NULL,
                              output_dir = "lmer_results") {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metadata <- seurat_obj@meta.data
  
  ## Run LMER for each module
  results_list <- list()
  
  for(module in module_names) {
    
    ## Build formula
    formula_str <- paste0(module, " ~ ", test_var, " + (1|", random_var, ")")
    
    ## Fit model
    model <- lmer(as.formula(formula_str), data = metadata)
    
    ## Extract coefficients
    coef_summary <- summary(model)$coefficients
    
    ## Get p-values for each comparison vs control
    if(!is.null(control_level)) {
      comparisons <- rownames(coef_summary)[grepl(test_var, rownames(coef_summary))]
      p_values <- coef_summary[comparisons, "Pr(>|t|)"]
      
      ## Adjust p-values
      p_adjusted <- p.adjust(p_values, method = "fdr")
      
      ## Create results table
      results <- data.frame(
        module = module,
        comparison = gsub(test_var, "", comparisons),
        p_value = p_values,
        p_adjusted = p_adjusted,
        stringsAsFactors = FALSE
      )
      
      results_list[[module]] <- results
      write.csv(results, file.path(output_dir, paste0(module, "_lmer_results.csv")), 
                row.names = FALSE)
    }
  }
  
  return(results_list)
}

## Visualization with LMER p-values
plot_lmer_results <- function(seurat_obj,
                              module_name,
                              test_var = "condition",
                              lmer_results,
                              output_file = NULL) {
  
  metadata <- seurat_obj@meta.data
  
  ## Prepare data for plotting
  plot_data <- data.frame(
    condition = metadata[[test_var]],
    value = metadata[[module_name]]
  )
  
  ## Create boxplot
  p <- ggplot(plot_data, aes(x = condition, y = value)) +
    geom_boxplot(aes(fill = condition), outlier.shape = NA) +
    labs(y = paste0(module_name, " score"),
         title = paste0("Module: ", module_name)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ## Add p-value annotations if available
  if(!is.null(lmer_results) && module_name %in% names(lmer_results)) {
    module_stats <- lmer_results[[module_name]]
    
    ## Create annotation dataframe
    y_max <- max(plot_data$value, na.rm = TRUE)
    module_stats$y_position <- y_max + seq_len(nrow(module_stats)) * 0.1 * diff(range(plot_data$value))
    module_stats$group1 <- levels(plot_data$condition)[1]
    module_stats$group2 <- module_stats$comparison
    
    p <- p + stat_pvalue_manual(module_stats,
                                label = "p_adjusted",
                                tip.length = 0,
                                step.increase = 0.05)
  }
  
  if(!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = 6)
  }
  
  return(p)
}


### OPTION 2: WILCOXON TEST
## Use for simple pairwise comparisons without batch effects
run_wilcox_analysis <- function(seurat_obj,
                                module_names,
                                group_var = "condition",
                                comparisons = "auto",
                                p_adjust_method = "BH",
                                output_dir = "wilcox_results") {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metadata <- seurat_obj@meta.data
  results_list <- list()
  
  ## Auto-generate comparisons if needed
  if(comparisons == "auto") {
    groups <- levels(factor(metadata[[group_var]]))
    comparisons <- list()
    for(i in 2:length(groups)) {
      comparisons[[i-1]] <- c(groups[1], groups[i])
    }
  }
  
  for(module in module_names) {
    
    ## Prepare data
    test_data <- data.frame(
      group = metadata[[group_var]],
      value = metadata[[module]]
    )
    
    ## Run comparisons
    stats_results <- compare_means(
      value ~ group,
      data = test_data,
      comparisons = comparisons,
      method = "wilcox.test",
      p.adjust.method = p_adjust_method
    )
    
    results_list[[module]] <- stats_results
    write.csv(stats_results, 
             file.path(output_dir, paste0(module, "_wilcox_results.csv")),
             row.names = FALSE)
  }
  
  return(results_list)
}

## Visualization with Wilcoxon p-values
plot_wilcox_results <- function(seurat_obj,
                               module_name,
                               group_var = "condition",
                               comparisons,
                               p_adjust_method = "BH",
                               output_file = NULL) {
  
  metadata <- seurat_obj@meta.data
  
  ## Prepare data
  plot_data <- data.frame(
    group = metadata[[group_var]],
    value = metadata[[module_name]]
  )
  
  ## Create plot with statistical comparisons
  p <- ggplot(plot_data, aes(x = group, y = value)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA) +
    stat_compare_means(comparisons = comparisons,
                      method = "wilcox.test",
                      p.adjust.method = p_adjust_method,
                      label = "p.signif",
                      tip.length = 0) +
    labs(y = paste0(module_name, " score"),
         title = paste0("Module: ", module_name)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  if(!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = 6)
  }
  
  return(p)
}


### USAGE EXAMPLE


## Load data and add module scores
seurat_obj <- readRDS("path/to/seurat_object.RDS")

## Define gene modules
avc_genes <- c("TBX2", "RSPO3", "BMP2", "MSX2", "HAND2")
vent_genes <- c("MYL2", "MYL4", "FHL2", "MYH7", "TBX20")

## Add module scores
seurat_obj <- AddModuleScore(seurat_obj, list(avc_genes), name = "AVC_module")
seurat_obj <- AddModuleScore(seurat_obj, list(vent_genes), name = "Vent_module")

module_names <- c("AVC_module1", "Vent_module1")

## METHOD 1: LMER (if you have batch effects)
if("batch" %in% colnames(seurat_obj@meta.data)) {
  
  lmer_results <- run_lmer_analysis(
    seurat_obj,
    module_names = module_names,
    test_var = "condition",
    random_var = "batch",
    control_level = "Control"
  )
  
  ## Plot results
  for(module in module_names) {
    plot_lmer_results(
      seurat_obj,
      module_name = module,
      test_var = "condition",
      lmer_results = lmer_results,
      output_file = paste0("plots/", module, "_lmer.pdf")
    )
  }
}

## METHOD 2: Wilcoxon (simple comparisons)
comparisons <- list(
  c("Control", "Treatment1"),
  c("Control", "Treatment2")
)

wilcox_results <- run_wilcox_analysis(
  seurat_obj,
  module_names = module_names,
  group_var = "condition",
  comparisons = comparisons
)

## Plot results
for(module in module_names) {
  plot_wilcox_results(
    seurat_obj,
    module_name = module,
    group_var = "condition",
    comparisons = comparisons,
    output_file = paste0("plots/", module, "_wilcox.pdf")
  )
}

print("Statistical analysis complete")