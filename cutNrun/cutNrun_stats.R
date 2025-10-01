library(reticulate)
library(multcompView)

## Load deeptools .npz
np  <- reticulate::import("numpy")
arr <- np$load("../../data/08_deeptools/coverage_summary.npz", allow_pickle = TRUE)
mat <- as.matrix(arr$f[["matrix"]])

## Get sample labels if available
keys <- names(arr$f)
labs <- if ("labels" %in% keys) as.character(arr$f[["labels"]]) else
  if ("bw_names" %in% keys) as.character(arr$f[["bw_names"]]) else
    paste0("sample_", seq_len(ncol(mat)))
colnames(mat) <- labs

## Long format for ANOVA
df <- stack(as.data.frame(mat))
colnames(df) <- c("coverage", "sample")
df <- df[is.finite(df$coverage), ]
df$sample <- factor(df$sample, levels = unique(labs))

## One-way ANOVA and Tukey HSD
fit <- aov(coverage ~ sample, data = df)
tk  <- TukeyHSD(fit)

## Compact letter display
cld <- multcompLetters4(fit, tk)
letters_df <- data.frame(
  sample = names(cld$sample$Letters),
  cld    = as.character(cld$sample$Letters),
  row.names = NULL
)

