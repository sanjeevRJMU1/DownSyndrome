library(Seurat)
library(dplyr)
library(SeuratWrappers)
library(biomaRt)
set.seed(123)
options(future.globals.maxSize = 90 * 1024^3)

### 1. DATA IMPORT
## Read 10X Cell Ranger aggregated output
data_10x <- Read10X(data.dir = "filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = data_10x, 
                                 project = "project_name",
                                 min.cells = 3, 
                                 min.features = 200)

## Add metadata from Cell Ranger aggregation sheet
aggr_sheet <- read.csv("aggregation.csv")
aggr_sheet$gem.group <- aggr_sheet$sample_id
aggr_sheet$cellID <- seq.int(nrow(aggr_sheet))

cellID <- as.numeric(gsub(".*-", "", colnames(seurat_obj)))
metadata_df <- data.frame(cellID = cellID, row.names = colnames(seurat_obj))
metadata_df <- merge(metadata_df, aggr_sheet, by = "cellID", all.x = TRUE)
rownames(metadata_df) <- colnames(seurat_obj)[match(metadata_df$cellID, cellID)]
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_df[,-1])

### 2. QC METRICS
## Human data patterns
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

## Mouse data patterns (use when applicable)
# seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
# seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")

## Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## Convert to mouse genes if needed
convertHumanGeneList <- function(x){
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                   host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", 
                   host = "https://dec2021.archive.ensembl.org/")
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                    values = x, mart = human, attributesL = c("mgi_symbol"), 
                    martL = mouse, uniqueRows=T)
  return(unique(genesV2[, 2]))
}

## For mouse: 
# s.genes <- convertHumanGeneList(s.genes)
# g2m.genes <- convertHumanGeneList(g2m.genes)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- CellCycleScoring(seurat_obj, 
                               s.features = s.genes,
                               g2m.features = g2m.genes)

### 3. QC FILTERING
## Standard thresholds from methods (adjust as needed for specific experiments)
seurat_obj <- subset(seurat_obj,
                    subset = nFeature_RNA > 2000 & 
                            nFeature_RNA < 8000 &
                            nCount_RNA < 30000 &
                            percent.mt < 15)

### 4. SCT NORMALIZATION WITH CELL CYCLE REGRESSION
seurat_obj <- SCTransform(seurat_obj, 
                         vars.to.regress = c("S.Score", "G2M.Score"),
                         return.only.var.genes = TRUE,
                         conserve.memory = FALSE)

### 5. PCA ANALYSIS
seurat_obj <- RunPCA(seurat_obj, npcs = 35)

### 6. BATCH CORRECTION USING FASTMNN
seurat_list <- SplitObject(seurat_obj, split.by = "gem.group")
seurat_obj <- RunFastMNN(object.list = seurat_list)

### 7. CLUSTERING
seurat_obj <- RunUMAP(seurat_obj, reduction = "mnn", dims = 1:35)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "mnn", dims = 1:35)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.25)

### 8. MARKER IDENTIFICATION
markers <- FindAllMarkers(seurat_obj)

### 9. OPTIONAL: FILTER X/Y CHROMOSOME GENES
## For sex-chromosome independent analysis
library(org.Hs.eg.db)  # Use org.Mm.eg.db for mouse
gene_ids <- rownames(seurat_obj)
chr_map <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "CHR",
                  keytype = "SYMBOL", multiVals = "first")
blacklist_xy <- names(chr_map)[chr_map %in% c("X", "Y")]
## Apply filter after differential expression analysis

### 10. CARDIAC CELL SELECTION (multiple valid approaches)
## Approach A: Marker-based selection
cardiac_subset <- subset(seurat_obj, TNNT2 > 0)

## Approach B: Cluster-based selection
# cardiac_subset <- subset(seurat_obj, idents = c("4", "5", "11"))

### 11. RE-PROCESS CARDIAC SUBSET
cardiac_subset <- SCTransform(cardiac_subset,
                             vars.to.regress = c("S.Score", "G2M.Score"),
                             return.only.var.genes = TRUE,
                             conserve.memory = FALSE)
cardiac_subset <- RunPCA(cardiac_subset, npcs = 35)
cardiac_list <- SplitObject(cardiac_subset, split.by = "gem.group")
cardiac_subset <- RunFastMNN(object.list = cardiac_list)
cardiac_subset <- RunUMAP(cardiac_subset, reduction = "mnn", dims = 1:35)
cardiac_subset <- FindNeighbors(cardiac_subset, reduction = "mnn", dims = 1:35)
cardiac_subset <- FindClusters(cardiac_subset, resolution = 0.25)

### 12. CELL TYPE ANNOTATION
## Manual annotation based on marker expression
cardiac_subset$cell_type <- "Unknown"
## Example assignments based on marker genes:
## AVC markers: TBX2, RSPO3, BMP2
## Ventricular markers: MYL2, IRX4, MYH7
# cardiac_subset$cell_type[cardiac_subset$seurat_clusters == "1"] <- "AVC"
# cardiac_subset$cell_type[cardiac_subset$seurat_clusters == "0"] <- "Ventricle"

### 13. SUBSET TO SPECIFIC POPULATIONS
Idents(cardiac_subset) <- "cell_type"
avc_subset <- subset(cardiac_subset, idents = "AVC")

## Optional: Further expression-based refinement
# RSPO3_subset <- subset(avc_subset, RSPO3 > 0.01)
# TBX2_subset <- subset(avc_subset, TBX2 > 0.01)

### 14. DIFFERENTIAL EXPRESSION ANALYSIS
Idents(avc_subset) <- "condition"
de_results <- FindMarkers(avc_subset,
                         ident.1 = "trisomic",
                         ident.2 = "disomic",
                         logfc.threshold = 0.25,
                         min.pct = 0.1)

## Optional: Filter X/Y genes from results
if(exists("blacklist_xy")) {
  de_results <- de_results[!(rownames(de_results) %in% blacklist_xy),]
}

### 15. MODULE SCORING
## AVC module
avc_genes <- c("TBX2", "RSPO3", "HAND2", "TGFB2", "SMAD7", "MSX2", "GATA6", "WNT11")
cardiac_subset <- AddModuleScore(cardiac_subset, 
                                features = list(avc_genes),
                                name = "AVCmodule")

## Ventricular module
vent_genes <- c("MYL2", "MYL4", "FHL2", "MYH7", "RBM20", "HOPX", "TBX20", "TBX5")
cardiac_subset <- AddModuleScore(cardiac_subset,
                                features = list(vent_genes),
                                name = "VentModule")

### 16. ITERATIVE QC (as per methods)
## "Iterative rounds of filtering poor quality clusters and re-running clustering"
## Identify and remove poor quality clusters based on:
## - High mitochondrial content
## - Low gene counts
## - Cell cycle markers
## Then repeat from step 4

### 17. DATA OUTPUT
saveRDS(seurat_obj, "seurat_full_dataset.RDS")
saveRDS(cardiac_subset, "cardiac_subset.RDS")
saveRDS(avc_subset, "avc_subset.RDS")
write.csv(de_results, "differential_expression_results.csv")

## Export cell proportions
metadata <- cardiac_subset@meta.data
cell_proportions <- metadata %>%
  group_by(condition, cell_type) %>%
  summarise(n = n()) %>%
  group_by(condition) %>%
  mutate(percent = n / sum(n) * 100)
write.csv(cell_proportions, "cell_proportions.csv")

################################################################################
### CROP-seq SPECIFIC ADDITIONS
################################################################################

## Process amplicon library to assign sgRNAs
## Output from get_barcodes.py script
sgRNA_mapping <- read.csv("sgRNA_cell_mapping.csv")
seurat_obj$sgRNA <- sgRNA_mapping$guide[match(colnames(seurat_obj), 
                                               sgRNA_mapping$cell_barcode)]
seurat_obj$target_gene <- sgRNA_mapping$target[match(colnames(seurat_obj),
                                                      sgRNA_mapping$cell_barcode)]

## Filter for single-guide cells
seurat_obj_cropseq <- subset(seurat_obj, !is.na(sgRNA))

## Export log-normalized expression for classifier
DefaultAssay(seurat_obj_cropseq) <- "RNA"
seurat_obj_cropseq <- NormalizeData(seurat_obj_cropseq)
expression_matrix <- as.matrix(GetAssayData(seurat_obj_cropseq, slot = "data"))
write.csv(expression_matrix, "log_normalized_expression_cropseq.csv")

## Calculate Gini coefficient for guide distribution
library(DescTools)
guide_counts <- table(seurat_obj_cropseq$sgRNA)
gini_coefficient <- Gini(as.numeric(guide_counts))