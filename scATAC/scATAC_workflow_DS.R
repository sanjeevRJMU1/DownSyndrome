### COMPREHENSIVE scATAC-seq ANALYSIS PIPELINE WITH ArchR
library(ArchR)
library(dplyr)
library(Seurat)
library(SeuratWrappers)
set.seed(123)
addArchRThreads(threads = 4)


### STEP 1: SETUP AND QUALITY CONTROL
## 1.1 Initialize ArchR
addArchRGenome("hg38")

## 1.2 Define input files and samples
input_files <- c("path/to/sample1/fragments.tsv.gz",
                 "path/to/sample2/fragments.tsv.gz",
                 "path/to/sample3/fragments.tsv.gz",
                 "path/to/sample4/fragments.tsv.gz")

sample_names <- c("Control_1", "Treatment_1", 
                  "Control_2", "Treatment_2")

## 1.3 Create Arrow Files
arrow_files <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = sample_names,
  minTSS = 6,
  minFrags = 10000,
  maxFrags = 3e+06,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

## 1.4 Infer Doublets
doublet_scores <- addDoubletScores(
  input = arrow_files,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  force = TRUE
)

## 1.5 Create ArchR Project
proj <- ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = "output/scATAC_analysis",
  copyArrows = TRUE
)

print(proj)

## 1.6 Quality Control Plots
p1 <- plotGroups(proj, groupBy = "Sample", name = "TSSEnrichment", plotAs = "ridges")
p2 <- plotGroups(proj, groupBy = "Sample", name = "log10(nFrags)", plotAs = "violin")
p3 <- plotFragmentSizes(proj)
p4 <- plotTSSEnrichment(proj)

plotPDF(p1, p2, p3, p4, name = "QC-Metrics.pdf", ArchRProj = proj)

## 1.7 Filter Doublets
proj <- filterDoublets(proj)
saveArchRProject(proj, outputDirectory = "output/scATAC_analysis")


### STEP 2: DIMENSION REDUCTION AND CLUSTERING
## 2.1 Iterative LSI
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(resolution = c(0.2), n.start = 10),
  varFeatures = 25000,
  dimsToUse = 1:25,
  LSIMethod = 2,
  force = TRUE
)

## 2.2 Batch Correction with Harmony (optional)
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:25,
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

## 2.3 Clustering
proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",  # or "IterativeLSI" if no batch correction
  dimsToUse = 1:25,
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6,
  force = TRUE
)

## 2.4 UMAP Embedding
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "UMAP",
  metric = "cosine",
  force = TRUE
)

## 2.5 Visualize Clusters
p5 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p6 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

plotPDF(p5, p6, name = "UMAP-Clusters.pdf", ArchRProj = proj)
saveArchRProject(proj)


### STEP 3A: GENE ACTIVITY SCORES
## 3A.1 Identify Marker Genes
markers_gs <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 1000
)

marker_list <- getMarkers(markers_gs, cutOff = "FDR <= 0.05 & Log2FC >= 1")

## 3A.2 Visualize Marker Genes
marker_genes <- c("TNNT2", "EPCAM", "TBX2", "RSPO3", "IRX4", 
                  "MYH7", "WT1", "COL1A1", "ISL1", "TAGLN")

## Add imputation weights for better visualization
proj <- addImputeWeights(proj, dimsToUse = 1:25)

p7 <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = marker_genes,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p7, name = "GeneScores-Markers.pdf", ArchRProj = proj)


### STEP 3B: scRNA-seq INTEGRATION
## 3B.1 Load scRNA-seq Reference
seurat_rna <- readRDS("path/to/scRNA_reference.RDS")
seurat_rna <- FindVariableFeatures(seurat_rna, nfeatures = 3000)
se_rna <- as.SingleCellExperiment(seurat_rna, assay = "SCT")

## 3B.2 Integrate with scATAC
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = se_rna,
  addToArrow = TRUE,
  groupRNA = "cell_type",  # Column in scRNA metadata
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  nGenes = 3000,
  useImputation = FALSE,
  force = TRUE
)

## 3B.3 Visualize Integration
p8 <- plotEmbedding(proj, colorBy = "cellColData", 
                    name = "predictedGroup", embedding = "UMAP")
plotPDF(p8, name = "RNA-Integration.pdf", ArchRProj = proj)

## 3B.4 Create Manual Annotations Based on Integration
## Add custom cluster labels based on predicted cell types
cluster_labels <- c("C1" = "CM-Ventricle", "C2" = "CM-AVC", 
                    "C3" = "Fibroblast", "C4" = "Endoderm")
proj$Clusters_Annotated <- cluster_labels[proj$Clusters]

saveArchRProject(proj)


### STEP 4: PEAK CALLING
## 4.1 Create Pseudo-bulk Replicates
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Clusters_Annotated",  # Or use treatment groups
  minCells = 40,
  maxCells = 5000,
  maxFragments = 50 * 10^6,
  force = TRUE
)

## 4.2 Call Peaks
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters_Annotated",
  pathToMacs2 = "/path/to/macs2",
  peaksPerCell = 2000,
  maxPeaks = 200000,
  minCells = 25,
  reproducibility = "2",
  force = TRUE
)

## 4.3 Add Peak Matrix
proj <- addPeakMatrix(proj, binarize = FALSE, ceiling = 4, force = TRUE)

## 4.4 Identify Marker Peaks
marker_peaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters_Annotated",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

## 4.5 Export Marker Peaks
marker_list <- getMarkers(marker_peaks, cutOff = "FDR <= 0.05 & Log2FC >= 1.5")

## Export peaks for each cluster
for(cluster in names(marker_list)) {
  peaks_df <- as.data.frame(marker_list[[cluster]])
  peaks_bed <- peaks_df[order(peaks_df$Log2FC, decreasing = TRUE), ]
  peaks_bed$rank <- seq_len(nrow(peaks_bed))
  
  write.table(peaks_bed[, c("seqnames", "start", "end", "rank", "Log2FC")],
              file = paste0("output/peaks/", cluster, "_markers.bed"),
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}

## 4.6 Motif Enrichment
proj <- addMotifAnnotations(proj, motifSet = "homer", name = "Motif", force = TRUE)

motif_enrichment <- peakAnnoEnrichment(
  seMarker = marker_peaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

heatmap_motifs <- plotEnrichHeatmap(motif_enrichment, n = 5, transpose = TRUE)
plotPDF(heatmap_motifs, name = "Motif-Enrichment.pdf", ArchRProj = proj)

saveArchRProject(proj)


### STEP 5: DIFFERENTIAL ACCESSIBILITY
## 5.1 Create Comparison Groups
## Add metadata for comparisons (e.g., cluster + treatment)
proj$Comparison_Group <- paste0(proj$Clusters_Annotated, "_", proj$Sample)

## 5.2 Run Differential Analysis Per Cluster
clusters <- unique(proj$Clusters_Annotated)
dar_results <- list()

for(cluster in clusters) {
  control_group <- paste0(cluster, "_Control")
  treatment_group <- paste0(cluster, "_Treatment")
  
  ## Check if both groups exist with sufficient cells
  if(sum(proj$Comparison_Group == control_group) >= 25 & 
     sum(proj$Comparison_Group == treatment_group) >= 25) {
    
    dar <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = "Comparison_Group",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = control_group,
      bgdGroups = treatment_group
    )
    
    dar_results[[cluster]] <- dar
    
    ## Plot MA
    p_ma <- plotMarkers(dar, name = control_group, 
                       cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
    plotPDF(p_ma, name = paste0("MA_", cluster, ".pdf"), ArchRProj = proj)
  }
}

## 5.3 Export DAR Results
for(cluster in names(dar_results)) {
  dar_gr <- getMarkers(dar_results[[cluster]], 
                      cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", 
                      returnGR = TRUE)
  
  if(length(dar_gr) > 0) {
    dar_df <- as.data.frame(dar_gr)
    write.csv(dar_df, paste0("output/DAR/", cluster, "_DAR.csv"))
  }
}

## 5.4 Create BigWigs for Visualization
getGroupBW(
  ArchRProj = proj,
  groupBy = "Comparison_Group",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000
)
saveArchRProject(proj)