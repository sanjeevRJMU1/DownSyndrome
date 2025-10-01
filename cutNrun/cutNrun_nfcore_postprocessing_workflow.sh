#!/bin/bash
set -euo pipefail

### 1. AVERAGE BIOLOGICAL REPLICATES 
## Average BigWig files from two biological replicates
BIGWIG_DIR="../../data/01_nfcore_outs/bigwig/"
AVERAGED_DIR="../../data/02_averaged_bigwigs/"
mkdir -p "$AVERAGED_DIR"

## Average replicates using bigWigCompare
bigwigCompare --bigwig1 "$BIGWIG_DIR/PB24_1_R1.bigWig" \
              --bigwig2 "$BIGWIG_DIR/PB24_2_R1.bigWig" \
              --operation mean \
              --outFileName "$AVERAGED_DIR/PB24_averaged.bigWig"

### 2. SUBTRACT IGG BACKGROUND 
## Subtract IgG signal from averaged sample
IGG_SUBTRACTED_DIR="../../data/03_IgG_subtracted/"
mkdir -p "$IGG_SUBTRACTED_DIR"

bigwigCompare --bigwig1 "$AVERAGED_DIR/PB24_averaged.bigWig" \
              --bigwig2 "$BIGWIG_DIR/PB_6_IgG_500k_R1.bigWig" \
              --operation subtract \
              --outFileName "$IGG_SUBTRACTED_DIR/PB24_IgG_subtracted.bigWig"

### 3. CONVERT TO BEDGRAPH FOR PEAK CALLING
## Convert BigWig to BedGraph for SEACR
BEDGRAPH_DIR="../../data/04_bedgraph/"
mkdir -p "$BEDGRAPH_DIR"

for bigwig_file in "$BIGWIG_DIR"/*.bigWig; do
    filename=$(basename "$bigwig_file" .bigWig)
    bedgraph_file="$BEDGRAPH_DIR/${filename}.bedGraph"
    
    echo "Converting $bigwig_file to $bedgraph_file..."
    bigWigToBedGraph "$bigwig_file" "$bedgraph_file"
done

### 4. PEAK CALLING WITH SEACR (FDR 4e-2 as per methods)
## Call peaks using FDR 0.04 with stringent threshold
PEAKS_DIR="../../data/05_SEACR_peaks/"
mkdir -p "$PEAKS_DIR"

for bedgraph_file in "$BEDGRAPH_DIR"/*.bedGraph; do
    base_name=$(basename "$bedgraph_file" .bedGraph)
    output_prefix="$PEAKS_DIR/${base_name}_peaks"
    
    echo "Running SEACR on $bedgraph_file with FDR 0.04..."
    SEACR "$bedgraph_file" 0.04 stringent "$output_prefix"
done

### 5. SUBTRACT IGG PEAKS AND scATAC-seq CONTAMINATING PEAKS
CLEANED_PEAKS_DIR="../../data/06_cleaned_peaks/"
mkdir -p "$CLEANED_PEAKS_DIR"

## Define input files
SAMPLE_PEAKS="$PEAKS_DIR/PB24_1_R1_peaks.stringent.bed"
IGG_PEAKS="$PEAKS_DIR/PB_6_IgG_500k_R1_peaks.stringent.bed"
SCATAC_CONTAM="../../data/scATAC_marker_peaks/fibroblast_endoderm_markers.bed"

## Subtract using bedtools subtract -A 
bedtools subtract -A \
    -a "$SAMPLE_PEAKS" \
    -b "$IGG_PEAKS" | \
bedtools subtract -A \
    -a stdin \
    -b "$SCATAC_CONTAM" \
    > "$CLEANED_PEAKS_DIR/PB24_1_cleaned.bed"

### 6. CREATE CONSENSUS PEAKS FROM REPLICATES
## Process both replicates and create consensus
CONSENSUS_DIR="../../data/07_consensus/"
mkdir -p "$CONSENSUS_DIR"

## Clean both replicates
for rep in 1 2; do
    SAMPLE_PEAKS="$PEAKS_DIR/PB24_${rep}_R1_peaks.stringent.bed"
    
    bedtools subtract -A \
        -a "$SAMPLE_PEAKS" \
        -b "$IGG_PEAKS" | \
    bedtools subtract -A \
        -a stdin \
        -b "$SCATAC_CONTAM" \
        > "$CLEANED_PEAKS_DIR/PB24_${rep}_cleaned.bed"
done

## Create union consensus
cat "$CLEANED_PEAKS_DIR"/PB24_*_cleaned.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - \
    > "$CONSENSUS_DIR/union_consensus.bed"

## Create intersection consensus (peaks in both replicates)
bedtools intersect \
    -a "$CLEANED_PEAKS_DIR/PB24_1_cleaned.bed" \
    -b "$CLEANED_PEAKS_DIR/PB24_2_cleaned.bed" | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - \
    > "$CONSENSUS_DIR/intersection_consensus.bed"

### 7. GENERATE COVERAGE MATRICES FOR TARGET GENES
## Using deepTools
MATRIX_DIR="../../data/08_deeptools/"
mkdir -p "$MATRIX_DIR"

## Compute matrix at gene regions
computeMatrix scale-regions \
    -S "$IGG_SUBTRACTED_DIR/PB24_IgG_subtracted.bigWig" \
    -R ../../data/target_genes_10kb_windows.bed \
    -b 10000 -a 10000 \
    --regionBodyLength 5000 \
    --skipZeros \
    -o "$MATRIX_DIR/target_genes_coverage_matrix.gz"

## Generate summary statistics
multiBigwigSummary BED-file \
    -b "$IGG_SUBTRACTED_DIR/PB24_IgG_subtracted.bigWig" \
    --BED ../../data/target_genes_10kb_windows.bed \
    -o "$MATRIX_DIR/coverage_summary.npz"
