### GENERATE BED FILES FOR GENOMIC REGIONS - METHODS ALIGNED
library(dplyr)
library(biomaRt)
library(rtracklayer)
library(GenomicRanges)

## DEFINE TARGET GENES
gene_list <- c("ROBO1", "ROBO2", "ADAMTS19", "SEMA3D", "EPHA7", 
               "SEMA6D", "SFRP1", "WNT2", "SEMA3A", "EPHA3", 
               "ADAMTS9", "VEGFC", "RSPO3", "ID4", "KAT6B", 
               "YAP1", "MSX2", "MEIS2", "ID1", "TBX2", "TBX3")

## CONNECT TO ENSEMBL
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## GET GENE COORDINATES
gene_coords <- getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                                    "start_position", "end_position", "strand"),
                     filters = "hgnc_symbol", 
                     values = gene_list, 
                     mart = human)

## CREATE WINDOWS AS PER METHODS (±10kb not ±100kb)
window_size <- 10000  # 10kb upstream and downstream

gene_windows <- gene_coords %>%
  filter(!chromosome_name %in% c("X", "Y", "MT")) %>%  # Remove non-autosomal
  mutate(
    # Create window coordinates around entire gene body
    window_start = start_position - window_size,
    window_end = end_position + window_size,
    # Format chromosome names for UCSC
    seqnames = paste0("chr", chromosome_name)
  ) %>%
  select(seqnames, window_start, window_end, hgnc_symbol, strand)

## CREATE GRANGES OBJECT
gene_regions <- GRanges(
  seqnames = gene_windows$seqnames,
  ranges = IRanges(start = gene_windows$window_start, 
                   end = gene_windows$window_end),
  strand = ifelse(gene_windows$strand == 1, "+", "-"),
  gene_name = gene_windows$hgnc_symbol
)

## EXPORT BED FILE
export.bed(gene_regions, con = "./data/target_genes_10kb_windows.bed")

## OPTIONAL: CREATE TSS-CENTERED WINDOWS
## If TSS-centered windows are needed instead of gene body
tss_windows <- gene_coords %>%
  filter(!chromosome_name %in% c("X", "Y", "MT")) %>%
  mutate(
    # TSS position based on strand
    tss = ifelse(strand == 1, start_position, end_position),
    # Windows around TSS
    window_start = tss - window_size,
    window_end = tss + window_size,
    seqnames = paste0("chr", chromosome_name)
  ) %>%
  select(seqnames, window_start, window_end, hgnc_symbol, strand)

tss_regions <- GRanges(
  seqnames = tss_windows$seqnames,
  ranges = IRanges(start = tss_windows$window_start, 
                   end = tss_windows$window_end),
  strand = ifelse(tss_windows$strand == 1, "+", "-"),
  gene_name = tss_windows$hgnc_symbol
)

export.bed(tss_regions, con = "./data/target_genes_TSS_10kb_windows.bed")