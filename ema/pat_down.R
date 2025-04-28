library(dplyr)
#library(highcharter)
library(tidyverse)
library(Seurat)
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host="https://asia.ensembl.org"))
print("biomart complete")

args <- commandArgs(TRUE)

# Access specific arguments by index
emadir <- args[1]
runid <- args[2]
#Get GTF
gtf_file <- args[3] # or gtf_file <- "~/Documents/data/PASseqData/humanSTARindex/Homo_sapiens.GRCh38.99.sorted.gtf"  # E.g., downloaded from Ensembl or Gencode

gene = 0# turn this to 1 if need gene exp seurat.

setwd(emadir)

expression_matrix <- ReadMtx(mtx = "filterdmatrix.mtx", # should be names as "matrix.mtx"
                             cells = "filterdcb.tsv", # should be names as "barcodes.tsv"
                             features = "pas_gene.tsv") # should be names as "features.tsv"
print("matrix read")
ann <- read.delim("annotatedpas.bed", header=F)
colnames(ann) <- c("chr", "start", "end", "PASid", "ensembl_gene_id", "hgnc_symbol", "strand")

ann$PASlabel <- paste("chr", ann$chr, ann$start, ann$end, ann$ensembl_gene_id, ann$hgnc_symbol,  sep="_")
rownames(expression_matrix) <- ann$PASlabel

print("annotation completed")
##################################################################
#Make a gene expression table
# gene_exp <- as.data.frame(expression_matrix)
# gene_exp$hgnc_symbol <- sapply(strsplit(as.character(rownames(gene_exp)), "_"), `[`, 5)
# gene_exp_merged <- aggregate(. ~ hgnc_symbol, gene_exp, sum)
# rownames(gene_exp_merged) <- gene_exp_merged$hgnc_symbol
# gene_exp_merged <- gene_exp_merged[, -c(1)]
##################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(gene == 1){
  
  library(data.table)
  library(stringi)
  
  print("1. Convert the matrix to data.table (directly from the matrix for efficiency)")
  #gene_exp <- as.data.table(expression_matrix, keep.rownames = TRUE)
  # Initialize an empty data.table for the result
  gene_exp <- data.table()
  # Process the matrix row-by-row in chunks
  chunk_size <- 5000  # Number of rows to process at a time
  total_rows <- nrow(expression_matrix)
  for (start_row in seq(1, total_rows, by = chunk_size)) {
    end_row <- min(start_row + chunk_size - 1, total_rows)
    # Subset the rows to work on
    temp_chunk <- as.data.table(expression_matrix[start_row:end_row, , drop = FALSE], keep.rownames = TRUE)
    # Extract 'hgnc_symbol' in the chunk
    temp_chunk[, hgnc_symbol := tstrsplit(rn, "_", fixed = TRUE)[[5]]]
    # Bind the processed chunk to the result
    gene_exp <- rbindlist(list(gene_exp, temp_chunk), use.names = TRUE)
    # Clear temporary objects and trigger garbage collection
    rm(temp_chunk)
    gc()
  }
  
  print("2. Extract hgnc_symbol efficiently using stringi")
  gene_exp[, hgnc_symbol := stri_split_fixed(rn, "_", simplify = TRUE)[, 5]]
  
  print("3. Process the aggregation in chunks to handle memory")
  chunk_size <- 2000  # Adjust based on your system's capacity
  unique_symbols <- unique(gene_exp$hgnc_symbol)
  num_chunks <- ceiling(length(unique_symbols) / chunk_size)
  
  result_list <- vector("list", num_chunks)
  
  for (i in seq_len(num_chunks)) {
    print(i)
    # Identify the chunk
    symbol_chunk <- unique_symbols[((i - 1) * chunk_size + 1):(i * chunk_size)]
    symbol_chunk <- na.omit(symbol_chunk)  # Ensure no NAs
    # Subset and aggregate only for the chunk
    result_list[[i]] <- gene_exp[hgnc_symbol %in% symbol_chunk, 
                                 lapply(.SD, sum), 
                                 by = hgnc_symbol, .SDcols = -c("rn")]
    # Remove intermediate objects and run garbage collection
    gc()
  }
  
  print("Combine the results from all chunks")
  gene_exp_merged <- rbindlist(result_list)
  
  print("4. Set rownames and drop the hgnc_symbol column")
  rownames(gene_exp_merged) <- gene_exp_merged$hgnc_symbol
  gene_exp_merged <- gene_exp_merged[, -1, with = FALSE]
  
  print("Cleanup")
  rm(result_list, gene_exp, unique_symbols)
  gc()
  saveRDS(gene_exp_merged, "gene_exp_merged.RDS")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print("a gene expression table made")
  
  seug <- CreateSeuratObject(counts = gene_exp_merged)
  seug <- NormalizeData(seug, normalization.method = "LogNormalize", scale.factor = 10000)
  #Variable features
  seug <- FindVariableFeatures(seug, selection.method = "vst", nfeatures = 2000)
  #Scale data
  seug <- ScaleData(seug, features = rownames(seug))
  #Linear dimension reduction - PCA
  seug <- RunPCA(seug, features = VariableFeatures(object = seug))
  #Finding the neighboring structure
  seug <- FindNeighbors(seug, dims = 1:15)
  #Cluster cells
  seug <- FindClusters(seug, resolution = 0.5)
  #Non linear dimension reduction
  seug <- RunUMAP(seug, dims = 1:15)
  print("gene exp seurat complete")
  saveRDS(seug, file = paste(runid, "seu_gene.RDS", sep = "_"))
  ##################################################################
}


#Create a seurat object for PAS expression
seu <- CreateSeuratObject(counts = expression_matrix)
seu

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
saveRDS(seu, file = paste(runid, "seu_pat.RDS", sep = "_"))

#Variable features
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
#Scale data
seu <- ScaleData(seu, features = rownames(seu))

saveRDS(seu, file = paste(runid, "seu_pat.RDS", sep = "_"))

#Linear dimension reduction - PCA
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
#Finding the neighboring structure
seu <- FindNeighbors(seu, dims = 1:15)
#Cluster cells
seu <- FindClusters(seu, resolution = 0.5)
#Non linear dimension reduction
seu <- RunUMAP(seu, dims = 1:15)
print("pas exp seurat complete")
saveRDS(seu, file = paste(runid, "seu_pat.RDS", sep = "_"))


p1 <- DimPlot(seug, label = T ) + ggtitle(paste("Gene exp based clusters (N=",length(unique(seug@meta.data$seurat_clusters)),")",sep = "")) + NoLegend()

p2 <- DimPlot(seu, label = T ) + ggtitle(paste("PAS exp based clusters (N=",length(unique(seu@meta.data$seurat_clusters)),")",sep = "")) + NoLegend()

pdf(paste(runid, "seurat_clusters.pdf", sep = "_"), width = 10, height = 5)
p1 + p2
dev.off()

library(GenomicRanges)
library(rtracklayer)

annotation <- import(gtf_file)

ann <- read.delim("annotatedpas.bed", header=F)
colnames(ann) <- c("chr", "start", "end", "PASid", "ensembl_gene_id", "hgnc_symbol", "strand")

ann$PASlabel <- paste("chr", ann$chr, ann$start, ann$end, ann$ensembl_gene_id, ann$hgnc_symbol,  sep="-")
print(head(ann))
dim(ann)

#select 3'-UTR regions only

# Convert to GRanges
roi_gr <- GRanges(
  seqnames = ann$chr,
  ranges = IRanges(start = as.numeric(ann$start), end = as.numeric(ann$end) ),
  strand = ann$strand,
  mcols = rownames(seu)# was ann$PASlabel
)
roi_gr

# Filter for 3'-UTR features
utr_3 <- annotation[annotation$type == "three_prime_utr"]

# Step 1: Extend the 3'-UTR regions
# Extend by 100 bp upstream and downstream
extended_utr <- resize(utr_3, width = width(utr_3) + 5000, fix = "start")

# Step 3: Find overlaps
overlaps <- findOverlaps(roi_gr, extended_utr)
# Get overlapping regions
overlapping_regions <- roi_gr[unique(queryHits(overlaps))]

# Display overlapping regions
print(overlapping_regions)

### Subset seurat object based on only PAS within 3'-UTRs AND worth keeping with expression
filtered.seu <- subset(seu, features = overlapping_regions$mcols)
dim(filtered.seu)

saveRDS(filtered.seu, file = paste(runid, "seu_3utr_pat.RDS", sep = "_"))

#Rscript ~/Dropbox/codes/ProjectEMA/test/pat_down.R /Users/yasinkaymaz/Documents/data/PASseqData/SRR8326029_EMA/emaout SRR8326029