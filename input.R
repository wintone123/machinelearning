# library
library(parallel)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(GenomicRanges)

# functions
PrepareChipseq <- function(reads){
  frag_len = median(chipseq::estimate.mean.fraglen(reads))
  reads_extended = SummarizedExperiment::resize(reads, width = frag_len)
  return(trim(reads_extended))
}

BinChipseq <- function(reads, input, bins){
  test <- countOverlaps(bins, reads) / length(reads)
  input <- countOverlaps(bins, input) / length(inpput)
  mcols(bins)$score = ifelse(test > input, test - input, 0)
  return(bins)
}

TileSequence <- function(seqname, start, end, tilewidth){
    start_list <- seq(start, end, by = tilewidth)
    end_list <- start_list + tilewidth -1
    if (start_list[length(start_list)] == end) {
        end_list <- end_list[c(1:length(end_list)-1)]
        end_list[length(end_list)] <- start_list[length(start_list)]
        start_list <- start_list[c(1:length(start_list)-1)]
    } else {
        end_list[length(end_list)] <- end
    }
    GRanges_temp <- GRanges(seqnames = Rle(seqname),
                          ranges = IRanges(start = start_list, end = end_list),
                          strand = Rle("*"))
    return(GRanges_temp)
}

mapping <- function(h3k4me3_temp, output_path) {
    mat <- matrix(h3k4me3_temp$score, nrow = 1)
    breaks <- c(min(mat), max(mat))
    pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
			 show_rownames = FALSE, show_colnames = FALSE, 
             col = col(length(breaks)), breaks = breaks, legend = FALSE, border_color = NA,
             filename = paste0(output_path, "/", h3k4me3@seqnames[1], "_h3k4me3", ".png"))
}

binning <- function(gene_name, anno, binsize, output_path, h3k4me3, h3k4me3_in) {
    gene <- anno[anno$gene_name == gene_name,][1,]
    tss_bins <- TileSequence(gene$gene_name, gene$TSS_up, gene$TSS_down, binsize)
    tts_bins <- TileSequence(gene$gene_name, gene$TTS_up, gene$TTS_down, binsize)
    inter_bins <- GRanges(seqnames = gene$gene_name,
                          ranges = IRanges(start = gene$TSS_down +1, end = gene$TTS_up - 1),
                          strand = Rle("*"))
    gene_bins <- Reduce("c", list(tss_bins, inter_bins, tts_bins))
    h3k4me3_temp <- BinChipseq(h3k4me3, h3k4me3_in ,gene_bins)
    mapping(h3k4me3_temp, output_path)
}

# parameters
test_path <- "/mnt/c/chipseq/test2/split_files/rep1_best_chr1.bed"
control_path <- "/mnt/c/chipseq/test2/split_files/control_best_chr1.bed"
anno_path <- "/mnt/c/others/anno_chr1.csv"
output_path <- "/mnt/c/machinelearning/test1"
TSS_up <- 2000
TSS_down <- 2000
TTS_up <- 2000
TTS_down <- 2000
binsize <- 100
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# load file
cat("-----------preparing file-----------", "\n")
h3k4me3 <- rtracklayer::import.bed(test_path)
h3k4me3_in <- rtracklayer::import.bed(control_path)

# bed arrange
h3k4me3 <- PrepareChipseq(h3k4me3)
h3k4me3_in <- PrepareChipseq(h3k4me3_in)

# load gemone info
cat("-----------preparing anno-----------", "\n")
anno_temp <- readr::read_delim(anno_path, delim = ",", col_names = TRUE)

# anno arrange
anno_temp <- filter(anno_temp, biotype %in% c("protein_coding"))
anno <- data.frame(gene_name = anno_temp$gene_name,
                   TSS = ifelse(anno_temp$strand == "+", anno_temp$TSS, anno_temp$TTS),
                   TTS = ifelse(anno_temp$strand == "+", anno_temp$TTS, anno_temp$TSS))
anno$TSS_up <- anno$TSS - TSS_up
anno$TSS_down <- anno$TSS + TSS_down
anno$TTS_up <- anno$TTS - TTS_up
anno$TTS_down <- anno$TTS + TTS_down
gene_list <- unique(anno$gene_name)

# procession
cat("-------------processing-------------", "\n")
# cl <- makeCluster(2)
# clusterExport(cl, c("anno", "binsize", "output_path", "h3k4me3", "h3k4me3_in", 
#                     "TileSequence", "mapping", "BinChipseq"), envir = globalenv())
# clusterEvalQ(cl, "GenomicRanges")
# parLapply(cl, gene_list, binning)
# stopCluster(cl)
for (gene_name in gene_list) {
    cat(gene_name, "\n")
    gene <- anno[anno$gene_name == gene_name,][1,]
    print(gene)
    tss_bins <- TileSequence("chr1", gene$TSS_up, gene$TSS_down, binsize)
    tts_bins <- TileSequence("chr1", gene$TTS_up, gene$TTS_down, binsize)
    inter_bins <- GRanges(seqnames = Rle("chr1"),
                          ranges = IRanges(start = gene$TSS_down +1, end = gene$TTS_up - 1),
                          strand = Rle("*"))
    gene_bins <- Reduce("c", list(tss_bins, inter_bins, tts_bins))
    print(gene_bins)
    h3k4me3_temp <- BinChipseq(h3k4me3, h3k4me3_in ,gene_bins)
    print(h3k4me3_temp)
    mapping(h3k4me3_temp, output_path)
}