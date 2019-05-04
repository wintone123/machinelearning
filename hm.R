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

# BinChipseq <- function(reads, input, bins){
#     test <- countOverlaps(bins, reads) / length(reads) * 1000000
#     input <- countOverlaps(bins, input) / length(inpput) * 1000000
#     mcols(bins)$score = ifelse(test > input, test - input, 0)
#     return(bins)
# }

BinChipseq <- function(reads, bins){
    mcols(bins)$score = countOverlaps(bins, reads) / length(reads) * 1000000
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

# mapping <- function(temp_list, gene_name, output_path) {
#     for (i in 1:length(temp_list)) {
#         if (i == 1) {
#             mat_df <- data.frame(matrix(eval(parse(text = temp_list[i]))$score, nrow = 1), row.names = temp_list[i])
#         } else {
#             mat_df <- rbind(mat_df ,data.frame(matrix(eval(parse(text = temp_list[i]))$score, nrow = 1), row.names = temp_list[i]))
#         }
#     }
#     mat <- as.matrix(mat_df)
#     mat <- t(mat)
#     breaks <- c(min(mat), max(mat))
#     pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
# 			 show_rownames = FALSE, show_colnames = FALSE, 
#              col = col(length(breaks)), breaks = breaks, legend = FALSE, border_color = NA,
#              filename = paste0(output_path, "/", gene_name, ".png"))
# }

# binning <- function(gene_name, anno, binsize, output_path, h3k4me3, h3k4me3_in) {
#     gene <- anno[anno$gene_name == gene_name,][1,]
#     tss_bins <- TileSequence(gene$gene_name, gene$TSS_up, gene$TSS_down, binsize)
#     tts_bins <- TileSequence(gene$gene_name, gene$TTS_up, gene$TTS_down, binsize)
#     inter_bins <- GRanges(seqnames = gene$gene_name,
#                           ranges = IRanges(start = gene$TSS_down +1, end = gene$TTS_up - 1),
#                           strand = Rle("*"))
#     gene_bins <- Reduce("c", list(tss_bins, inter_bins, tts_bins))
#     h3k4me3_temp <- BinChipseq(h3k4me3, h3k4me3_in ,gene_bins)
#     mapping(h3k4me3_temp, output_path)
# }

# parameters
HM_list <- c("h3k4me3", "h3k4me1", "h3k9me3", "h3k36me3")
data_path <- "/mnt/c/chipseq/test8"
anno_path <- "/mnt/c/others/anno_chr1.csv"
output_path <- "/mnt/c/machinelearning/test1/gene_image"
TSS_up <- 2000
TSS_down <- 2000
TTS_up <- 2000
TTS_down <- 2000
binsize <- 100
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# load file
cat("-----------preparing file-----------", "\n")
for (name in HM_list) {
    assign(name, PrepareChipseq(rtracklayer::import.bed(paste0(data_path, "/", name, "/", name, "_chr1.bed"))))
}

# load gemone info
cat("-----------preparing anno-----------", "\n")
anno_temp <- readr::read_delim(anno_path, delim = ",", col_names = TRUE)

# anno arrange
anno <- filter(anno_temp, biotype %in% c("protein_coding") & abs(TTS - TSS) > TSS_down + TTS_up)
# anno <- data.frame(gene_name = anno_temp$gene_name,
#                    TSS = ifelse(anno_temp$strand == "+", anno_temp$TSS, anno_temp$TTS),
#                    TTS = ifelse(anno_temp$strand == "+", anno_temp$TTS, anno_temp$TSS))
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
HM_temp_list <- paste0(HM_list, "_temp")
for (i in 1:length(gene_list)) {
    gene <- anno[anno$gene_name == gene_list[i],][1,]
    tss_bins <- TileSequence("1", gene$TSS_up, gene$TSS_down, binsize)
    tts_bins <- TileSequence("1", gene$TTS_up, gene$TTS_down, binsize)
    inter_bins <- GRanges(seqnames = Rle("1"),
                          ranges = IRanges(start = gene$TSS_down +1, end = gene$TTS_up - 1),
                          strand = Rle("*"))
    gene_bins <- Reduce("c", list(tss_bins, inter_bins, tts_bins))
    gene_bins$name <- rep(gene_list[i], length(gene_bins))
    if (i == 1) {
        out_bins <- gene_bins
    } else {
        out_bins <- Reduce("c", list(out_bins, gene_bins))
    }
}
for (name in HM_list) {
    assign(paste0(name, "_temp"), BinChipseq(eval(parse(text = name)), out_bins))
}
for (i in 1:length(HM_list)) {
    out_df_temp <- data.frame(t(matrix(eval(parse(text = HM_temp_list[i]))$score, nrow = 81)))
    out_df_temp <- cbind(HM = rep(HM_list[i], nrow(out_df_temp)), out_df_temp)
    out_df_temp <- cbind(gene = gene_list, out_df_temp)
    if (i == 1) {
        out_df <- out_df_temp
    } else {
        out_df <- rbind(out_df, out_df_temp)
    }    
}
colnames(out_df)[1] <- "gene" 
colnames(out_df)[2] <- "HM" 
out_df <- arrange(out_df, gene, HM)

# output
cat("------------writing csv-------------", "\n")
write.csv(out_df, file.path(output_path, "output.csv"), row.names = FALSE)