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

binning <- function(chrom) {
    # load file
    for (name in HM_list) {
        assign(name, PrepareChipseq(rtracklayer::import.bed(paste0(data_path, "/", name, "/", name, "_chr", chrom, ".bed"))))
    }

    # load gemone info
    anno_temp <- readr::read_delim(paste0(anno_path, "/anno_chr", chrom, ".csv"), delim = ",")

    # anno arrange
    anno <- filter(anno_temp, gene_biotype %in% c("protein_coding") & abs(start - end) > TSS_down + TTS_up)
    anno_po <- filter(anno, strand == "+")
    anno_po$TSS_up <- anno_po$start - TSS_up
    anno_po$TSS_down <- anno_po$start + TSS_down
    anno_po$TTS_up <- anno_po$end - TTS_up
    anno_po$TTS_down <- anno_po$end + TTS_down
    gene_po_list <- unique(anno_po$gene_name)
    anno_ne <- filter(anno, strand == "-")
    anno_ne$TSS_up <- anno_ne$end + TSS_up
    anno_ne$TSS_down <- anno_ne$end - TSS_down
    anno_ne$TTS_up <- anno_ne$start + TTS_up
    anno_ne$TTS_down <- anno_ne$start - TTS_down
    gene_ne_list <- unique(anno_ne$gene_name)

    # positive strand
    for (i in 1:length(gene_po_list)) {
        gene <- anno_po[anno_po$gene_name == gene_po_list[i],][1,]
        tss_bins <- TileSequence("1", gene$TSS_up, gene$TSS_down, binsize)
        tts_bins <- TileSequence("1", gene$TTS_up, gene$TTS_down, binsize)
        inter_bins <- GRanges(seqnames = Rle("1"),
                            ranges = IRanges(start = gene$TSS_down +1, end = gene$TTS_up - 1),
                            strand = Rle("*"))
        gene_bins <- Reduce("c", list(tss_bins, inter_bins, tts_bins))
        gene_bins$name <- rep(gene_po_list[i], length(gene_bins))
        if (i == 1) {
            out_po_bins <- gene_bins
        } else {
            out_po_bins <- Reduce("c", list(out_po_bins, gene_bins))
        }
    }
    # negative strand
    for (i in 1:length(gene_ne_list)) {
        gene <- anno_ne[anno_ne$gene_name == gene_ne_list[i],][1,]
        tss_bins <- TileSequence("1", gene$TSS_down, gene$TSS_up, binsize)
        tts_bins <- TileSequence("1", gene$TTS_down, gene$TTS_up, binsize)
        inter_bins <- GRanges(seqnames = Rle("1"),
                            ranges = IRanges(start = gene$TTS_up + 1, end = gene$TSS_down - 1),
                            strand = Rle("*"))
        gene_bins <- Reduce("c", list(tts_bins, inter_bins, tss_bins))
        gene_bins$name <- rep(gene_ne_list[i], length(gene_bins))
        if (i == 1) {
            out_ne_bins <- gene_bins
        } else {
            out_ne_bins <- Reduce("c", list(out_ne_bins, gene_bins))
        }
    }

    HM_po_temp_list <- paste0(HM_list, "_po_temp")
    for (name in HM_list) {
        assign(paste0(name, "_po_temp"), BinChipseq(eval(parse(text = name)), out_po_bins))
    }
    HM_ne_temp_list <- paste0(HM_list, "_ne_temp")
    for (name in HM_list) {
        assign(paste0(name, "_ne_temp"), BinChipseq(eval(parse(text = name)), out_ne_bins))
    }

    # positive strand
    for (i in 1:length(HM_list)) {
        out_df_temp <- data.frame(t(matrix(eval(parse(text = HM_po_temp_list[i]))$score, nrow = 81)))
        out_df_temp <- cbind(HM = rep(HM_list[i], nrow(out_df_temp)), out_df_temp)
        out_df_temp <- cbind(gene = gene_po_list, out_df_temp)
        if (i == 1) {
            out_df <- out_df_temp
        } else {
            out_df <- rbind(out_df, out_df_temp)
        } 
    }
    # negative strand
    for (i in 1:length(HM_list)) {
        out_df_temp <- data.frame(t(matrix(eval(parse(text = HM_ne_temp_list[i]))$score, nrow = 81)))
        out_df_temp <- out_df_temp[,ncol(out_df_temp):1]
        out_df_temp <- cbind(HM = rep(HM_list[i], nrow(out_df_temp)), out_df_temp)
        out_df_temp <- cbind(gene = gene_ne_list, out_df_temp)
        out_df <- rbind(out_df, out_df_temp)
    }
    return(out_df)
}

# parameters
HM_list <- c("h3k4me3", "h3k4me1", "h3k9me3", "h3k36me3")
chr_list <- c(1:19, "X", "Y")
data_path <- "/mnt/c/chipseq/test8"
anno_path <- "/mnt/c/others/anno_mouse"
output_path <- "/mnt/c/machinelearning/test2"
TSS_up <- 2000
TSS_down <- 2000
TTS_up <- 2000
TTS_down <- 2000
binsize <- 100
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# procession
cat("---------------processing---------------", "\n")
cl <- makeForkCluster(2)
out_df_list <- parLapply(cl, chr_list, binning)
output <- Reduce("rbind", out_df_list)
stopCluster(cl)

colnames(output)[1] <- "gene" 
colnames(output)[2] <- "HM" 
output <- arrange(output, gene, HM)

# output
cat("-------------writing csv-------------", "\n")
write.csv(output, file.path(output_path, "output.csv"), row.names = FALSE)