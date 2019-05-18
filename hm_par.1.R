# library
library(parallel)
library(dplyr)
library(tidyr)
library(GenomicRanges)

# functions
PrepareChipseq <- function(reads){
    frag_len = median(chipseq::estimate.mean.fraglen(reads))
    reads_extended = SummarizedExperiment::resize(reads, width = frag_len)
    return(trim(reads_extended))
}

BinChipseq <- function(reads, bins){
    mcols(bins)$score = countOverlaps(bins, reads) / length(reads) * 1000000
    return(bins)
}

TileSequence <- function(start, end, tilewidth){
    start_list <- seq(start, end, by = tilewidth)
    end_list <- start_list + tilewidth -1
    if (start_list[length(start_list)] == end) {
        end_list <- end_list[c(1:length(end_list)-1)]
        end_list[length(end_list)] <- start_list[length(start_list)]
        start_list <- start_list[c(1:length(start_list)-1)]
    } else {
        end_list[length(end_list)] <- end
    }
    out <- list("start" = start_list, "end" = end_list)
    return(out)
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
    gene_po_list <- unique(anno_po$gene_name)
    anno_ne <- filter(anno, strand == "-")
    gene_ne_list <- unique(anno_ne$gene_name)

    # positive strand
    for (i in 1:length(gene_po_list)) {
        gene <- anno_po[anno_po$gene_name == gene_po_list[i],][1,]
        up_list <- TileSequence(gene$start-TSS_up, gene$start+TSS_down, binsize)
        down_list <- TileSequence(gene$end-TTS_up, gene$end+TTS_down, binsize)
        po_temp_df <- data.frame(start = c(up_list$start, gene$start, down_list$start),
                                 end = c(up_list$end, gene$end, down_list$end))
        if (i == 1) {
            out_po_df <- po_temp_df
        } else {
            out_po_df <- rbind(out_po_df, po_temp_df)
        }
    }
    out_po_bins <- GRanges(seqnames = Rle(chrom),
                           ranges = IRanges(start = out_po_df$start, end = out_po_df$end),
                           strand = Rle("*"))
    HM_po_temp_list <- paste0(HM_list, "_po_temp")
    for (name in HM_list) {
        assign(paste0(name, "_po_temp"), BinChipseq(eval(parse(text = name)), out_po_bins))
    }
    for (i in 1:length(HM_list)) {
        out_df_temp <- data.frame(t(matrix(eval(parse(text = HM_po_temp_list[i]))$score, nrow = binnum)))
        out_df_temp <- cbind(HM = rep(HM_list[i], nrow(out_df_temp)), out_df_temp)
        out_df_temp <- cbind(gene = gene_po_list, out_df_temp)
        if (i == 1) {
            out_df <- out_df_temp
        } else {
            out_df <- rbind(out_df, out_df_temp)
        } 
    }

    # negative strand
    for (i in 1:length(gene_ne_list)) {
        gene <- anno_ne[anno_ne$gene_name == gene_ne_list[i],][1,]
        up_list <- TileSequence(gene$start-TTS_up, gene$start+TTS_down, binsize)
        down_list <- TileSequence(gene$end-TSS_up, gene$end+TSS_down, binsize)
        ne_temp_df <- data.frame(start = c(down_list$start, gene$start, up_list$start),
                                 end = c(down_list$end, gene$end, up_list$end))
        if (i == 1) {
            out_ne_df <- ne_temp_df
        } else {
            out_ne_df <- rbind(out_ne_df, ne_temp_df)
        } 
    }
    out_ne_bins <- GRanges(seqnames = Rle(chrom),
                           ranges = IRanges(start = out_ne_df$start, end = out_ne_df$end),
                           strand = Rle("*"))
    HM_ne_temp_list <- paste0(HM_list, "_ne_temp")
    for (name in HM_list) {
        assign(paste0(name, "_ne_temp"), BinChipseq(eval(parse(text = name)), out_ne_bins))
    }
    for (i in 1:length(HM_list)) {
        out_df_temp <- data.frame(t(matrix(eval(parse(text = HM_ne_temp_list[i]))$score, nrow = binnum)))
        out_df_temp <- out_df_temp[,ncol(out_df_temp):1]
        out_df_temp <- cbind(HM = rep(HM_list[i], nrow(out_df_temp)), out_df_temp)
        out_df_temp <- cbind(gene = gene_ne_list, out_df_temp)
        out_df <- rbind(out_df, out_df_temp)
    }
    return(out_df)
}

# parameters
HM_list <- c("h3k4me1", "h3k4me1", "h3k9me3", "h3k36me3")
chr_list <- c(1:19, "X", "Y")
data_path <- "/mnt/c/chipseq/test8"
anno_path <- "/mnt/c/others/anno_mouse"
output_path <- "/mnt/c/machinelearning/test3"
TSS_up <- 2000
TSS_down <- 2000
TTS_up <- 2000
TTS_down <- 2000
binsize <- 100
binnum <- length(seq(-TSS_up, TSS_down, by = binsize)) -1 + length(seq(-TTS_up, TTS_down, by = binsize)) -1 + 1

# procession
cat("-------------processing--------------", "\n")
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