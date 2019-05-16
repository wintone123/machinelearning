# library
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(parallel)

# functions
max_legend <- function(a) {
    cond <- TRUE
    n <- 0
    while (cond) {
        n <- n + 1
        b = a %/% 2
        if (b == 1) {
            cond <- FALSE
            return(n + 1)
        } else {
            a <- b
        }
    }
}

# imaging <- function(mt) {
#     pheatmap(mt, cluster_rows = FALSE, cluster_cols = FALSE,
#              show_rownames = FALSE, show_colnames = FALSE, 
#              col = col(length(breaks)), breaks = breaks, legend = FALSE, 
#              border_color = NA,
#              filename = paste0(path_out, "/", gene_name, ".png"))
# }

imaging_par <- function(gene_name) {
    temp_mt <- as.matrix(filter(imput, gene == gene_name)[1:4,3:83])
    if (sum(temp_mt) != 0) {
        pheatmap(temp_mt, cluster_rows = FALSE, cluster_cols = FALSE,
                 show_rownames = FALSE, show_colnames = FALSE, 
                 col = col(length(breaks)), breaks = breaks, legend = FALSE, 
                 border_color = NA,
                 filename = paste0(path_out, "/", gene_name, ".png"))
    } 
}

# parameter
path_in <- "/mnt/c/machinelearning/test2/output.csv"
path_out <- "/mnt/c/machinelearning/test2/gene_image"
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# load file
cat("-------------loading csv-------------", "\n")
imput <- readr::read_csv(path_in, col_names = TRUE)
imput_mt <- as.matrix(imput[,3:83])
breaks <- seq(0, max_legend(max(imput_mt)))

# process
cat("--------------processing--------------", "\n")
if (dir.exists(path_out) == FALSE) {
    dir.create(path_out)
}
cl <- makeForkCluster(2)
parLapply(cl, unique(imput$gene), imaging_par)
stopCluster(cl)
# for (gene_name in unique(imput$gene)) {
#     temp_mt <- as.matrix(filter(imput, gene == gene_name)[,3:83])
#     if (sum(temp_mt) != 0) {
#         imaging(temp_mt)
#     }
# }