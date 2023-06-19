###NETWORK ANNOTATION ###
library(org.Hs.eg.db)
library(tibble)
require(readr)  # for read_csv()
require(dplyr) #data frame handling 
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
network_file <- args[1]
out_file <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/disease_networks/",network_file,"_symbol.tsv")

network <- readr::read_tsv(network_file, col_names = FALSE)
annot.df <- data.frame("Symbols" = mapIds(org.Hs.eg.db, keys = as.character(network$X2), column = "SYMBOL", keytype = "ENTREZID"), network)
annot.df <- annot.df %>% dplyr::select(-X2)
colnames(annot.df) <- c("Gene2", "X1", "X3")
annot.df <- data.frame("Symbols" = mapIds(org.Hs.eg.db, keys = as.character(annot.df$X1), column = "SYMBOL", keytype = "ENTREZID"), annot.df)
annot.df <- annot.df %>% dplyr::select(-X1)
colnames(annot.df) <- c("Gene1","Gene2","Connection")
rownames(annot.df) <- NULL

write_delim(annot.df, file = out_file, 
            delim = '\t', col_names = TRUE)


