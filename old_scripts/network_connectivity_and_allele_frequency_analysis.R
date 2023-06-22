#### AFR EURO ALLELE AND NETWORK COMPARISON ####

exomeAF <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/mean_enhancedAF_exome.tsv", col_names = TRUE)
network <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/disease_networks/mammary_epithelium_symbol.tsv")
network <- network %>% dplyr::select(-4)

## read in breast TCGA runs for afr and eur ## 

afr <-  readr::read_delim("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/paper_runs/ancestry_files/estimated_ancestry_afr_importance_scores.cut.txt", delim =" ")
afr.sig <- afr$Importance

eur <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/paper_runs/ancestry_files/estimated_eur_importance_scores.cuttxt") 
eur.sig <- eur$Importance


common.network1 <- network[((network$Gene1 %in% afr.sig) & (network$Gene2 %in% eur.sig)),]
common.network2 <- network[((network$Gene1 %in% eur.sig) & (network$Gene2 %in% afr.sig)),]
all.network <- rbind(common.network1, common.network2)

all.genes <- unique(c(all.network$Gene1, all.network$Gene2))

exome.in.net <- exomeAF[exomeAF$gene %in% all.genes,]


dat <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(dat) <- c("gene", "afr_AF", "eur_AF")

max.allele <- function(gene){
  temp <- exome.in.net[exome.in.net$gene == gene,]
  
  max.eur <- max(temp$nfe)
  max.afr <- max(temp$afr)
  the.row <- c(gene, max.afr, max.eur)
  dat[nrow(dat) + 1,] <- the.row
  return(dat)
}

my.d <- lapply(all.genes, max.allele)
my.d <- do.call(rbind,my.d)


write.table(all.network, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/afr_eur_sig.network.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(my.d, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/afr_eur_sig.enhanced_allele_freq.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)



all.network
conn.dat <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(conn.dat) <- c("Gene1", "Gene2", "Connection")

for(i in 1:length(all.genes)){
  gene <- all.genes[i]
  temp <- all.network[(all.network$Gene1 == gene | all.network$Gene2 == gene),]
  sorted.conn <- temp[order(temp$Connection, decreasing = TRUE),]
  keep.conn <- sorted.conn[1:3,]
  conn.dat <- rbind(conn.dat, keep.conn)
  #conn.dat[nrow(conn.dat) + 1,] <- the.row
}
write.table(conn.dat, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/afr_eur_sig.network_max_conn_3.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)

duplicate <- data.frame("Gene1" = conn.dat$Gene2, "Gene2" = conn.dat$Gene1, "Connection" = conn.dat$Connection) #doing this to makee sure all genes are represented in Gene1

alls <- merge(my.d, duplicate, by.x = "gene", by.y = "Gene1")
write.table(alls, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/afr_eur_sig.enhanced_allele_freq_and_max_conn_3.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)


