######################## NETWORK FIGURE #2 WITH BETWEENESS CENTRALITY AND HOW OFTEN IN SIG ######################## 

disease <- "breast"

#dir <- "/pf08_03thyroid_10000_V2_expr700_Netowrk_and_allele_change_ANCVAR_poster_version"

dir <- "/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal"
################ define lists for number of models in each ancestry for help reading in files ################ 
## BREAST 
eur.num <- list(1:17)
eur.num <- unlist(eur.num)
afr.num <- list(1:2)
afr.num <- unlist(afr.num)
eas.num <- list(1:1)
eas.num <- unlist(eas.num)
all.num <- 20
admixed.num <- list(1:1)
admixed.num <- unlist(admixed.num)
mixed.num <- list(1:6)
mixed.num <- unlist(mixed.num)

## THYROID:
eur.num <- list(1:23)
eur.num <- unlist(eur.num)
afr.num <- list(1:1)
afr.num <- unlist(afr.num)
eas.num <- list(1:3)
eas.num <- unlist(eas.num)
all.num <- 27
admixed.num <- list(1:1)
admixed.num <- unlist(admixed.num)
mixed.num <- list(1:9)
mixed.num <- unlist(mixed.num)

## UTERINE
eur.num <- list(1:12)
eur.num <- unlist(eur.num)
afr.num <- list(1:2)
afr.num <- unlist(afr.num)
eas.num <- 1
all.num <- 15
admixed.num <- list(1:1)
admixed.num <- unlist(admixed.num)
mixed.num <- list(1:6)
mixed.num <- unlist(mixed.num)

eur.count <- length(eur.num)
afr.count <- length(afr.num)
eas.count <- length(eas.num)
admixed.count <- length(admixed.num)
mixed.count <- length(mixed.num)

eur.model.names <- paste0("eur_model_",eur.num)
afr.model.names <- paste0("afr_model_",afr.num)
eas.model.names <- paste0("eas_model_",eas.num)
admixed.names <- paste0("admixed_model_",admixed.num) # only one model 
mixed.names <- paste0("mixed_model_", mixed.num)

# disease <- "thyroid/pf4_thyroid_new_node2_5000V2_rescale500"

##################  read in signatures from phyloframe and bechmark ancestry models  ##################
pf.eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/eur/")
pf.afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/afr/")
pf.admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/admixed/")
pf.eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/eas/")
pf.mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/mixed/")

pf.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
pf.eur.myfiles <- lapply(pf.eur.names, readr::read_tsv)
pf.afr.names <- paste0(pf.afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
pf.afr.myfiles <- lapply(pf.afr.names, readr::read_tsv)
pf.eas.names <- paste0(pf.eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt") 
pf.eas.myfiles <- lapply(pf.eas.names, readr::read_tsv)
pf.admixed.names <- paste0(pf.admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt") 
pf.admixed.myfiles <- lapply(pf.admixed.names, readr::read_tsv)
pf.mixed.names <- paste0(pf.mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt") 
pf.mixed.myfiles <- lapply(pf.mixed.names, readr::read_tsv)

################## read in network ##################
network <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/disease_networks/mammary_epithelium_symbol.tsv")
network <- network %>% dplyr::select(-4)

################## read in cosmic genes and bind ancestry signatures ##################
ccgenes <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/Cancer_census_genes.tsv")
eur.sig <- do.call(rbind,pf.eur.myfiles)
afr.sig <- do.call(rbind,pf.afr.myfiles)
eas.sig <- do.call(rbind,pf.eas.myfiles)
# admixed.sig <- do.call(rbind, pf.admixed.myfiles)
# mixed.sig <- do.call(rbind, pf.mixed.myfiles)

eur.cc.sig <- eur.sig[eur.sig$Variable %in% ccgenes$`Gene Symbol`,]
eur.cc.sig.u <- unique(eur.cc.sig$Variable)
afr.cc.sig <- afr.sig[afr.sig$Variable %in% ccgenes$`Gene Symbol`,]
afr.cc.sig.u <- unique(afr.cc.sig$Variable)
eas.cc.sig <- eas.sig[eas.sig$Variable %in% ccgenes$`Gene Symbol`,]
eas.cc.sig.u <- unique(eas.cc.sig$Variable)
# admixed.cc.sig <- admixed.sig[admixed.sig$Variable %in% ccgenes$`Gene Symbol`,]
# admixed.cc.sig.u <- unique(admixed.cc.sig$Variable)
# mixed.cc.sig <- mixed.sig[mixed.sig$Variable %in% ccgenes$`Gene Symbol`,]
# mixed.cc.sig.u <- unique(mixed.cc.sig$Variable)

all.genes.pfmods <- c(eur.sig, afr.sig, eas.sig)
all.genes.pfmods <- all.genes.pfmods$Variable

all.genes <- c(eas.cc.sig.u,eur.cc.sig.u,afr.cc.sig.u)

c.net <- network[network$Gene1 %in% ccgenes$ | network$Gene2 %in% all.genes,]
cc.network2 <- network[network$Gene1 %in% all.genes | network$Gene2 %in% all.genes,]
#write.table(cc.network2,"/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/figures/network_BRCA_Cosmic_sig_asc_networkConn1.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
cc.network2 <- cc.network2[cc.network2$Connection > 0.99,]


length(unique(c(cc.network2$Gene1,cc.network2$Gene2)))
head(cc.network2)
table(cc.network2$Gene2)

## -- some of the dz genes are not in network -> trim them out -- ##
dat.graph <- graph.data.frame(cc.network2, directed =  FALSE)
idk <- betweenness(dat.graph, v = V(dat.graph), directed = FALSE,  weights = NULL)

idk.dat <- as.data.frame(idk)
idk.dat <- rownames_to_column(idk.dat, "gene")


dat <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(dat) <- c("gene","eur", "afr","eas")
for(i in 1:length(all.genes)){
  current.gene <- all.genes[i]
  eur.total <- eur.cc.sig[eur.cc.sig$Variable == current.gene,]
  included <- nrow(eur.total)
  eur.percent <- included/length(eur.num) * 100
  
  afr.total <- afr.cc.sig[afr.cc.sig$Variable == current.gene,]
  included <- nrow(afr.total)
  afr.percent <- included/length(afr.num) * 100
  
  eas.total <- eas.cc.sig[eas.cc.sig$Variable == current.gene,]
  included <- nrow(eas.total)
  eas.percent <- included/length(eas.num) * 100
  
  new.row <- c(current.gene, eur.percent, afr.percent, eas.percent)
  dat <- rbind(dat,new.row)
}


eur.num <- list(1:17)
eur.num <- unlist(eur.num)
afr.num <- list(1:2)
afr.num <- unlist(afr.num)
eas.num <- list(1:1)

write.table(dat,"/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/figures/network_BRCA_Cosmic_sig.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)




