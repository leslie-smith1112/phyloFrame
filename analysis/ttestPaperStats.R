### Code for all t.tests performed in the PhyloFrame paper. The paragraph they are in and their results is before the chunk of relevant code ### 

# Because EAF is calculated from healthy tissue, this approach can integrate information from 
# under-represented populations that are not present in a smaller disease specific databases, including TCGA.
# To confirm this, we compared EAF in COSMIC cancer-related genes to non-COSMIC genes and did not find an
# enrichment in EAF (Supplementary Fig.~\ref{supp_fig:EAF_density}; t.test, df = 16863435, two-tailed (Welch), p-value = 1), 
# as expected, given EAF is calculated from healthy tissue.
ccgenes <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/Cancer_census_genes.tsv")
exomeAF <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/mean_enhancedAF_exome.tsv", col_names = TRUE)
cc.exome <- exomeAF[exomeAF$gene %in% ccgenes$`Gene Symbol`,]
temp.exome <- exomeAF %>% select(-c("chrom", "position","rs_id","ref_allele", "alt_allele", "consequence", "type","distance"))

#seperate cosmic and not cosmic genes and make melt them to get all EAFs in value column 
cc.edit <- temp.exome[temp.exome$gene %in% ccgenes$`Gene Symbol`,]
long.cc <- melt(cc.edit, id.vars = "gene",variable.name = "ancestry")
head(long.cc)
not.cc <- temp.exome[!(temp.exome$gene %in% ccgenes$`Gene Symbol`),]
long.not <- melt(not.cc, id.vars = "gene",variable.name = "ancestry")
head(long.not)

t.test(long.not$value, long.cc$value) 

########################################################################################################
# COSMIC genes with high EAFs are more frequently enriched in African and East Asian but not European 
# ancestries (Supplementary Fig.~\ref{supp_fig:EAF_density}B; t-test, df = 255930, two-tailed (Welch) p-value < 2.2e-16 ) .
afr.eeas <- c(cc.exome$eas,cc.exome$afr) 
t.test(afr.eeas, cc.exome$nfe)

########################################################################################################
# For example, while all of the BRCA models have high AUC (0.99-1), there is significantly less overlap 
# in the disease signatures identified by the benchmark model compared to the PhyloFrame models trained 
# on the same data (Fig.~\ref{fig:grid_AUC}G; t-test, df = 193.46, two-tailed (Welch) $p-value = 6.494e-06$).

set.seed(123)
disease <- "breast" #define disease and directory of results 
dir <- "/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal"

eur.num <- 1:17
afr.num <- 1:2
eas.num <- 1:1
all.num <- 20
admixed.num <- 1:1
mixed.num <- 1:6

eur.model.names <- paste0("eur_model_",eur.num)
afr.model.names <- paste0("afr_model_",afr.num)
eas.model.names <- paste0("eas_model_",eas.num)
admixed.names <- paste0("admixed_model_",admixed.num) # only one model 
mixed.names <- paste0("mixed_model_", mixed.num)

## read in signatures from phyloframe and bechmark ancestry models 
pf.eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/eur/")
pf.afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/afr/")
pf.admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/admixed/")
pf.eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/eas/")
pf.mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/mixed/")

eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/eur/")
afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/afr/")
admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/admixed/")
eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/eas/")
mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/mixed/")

pf.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
pf.eur.myfiles <- lapply(pf.eur.names, readr::read_tsv)
pf.afr.names <- paste0(pf.afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
pf.afr.myfiles <- lapply(pf.afr.names, readr::read_tsv)
pf.eas.names <- paste0(pf.eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt") 
pf.eas.myfiles <- lapply(pf.eas.names, readr::read_tsv)

df.eur.names <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
eur.myfiles <- lapply(df.eur.names, readr::read_tsv)
df.afr.names <- paste0(afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
afr.myfiles <- lapply(df.afr.names, readr::read_tsv)
df.eas.names <- paste0(eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt") 
eas.myfiles <- lapply(df.eas.names, readr::read_tsv)

### SIGNAURE SIGNATURE CORRELATION ## 
# ------------------- EUR ------------------- # 
# data frame for phyloFrame and benchmark signtature percentages - 20 models total (17 eur, 2 afr, 1 eas)
eur.df <- data.frame(matrix(ncol = all.num,nrow = 0))
pf.eur.df <- data.frame(matrix(ncol = all.num,nrow = 0))
#for every ancestry signature from every model - divide intersection/union with signatures * 100 to get
#% overlap of each signature and build matrix
# do for each ancestry model
#*.myfiles is the list of the ancestry model signatures
for(i in 1:eur.count){
  temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  pf.temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  for(j in 1:eur.count){
    ### check the intersection of model i with model j 
    the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable))
    ### check how many genes total from the 2 models 
    gene.union <- union(eur.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100  
    temp[1,j] <- overlap
    ## phyloframe ## 
    pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j] <- pf.overlap
  }
  for(j in 1:afr.count){
    the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable))
    gene.union <- union(eur.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+eur.count] <- overlap 
    ## phyloframe ## 
    pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j+eur.count] <- pf.overlap  
  }
  for(j in 1:eas.count){
    the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable))
    gene.union <- union(eur.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+ (eur.count + afr.count)] <- overlap 
    ## phyloframe ## 
    pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j + (eur.count + afr.count) ] <- pf.overlap 
  }
  eur.df <- rbind(eur.df, temp)
  pf.eur.df <- rbind(pf.eur.df, pf.temp)
}
rownames(eur.df) <- eur.model.names
rownames(pf.eur.df) <- eur.model.names
rownames(eur.df) <- paste0(rownames(eur.df))
rownames(pf.eur.df) <- paste0(rownames(pf.eur.df))
c.name <- c(eur.model.names, afr.model.names,eas.model.names)
colnames(eur.df) <- c.name
colnames(pf.eur.df) <- c.name

afr.df <- data.frame(matrix(ncol = all.num,nrow = 0))
pf.afr.df <- data.frame(matrix(ncol = all.num,nrow = 0))
#### AFR
for(i in 1:afr.count){
  temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  pf.temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  for(j in 1:eur.count){
    the.int <- length(intersect(afr.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable))
    gene.union <- union(afr.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j] <- overlap
    ## phyloframe 
    pf.the.int <- length(intersect(pf.afr.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.afr.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j] <- pf.overlap
  }
  for(j in 1:afr.count){
    the.int <- length(intersect(afr.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable))
    gene.union <- union(afr.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+ eur.count] <- overlap
    ## phyloFrame 
    pf.the.int <- length(intersect(pf.afr.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.afr.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j+eur.count] <- pf.overlap
  }
  for(j in 1:eas.count){
    the.int <- length(intersect(afr.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable))
    gene.union <- union(afr.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+ (eur.count + afr.count)] <- overlap
    ## phyloframe ## 
    pf.the.int <- length(intersect(pf.afr.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.afr.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j + (eur.count + afr.count) ] <- pf.overlap
  }
  afr.df <- rbind(afr.df, temp)
  pf.afr.df <- rbind(pf.afr.df, pf.temp)
}
rownames(afr.df) <- afr.model.names
rownames(pf.afr.df) <- afr.model.names
rownames(afr.df) <- paste0(rownames(afr.df))
rownames(pf.afr.df) <- paste0(rownames(pf.afr.df))
colnames(afr.df) <- c.name
colnames(pf.afr.df) <- c.name

#### EAS
eas.df <- data.frame(matrix(ncol = all.num,nrow = 0))
pf.eas.df <- data.frame(matrix(ncol = all.num,nrow = 0))

for(i in 1:eas.count){
  temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  pf.temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  for(j in 1:eur.count){
    ### checking model1 of phyloframe against all other models of benchamrk
    the.int <- length(intersect(eas.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable))
    gene.union <- union(eas.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j] <- overlap
    ## phyloFrame ## 
    pf.the.int <- length(intersect(pf.eas.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eas.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j] <- pf.overlap
  }
  for(j in 1:afr.count){
    ### checking model1 of phyloframe against all other models of benchamrk
    the.int <- length(intersect(eas.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable))
    gene.union <- union(eas.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+eur.count] <- overlap
    # phyloframe 
    pf.the.int <- length(intersect(pf.eas.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eas.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j+ eur.count] <- pf.overlap
  }
  for(j in 1:eas.count){
    ### checking model1 of phyloframe against all other models of benchamrk
    the.int <- length(intersect(eas.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable))
    gene.union <- union(eas.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+ (eur.count + afr.count)] <- overlap
    ## phyloframe ## 
    pf.the.int <- length(intersect(pf.eas.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eas.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j + (eur.count + afr.count) ] <- pf.overlap
  }
  eas.df <- rbind(eas.df, temp)
  pf.eas.df <- rbind(pf.eas.df, pf.temp)
}
# 
rownames(eas.df) <- eas.model.names
rownames(pf.eas.df) <- eas.model.names
colnames(eas.df) <- c.name
colnames(pf.eas.df) <- c.name

#bind all ancestry dataframes and create row labels 
benchmark <-rbind(eur.df, afr.df,eas.df)
rownames(benchmark) <- paste0("benchmark_", rownames(benchmark))
phyloFrame <- rbind(pf.eur.df, pf.afr.df, pf.eas.df)
rownames(phyloFrame) <- paste0("PF_", rownames(phyloFrame))


## correlation matrix 
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
bench.mat <- as.matrix(benchmark)
pf.mat <- as.matrix(phyloFrame)
# # get correlation matrix
res.bm <- cor(bench.mat)
res.pf <- cor(pf.mat)
round(res.bm, 2)
round(res.pf, 2)
#
# rcorr(x, type = c("pearson","spearman"))
# install.packages("Hmisc")
library("Hmisc")
res2.pf <- rcorr(as.matrix(pf.mat))
res2.bm <- rcorr(as.matrix(bench.mat))
# #flatted to keep only cor coefficient and p value
pf  <- flattenCorrMatrix(res2.pf$r, res2.pf$P)
bm  <- flattenCorrMatrix(res2.bm$r, res2.bm$P)
t.test(pf$cor, bm$cor)

########################################################################################################
# COSMIC genes captured by PhyloFrames BRCA model 
# (Supplementary Fig.~\ref{supp_fig:COSMIC_dotPlots}) consistently 
# have a higher EAF in ancestries not found in the training data, most
# significantly in South Asians (two-sided (Welch) t.test, $df = 137,504$, $p-value = 1.126e-10$).

## asssumes data loaded from previous statistic
## read in the additional admixed and mixed signatures 
admixed.num <- 1
pf.admixed.names <- paste0(pf.admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt")
pf.admixed.myfiles <- lapply(pf.admixed.names, readr::read_tsv)
mixed.num <- 1:6
pf.mixed.names <- paste0(pf.mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt")
pf.mixed.myfiles <- lapply(pf.mixed.names, readr::read_tsv)

df.admixed.names <- paste0(admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt")
admixed.myfiles <- lapply(df.admixed.names, readr::read_tsv)
df.mixed.names <- paste0(mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt")
mixed.myfiles <- lapply(df.mixed.names, readr::read_tsv)

## PF ## 
eur.sig <- do.call(rbind,pf.eur.myfiles)
afr.sig <- do.call(rbind,pf.afr.myfiles)
eas.sig <- do.call(rbind,pf.eas.myfiles)
admixed.sig <- do.call(rbind, pf.admixed.myfiles)
mixed.sig <- do.call(rbind, pf.mixed.myfiles)
all.pf.sig <- rbind(eur.sig, afr.sig, eas.sig, admixed.sig, mixed.sig)

## BM ## 
bm.eur.sig <- do.call(rbind,eur.myfiles)
bm.afr.sig <- do.call(rbind,afr.myfiles)
bm.eas.sig <- do.call(rbind,eas.myfiles)
bm.admixed.sig <- do.call(rbind, admixed.myfiles)
bm.mixed.sig <- do.call(rbind, mixed.myfiles)
all.bm.sig <- rbind(bm.eur.sig, bm.afr.sig, bm.eas.sig,bm.admixed.sig, bm.mixed.sig)

# combine all gene signature genes and look at mutations in the included COSMIC genes 
pf.genes <- exomeAF[exomeAF$gene %in% all.pf.sig$Variable,]
bm.genes <- exomeAF[exomeAF$gene %in% all.bm.sig$Variable,]
bm.genes <- bm.genes[bm.genes$gene %in% ccgenes$`Gene Symbol`,]
pf.genes <- pf.genes[pf.genes$gene %in% ccgenes$`Gene Symbol`,]
t.test(pf.genes$sas, bm.genes$sas)

########################################################################################################
# These genes are shared at much higher rates between models in PhyloFrame; 76\% COSMIC 
# genes identified in PhyloFrame models are found in multiple signatures compared to 
# only 21\% COSMIC genes identified in more than one benchmark signature.
# Each PhyloFrame BRCA model signature identifies more COSMIC genes than its benchmark counterpart 
# (13 vs 8 COSMIC genes on average; Supplementary Fig.~\ref{supp_fig:COSMIC_dotPlots}, two-sided (Welch) t.test, $df = 31.268$, $p-value = 1.079e-06$).

## asssumes data loaded from previous statistic
# get mean of cosmic genes in pf and benchmark signatures for cosmic genes
dat <- data.frame(matrix(ncol = 2, nrow = 0))
##EUR
for(i in 1:length(pf.eur.myfiles)){
  current <- pf.eur.myfiles[i][[1]]$Variable
  common <- length(intersect(current, ccgenes$`Gene Symbol`))
  current.bm <- eur.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, ccgenes$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

##AFR
for(i in 1:length(pf.afr.myfiles)){
  current <- pf.afr.myfiles[i][[1]]$Variable
  common <- length(intersect(current, ccgenes$`Gene Symbol`))
  current.bm <- afr.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, ccgenes$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

##EAS
for(i in 1:length(pf.eas.myfiles)){
  current <- pf.eas.myfiles[i][[1]]$Variable
  common <- length(intersect(current, ccgenes$`Gene Symbol`))
  current.bm <- eas.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, ccgenes$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

t.test(dat$phyloFrame, dat$benchmark)
mean(dat$phyloFrame)
mean(dat$benchmark)

