### - FIGURE 3 G - ### 

############ SIGNATURE SIGNATURE CORRELATION PLOT  ##############
library(ComplexHeatmap)
set.seed(123)
disease <- "breast" #define disease and directory of results 

dir <- "/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal"

############ EUR ############
## define lists for number of models in each ancestry for help reading in files ## 

if(disease == "breast"){ #TODO this is the worst possible way to do this - needs
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
}else if(disease == "thyroid"){
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
}else if(disease == "uterine"){
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
}else{
  stop("must enter valid disease")
}
# get number of models for each ancestry to help read in signature results 
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
# note: we dont include the admixed and mixed signatures in this plot 
# pf.admixed.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
# pf.admixed.myfiles <- lapply(pf.eur.names, readr::read_tsv)
# pf.mixed.names <- paste0(pf.afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
# pf.afr.myfiles <- lapply(pf.afr.names, readr::read_tsv)

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

all <- rbind(phyloFrame, benchmark)
#all <- rownames_to_column(all, "model_num")
all.mat <- as.matrix(all)
df.row <- c(rep("eur",eur.count), rep("afr",afr.count), rep("eas",eas.count),rep("eur",eur.count), rep("afr",afr.count), rep("eas",eas.count))
df <- c(rep("eur",eur.count), rep("afr",afr.count), rep("eas",eas.count))

#column annotation
for.annotation <- data.frame(Ancestry = df)
rownames(for.annotation) <- c.name
for.annotation$Ancestry <- as.factor(for.annotation$Ancestry)
#row annotation
for.row.annotation <- data.frame(Ancestry = df.row)
rownames(for.row.annotation) <- rownames(all)
for.row.annotation$Ancestry <- as.factor(for.row.annotation$Ancestry)

## correlation matrix 
res <- cor(all.mat)
round(res, 2)

anc.color <- c("#CA4136","#1F619E","#496849")
names(anc.color) <- c("eur","afr","eas")
the.color <- list(Ancestry = anc.color)
library(RColorBrewer)
display.brewer.all() 
# display.brewer.pal(n = 9, name = 'YlOrBr')
# brewer.pal(n = 9, name = "Reds")
library(pheatmap)
heatmap <- pheatmap(
  all.mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_colors = the.color,
  annotation_col = for.annotation,
  annotation_row = for.row.annotation,
  main = "Model Signature Overlap",
  colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(25)
  
)

heatmap

## correlation matrix 
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#   )
# }
# bench.mat <- as.matrix(benchmark)
# pf.mat <- as.matrix(phyloFrame)
# # get correlation matrix 
# res.bm <- cor(bench.mat)
# res.pf <- cor(pf.mat)
# round(res.bm, 2)
# round(res.pf, 2)
# 
# rcorr(x, type = c("pearson","spearman"))
# install.packages("Hmisc")
# library("Hmisc")
# res2.pf <- rcorr(as.matrix(pf.mat))
# res2.bm <- rcorr(as.matrix(bench.mat))
# #flatted to compare cor coefficient and p value 
# pf  <- flattenCorrMatrix(res2.pf$r, res2.pf$P)
# bm  <- flattenCorrMatrix(res2.bm$r, res2.bm$P)
# 
# t.test(pf$cor, bm$cor)

#write.table(all, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/figures/heatmap_dataframe.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
#all <- read.table(file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/figures/heatmap_dataframe.tsv", sep = "\t")

#write.table(cc.df, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/model_sig_cancer_consensus.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)




