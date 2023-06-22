### 
#### GRAPHING MODEL PEFORMANCE #### 
#### THESE ARE THE SCATTERPLOT PHYLOFRAME VS BENCHMARK FOR ALL BATCHES IN EACH MODEL ## 
require(tibble)
require(readr)  # for read_csv()
require(dplyr) #data frame handling 
require(tidyr)
require(stringi)
library(ggplot2)
set.seed(1234)
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/scaled_importance_matrix.R")
#disease <- "thyroid"
#ancestry.sample.mixed <- "e"
#ancestry.sample <- "eur"
#ancestry.model <- "eur"
#dir <- "pen1_version1_3"
#run(disease, dir, ancestry.sample, ancestry.model, 0.5)
# 
# run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eur","eur",0.8)
# 
# 
#run(disease, dir, ancestry.sample, ancestry.model, mixture)

run <- function(disease, dir, ancestry.sample, ancestry.model, mixture){
  temp.list <- list.files(paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/samples/",ancestry.sample,"/"))
  num <- length(temp.list)
  num
  
  all <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(all) <- c("model.num", "ancestry", "benchmark", "phyloFrame")
  
  df.path <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",dir,"/model_runs/",ancestry.model)
  pf.path <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",dir,"/model_runs/phyloFrame/",ancestry.model)
  
  mod.num <- list(1:num)
  mod.num <- unlist(mod.num)
  

  df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"afr_metrics.tsv") #CHANGE HERE
  pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"afr_metrics.tsv")#CHANGE HERE
  myfiles <- lapply(df.names, readr::read_tsv)
  pffiles <- lapply(pf.names, readr::read_tsv)
  
  # - for benchmark model - #
  afr.mod <- do.call("rbind", myfiles)
  afr.mod <- afr.mod[afr.mod$.metric == "roc_auc",]
  
  # - phyloframe model - # 
  pf.afrmod <- do.call("rbind", pffiles)
  pf.afrmod <- pf.afrmod[pf.afrmod$.metric == "roc_auc",]
  # 
  # model.num <- paste0("model_",mod.num)
  # ancestry <- rep("afr", num) #CHANGE HERE AND BELOW 
  # new.afr.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
  #                          "benchmark" = afr.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
  # all <- rbind(all, new.afr.df)

    ###ADMIXED 
  df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"admixed_metrics.tsv") #CHANGE HERE
  pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"admixed_metrics.tsv")#CHANGE HERE
  myfiles <- lapply(df.names, readr::read_tsv)
  pffiles <- lapply(pf.names, readr::read_tsv)
  
  #- for benchmark model - #
  admixed.mod <- do.call("rbind", myfiles)
  admixed.mod <- admixed.mod[admixed.mod$.metric == "roc_auc",]
  
  #- phyloframe model -# 
  pf.admixedmod <- do.call("rbind", pffiles)
  pf.admixedmod <- pf.admixedmod[pf.admixedmod$.metric == "roc_auc",]


  ####EAS 
  df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eas_metrics.tsv") #CHANGE HERE
  pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"eas_metrics.tsv")#CHANGE HERE
  
  myfiles <- lapply(df.names, readr::read_tsv)
  pffiles <- lapply(pf.names, readr::read_tsv)
  
  #- for benchmark model - #
  eas.mod <- do.call("rbind", myfiles)
  eas.mod <- eas.mod[eas.mod$.metric == "roc_auc",]
  #- phyloframe model -# 
  pf.easmod <- do.call("rbind", pffiles)
  pf.easmod <- pf.easmod[pf.easmod$.metric == "roc_auc",]


  #### EUR 
  df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eur_metrics.tsv") #CHANGE HERE
  pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"eur_metrics.tsv")#CHANGE HERE
  
  myfiles <- lapply(df.names, readr::read_tsv)
  pffiles <- lapply(pf.names, readr::read_tsv)
  
  #- for benchmark model - #
  eur.mod <- do.call("rbind", myfiles)
  eur.mod <- eur.mod[eur.mod$.metric == "roc_auc",]
  #- phyloframe model -# 
  pf.eurmod <- do.call("rbind", pffiles)
  pf.eurmod <- pf.eurmod[pf.eurmod$.metric == "roc_auc",]

  ####MIXED 
  df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"mixed_metrics.tsv") #CHANGE HERE
  pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"mixed_metrics.tsv")#CHANGE HERE
  
  myfiles <- lapply(df.names, readr::read_tsv)
  pffiles <- lapply(pf.names, readr::read_tsv)
  
  #- for benchmark model - #
  mixed.mod <- do.call("rbind", myfiles)
  mixed.mod <- mixed.mod[mixed.mod$.metric == "roc_auc",]
  #- phyloframe model -# 
  pf.mixedmod <- do.call("rbind", pffiles)
  pf.mixedmod <- pf.mixedmod[pf.mixedmod$.metric == "roc_auc",]

  
  ### performance t.test ### 
  ## afr 
  print("AFR")
  head(pf.mod)
  head(afr.mod)
  print(t.test(pf.afrmod$.estimate, afr.mod$.estimate) )# eur thyroid not significant 
  
## admixed 
  print("ADMIXED")
  print(t.test(pf.admixedmod$.estimate, admixed.mod$.estimate)) #eur thyroid not significant 

## eas 
  print("EAS")
  print(t.test(pf.easmod$.estimate, eas.mod$.estimate)) #eur thyroid not significant 
## eur 
  print("EUR")
  print(t.test(pf.eurmod$.estimate, eur.mod$.estimate)) #eur thyroid not significant 
## mixed 
  print("MIXED")
  print(t.test(pf.mixedmod$.estimate, mixed.mod$.estimate)) #eur thyroid not significant 
  
  ## all 
  pf.all <- rbind(pf.admixedmod, pf.afrmod, pf.easmod, pf.eurmod, pf.mixedmod)
  bm.all <- rbind(admixed.mod, afr.mod, eas.mod, eur.mod, mixed.mod)
  print("ALL")
  print(t.test(pf.all$.estimate, bm.all$.estimate))
  to.return <- list("benchmark" = afr.mod, "phyloFrame" = pf.afrmod)
  return(to.return)
} 



run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "admixed","admixed",0.8)
run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "afr","afr",0.8)
run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eas","eas",0.8)
run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eur","eur",0.8)
thyroid.mix <- run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "mixed_ancestry","mixed",0.8)
thyroid.bm <- thyroid.mix$benchmark
thyroid.pf <- thyroid.mix$phyloFrame

run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "admixed","admixed",0.8)
run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "afr","afr",0.8)
run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eas","eas",0.8)
run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eur","eur",0.8)
brca.mix <- run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "mixed_ancestry","mixed",0.8)
brca.bm <- brca.mix$benchmark
brca.pf <- brca.mix$phyloFrame


run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "admixed","admixed",0.8)
run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "afr","afr",0.8)
run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eas","eas",0.8)
run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eur","eur",0.8)
uterine.mix <- run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "mixed_ancestry","mixed",0.8)
uterine.bm <- uterine.mix$benchmark
uterine.pf <- uterine.mix$phyloFrame

temp.bm <- rbind(uterine.bm, thyroid.bm)
temp.pf <- rbind(thyroid.pf, thyroid.pf)
t.test(temp.bm$.estimate, temp.pf$.estimate)
