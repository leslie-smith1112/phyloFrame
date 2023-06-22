############## SUPP FIG 6 ###### 

source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breast.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/thyroid.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/uterine.R")

require(readr)  # for read_csv()
require(dplyr) #data frame handling 
library(tidyr)
library(stringi)
library(ggplot2)
####################################  DEFINE DISEASE AND MODEL NUMBERS FOR EACH ANCESTRY #################################### 

set.seed(123)


to.return <- run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal","eur")
to.return <- run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal","eas")
to.return <- run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal","admixed")
to.return <- run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal","afr")
to.return <- run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal","mixed")

to.return <- run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal","eur")
to.return <- run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal","eas")
to.return <- run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal","admixed")
to.return <- run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal","afr")
to.return <- run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal","mixed")

to.return <- run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal","eur")
to.return <- run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal","admixed")
to.return <- run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal","afr")
to.return <- run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal","mixed")
head(to.return)

##### BRIEF ADD FOR CHECKING UTERINE SEROUS SAMPLES #### 
# clinical <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/ucec_tcga_pan_can_atlas_2018_clinical_data.tsv")
# clin <- clinical %>% dplyr::select(`Sample ID`,`Tumor Type`)
# table(clin$`Tumor Type`)
# keeping <- to.return %>% dplyr::select(sample_id,percent_correct)
# test <- merge(keeping, clin, by.x = "sample_id", by.y = "Sample ID")
# head(test)
# table(test$percent_correct,test$`Tumor Type`)

# estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
# samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == "UCEC",]
# samples.ancestry <- samples.ancestry %>% dplyr::select(patient, consensus_ancestry)
# samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
# ancestry.ucec <- merge(samples.ancestry, clin, by.x = "patient", by.y = "Sample ID")
# NOTE: EUR model each sample predicted N times, 760 Endometrial predicted 100 a few at 70-90, Serous: between 10-45 % correctly predicted 
# ADM one model endometrical 50 wrong 674 right, serous 29 wrong 151 right 
#AFr endometrial 757 haf 100% right, 6 had 0 % right and 9 had 50%, Serous had 71 had 0% right, 60 had 100% right, 79 had 50% right 
#MIXED endometiral had 756 correct and only 16 below that, serous had 49 wrong and all over tha map basically, only 5 got 100% 

to.return$percent_correct <- as.numeric(to.return$percent_correct)
to.return$consensus_ancestry <- factor(to.return$consensus_ancestry, levels = c("admixedPF", "admixedBM", "afrPF","afrBM", "easPF","easBM","eurPF","eurBM"))

p <- ggplot(to.return, aes(x=consensus_ancestry, y=percent_correct, fill = consensus_ancestry, color = consensus_ancestry)) + geom_boxplot() + 
 scale_fill_manual(values = c("#FDC652","#ffe943","#1F619E", "#9fccfa","#496849","#63b34c","#CA4136","#f5a2a5"
                                                        )) + scale_color_manual(values = c("#ca9e43","#cab836","#153e64", "#7b9dbf","#2b3d2b","#427533","#7c2a23","#b97d7f"))+ theme_minimal() 
p

## FOR MODELS WITH ONLY ONE MODEL ### 
## for eas 
to.return$percent_correct <- as.numeric(to.return$percent_correct)
to.return$consensus_ancestry <- factor(to.return$consensus_ancestry, levels = c("admixedPF", "admixedBM", "afrPF","afrBM", "eurPF","eurBM"))

p <- ggplot(to.return, aes(x=consensus_ancestry, y=percent_correct, fill = consensus_ancestry, color = consensus_ancestry)) + geom_boxplot() +
  scale_fill_manual(values = c("#FDC652","#ffe943","#1F619E", "#9fccfa","#CA4136","#f5a2a5"
  )) + scale_color_manual(values = c("#ca9e43","#cab836","#153e64", "#7b9dbf","#7c2a23","#b97d7f"))+ theme_minimal() 
p

## for admixed
to.return$percent_correct <- as.numeric(to.return$percent_correct)
to.return$consensus_ancestry <- factor(to.return$consensus_ancestry, levels = c("afrPF","afrBM", "easPF","easBM","eurPF","eurBM"))

p <- ggplot(to.return, aes(x=consensus_ancestry, y=percent_correct, fill = consensus_ancestry, color = consensus_ancestry)) + geom_boxplot() +
  scale_fill_manual(values = c("#1F619E", "#9fccfa","#496849","#63b34c","#CA4136","#f5a2a5"
  )) + scale_color_manual(values = c("#153e64", "#7b9dbf","#2b3d2b","#427533","#7c2a23","#b97d7f"))+ theme_minimal() 
p

## afr 
to.return$percent_correct <- as.numeric(to.return$percent_correct)
to.return$consensus_ancestry <- factor(to.return$consensus_ancestry, levels = c("admixedPF", "admixedBM", "easPF","easBM","eurPF","eurBM"))

p <- ggplot(to.return, aes(x=consensus_ancestry, y=percent_correct, fill = consensus_ancestry, color = consensus_ancestry)) + geom_boxplot() +
  scale_fill_manual(values = c("#FDC652","#ffe943","#496849","#63b34c","#CA4136","#f5a2a5"
  )) + scale_color_manual(values = c("#ca9e43","#cab836","#2b3d2b","#427533","#7c2a23","#b97d7f"))+ theme_minimal() 
p




####################################  read in ancestry samples from different batches (each ancestry has different # of batches) #################################### 
## function returns data frame with prediction % per sample for the ancestry model for benchmark and phyloFrame to be plotted  
run <- function(disease, cur.dir, ancestry){
  if (disease == "breast"){
    cancer.type <- "BRCA"
    subtype1 <- "Basal"
    subtype2 <- "Luminal"
    eur.num <- list(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4), 
                    rep(10,4), rep(11,4), rep(12,4), rep(13,4), rep(14,4), rep(15,4), rep(16,4), rep(17,4))
    eur.num <- unlist(eur.num)
    afr.num <- list(rep(1,4), rep(2,4))
    afr.num <- unlist(afr.num)
    eas.num <- list(rep(1,3))
    eas.num <- unlist(eas.num)
    admixed.num <- list(rep(1,3))
    admixed.num <- unlist(admixed.num)
    mixed.num <- list(rep(1,4),rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4))
    mixed.num <- unlist(mixed.num)
    eur.anc <- list("eur", "afr", "eas", "admixed")#17 models 
    afr.anc <- list("eur", "afr", "eas", "admixed")#2 models 
    eas.anc <- list("eur", "afr", "admixed") #1 model, not enought eas to run eas model on eas samples 
    admixed.anc <- list("eur", "afr", "eas") #1 model 
    mixed.anc <- list("eur", "afr", "eas", "admixed") # 6 models 
  }else if(disease == "thyroid"){
    cancer.type <- "THCA"
    subtype1 <- "M0"
    subtype2 <- "MX"
    eur.num <- list(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4), 
                    rep(10,4), rep(11,4), rep(12,4), rep(13,4), rep(14,4), rep(15,4), rep(16,4), rep(17,4), rep(18,4), rep(19,4), rep(20,4), rep(21,4), rep(22,4), rep(23,4))
    eur.num <- unlist(eur.num)
    afr.num <- list(rep(1,3))
    afr.num <- unlist(afr.num)
    eas.num <- list(rep(1,4),rep(2,4), rep(3,4))
    eas.num <- unlist(eas.num)
    admixed.num <- list(rep(1,3))
    admixed.num <- unlist(admixed.num)
    mixed.num <- list(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4))
    mixed.num <- unlist(mixed.num)
    eur.anc <- list("eur", "afr", "eas", "admixed")# 23 models 
    afr.anc <- list("eur", "eas", "admixed") #1 model 
    eas.anc <- list("eur", "afr", "eas", "admixed") # 3 models 
    admixed.anc <- list("eur", "eas", "afr") #1 model 
    mixed.anc <- list("eur", "eas", "afr","admixed") #1 model
    
  }else if(disease == "uterine"){
    cancer.type <- "UCEC"
    subtype1 <- "Endometrioid"
    subtype2 <- "Serous"
    ## UTERINE
    eur.num <- list(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4), 
                    rep(10,4), rep(11,4), rep(12,4))
    eur.num <- unlist(eur.num)
    afr.num <- list(rep(1,4), rep(2,4))
    afr.num <- unlist(afr.num)
    admixed.num <- list(rep(1,3))
    admixed.num <- unlist(admixed.num)
    mixed.num <- list(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4))
    mixed.num <- unlist(mixed.num)
    #not ennough samples for uterine east asians for model 
    #uterine 
    eur.anc <- list("eur", "afr", "eas", "admixed")#12 models 
    afr.anc <- list("eur", "afr", "eas", "admixed") #2 models 
    admixed.anc <- list("eur", "eas", "afr") #1 model 
    mixed.anc <- list("eur", "eas", "afr","admixed") #6 model

    
  }else{
    print("Please enter a valid disease, currently they are: 1. breast 2. thyroid 3. uterine")
    continue <- 0
  }
  
  
  ####################################  READ IN ANCESTRY INFORMATION #################################### 
  estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
  samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == cancer.type,] 
  samples.ancestry <- samples.ancestry[,-4:-8]
  samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
  ####################################  ####################################
  
  ## do once for phyloFrame then once for benchmark (below)
  ## PHYLOFRAME
  if(ancestry == "admixed"){
    pf.admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/admixed/")
    admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/admixed/")
    #################################### EUR MODEL  PHYLOFRAME####################################
    ## read in results for every ancestry  that admixed model was tested on##
    pf.df.admixed.names <- paste0(pf.admixed.dir,"model_",admixed.num,"/model_",admixed.num, admixed.anc,"_results.tsv") 
    pf.admixed.myfiles <- lapply(pf.df.admixed.names, readr::read_tsv)
    ## this gives every models prediction of every sample tested no ## 
    pf.admixed <- do.call(rbind,pf.admixed.myfiles) 
    
    ## add ancestry percentages information to samples ##
    all.admixed <- merge(pf.admixed, samples.ancestry, by.x="sample_id", by.y="patient")
    all.admixed$Model <- "PF.Eur"
    table(all.admixed$consensus_ancestry)
    all.admixed$percent_correct <- NA
    admix <- c("afr_admix","eas_admix","eur_admix", "admix")
    all.admixed$consensus_ancestry[all.admixed$consensus_ancestry %in% admix] <- "admixed"
    ## new dataframe for only keeping one instance of each sample ## 
    admixed.new <- data.frame(matrix(ncol = 9, nrow = 0))
    
    colnames(admixed.new) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                               "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
    ## for each sample, find % correctly idenfied by the model and add to dataframes ## 
    samples <- unique(all.admixed$sample_id)
    for(i in 1:length(samples)){
      temp.dat <- all.admixed[all.admixed$sample_id == samples[i],] #get all predictions of sample
      total <- nrow(temp.dat)
      correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
      percent <- (correct/total) * 100
      all.admixed[all.admixed$sample_id == samples[i],]$percent_correct <- percent
      new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],
                    temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], 
                    temp.dat$`Admixture % SAS`[1], temp.dat$Model[1],percent )
      admixed.new[nrow(admixed.new) + 1,]<-new.list
    }
    
    
    #### BENCHMARK ### 
    df.admixed.names <- paste0(admixed.dir,"model_",admixed.num,"/model_",admixed.num, admixed.anc,"_results.tsv") 
    admixed.myfiles <- lapply(df.admixed.names, readr::read_tsv)
    admixed.bm <- do.call(rbind,admixed.myfiles)
    
    all.admixed.bm <- merge(admixed.bm, samples.ancestry, by.x="sample_id", by.y="patient")
    all.admixed.bm$Model <- "Eur"
    all.admixed.bm$percent_correct <- NA
    all.admixed.bm$consensus_ancestry[all.admixed.bm$consensus_ancestry %in% admix] <- "admixed"
    admixed.new.bm <- data.frame(matrix(ncol = 9, nrow = 0))
    colnames(admixed.new.bm) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                                  "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
    
    samples <- unique(all.admixed.bm$sample_id)
    for(i in 1:length(samples)){
      temp.dat <- all.admixed.bm[all.admixed.bm$sample_id == samples[i],]
      total <- nrow(temp.dat)
      correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
      percent <- (correct/total) * 100
      all.admixed.bm[all.admixed.bm$sample_id == samples[i],]$percent_correct <- percent
      new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],temp.dat$`Admixture % AMR`[1],
                    temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], temp.dat$`Admixture % SAS`[1], 
                    temp.dat$Model[1],percent)
      admixed.new.bm[nrow(admixed.new.bm) + 1,]<-new.list
    }
    admixed.new.bm$consensus_ancestry <- paste0(admixed.new.bm$consensus_ancestry, "BM")
    admixed.new$consensus_ancestry <- paste0(admixed.new$consensus_ancestry, "PF")
    to.return <- rbind(admixed.new.bm, admixed.new)
    
  }else if(ancestry == "afr"){##############################AFRICAN ########################################## 
    pf.afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/afr/")
    afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/afr/")
    pf.df.afr.names <- paste0(pf.afr.dir,"model_",afr.num,"/model_",afr.num, afr.anc,"_results.tsv") 
    pf.afr.myfiles <- lapply(pf.df.afr.names, readr::read_tsv)
    pf.afr <- do.call(rbind,pf.afr.myfiles)
    
    ## add ancestry information to samples 
    all.afr <- merge(pf.afr, samples.ancestry, by.x="sample_id", by.y="patient")
    all.afr$Model <- "PF.Afr"
    all.afr$percent_correct <- NA
    admix <- c("afr_admix","eas_admix","eur_admix", "admix")
    all.afr$consensus_ancestry[all.afr$consensus_ancestry %in% admix] <- "admixed"
    ## new dataframe for keeping only 1 instance of each sample 
    afr.new <- data.frame(matrix(ncol = 9, nrow = 0))
    colnames(afr.new) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                           "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
    samples <- unique(all.afr$sample_id)
    for(i in 1:length(samples)){
      temp.dat <- all.afr[all.afr$sample_id == samples[i],]
      total <- nrow(temp.dat)
      correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
      percent <- (correct/total) * 100
      all.afr[all.afr$sample_id == samples[i],]$percent_correct <- percent
      new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],
                    temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], 
                    temp.dat$`Admixture % SAS`[1], temp.dat$Model[1],percent )
      afr.new[nrow(afr.new) + 1,]<-new.list
    }
    
    #### BENCHMARK #### 
    df.afr.names <- paste0(afr.dir,"model_",afr.num,"/model_",afr.num, afr.anc,"_results.tsv") 
    afr.myfiles <- lapply(df.afr.names, readr::read_tsv)
    afr.bm <- do.call(rbind,afr.myfiles)
    
    all.afr.bm <- merge(afr.bm, samples.ancestry, by.x="sample_id", by.y="patient")
    all.afr.bm$Model <- "Afr"
    all.afr.bm$percent_correct <- NA
    all.afr.bm$consensus_ancestry[all.afr.bm$consensus_ancestry %in% admix] <- "admixed"
    afr.new.bm <- data.frame(matrix(ncol = 9, nrow = 0))
    
    colnames(afr.new.bm) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                              "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
    samples <- unique(all.afr.bm$sample_id)
    for(i in 1:length(samples)){
      temp.dat <- all.afr.bm[all.afr.bm$sample_id == samples[i],]
      total <- nrow(temp.dat)
      correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
      percent <- (correct/total) * 100
      all.afr.bm[all.afr.bm$sample_id == samples[i],]$percent_correct <- percent
      new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],
                    temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], 
                    temp.dat$`Admixture % SAS`[1], temp.dat$Model[1],percent )
      afr.new.bm[nrow(afr.new.bm) + 1,]<-new.list
    }
    afr.new.bm$consensus_ancestry <- paste0(afr.new.bm$consensus_ancestry, "BM")
    afr.new$consensus_ancestry <- paste0(afr.new$consensus_ancestry, "PF")
    to.return <- rbind(afr.new, afr.new.bm)
    
  }else if (ancestry == "eas"){#############################EAST ASIAN ###########################################
    pf.eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/eas/")
    eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/eas/")
    pf.df.eas.names <- paste0(pf.eas.dir,"model_",eas.num,"/model_",eas.num, eas.anc,"_results.tsv") 
    pf.eas.myfiles <- lapply(pf.df.eas.names, readr::read_tsv)
    pf.eas <- do.call(rbind,pf.eas.myfiles)
    
    ## add ancestry information to samples 
    all.eas <- merge(pf.eas, samples.ancestry, by.x="sample_id", by.y="patient")
    all.eas$Model <- "PF.Eas"
    all.eas$percent_correct <- NA 
    admix <- c("afr_admix","eas_admix","eur_admix", "admix")
    all.eas$consensus_ancestry[all.eas$consensus_ancestry %in% admix] <- "admixed"
    eas.new <- data.frame(matrix(ncol = 9, nrow = 0))
    colnames(eas.new) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                           "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
    
    samples <- unique(all.eas$sample_id)
    for(i in 1:length(samples)){
      temp.dat <- all.eas[all.eas$sample_id == samples[i],]
      total <- nrow(temp.dat)
      correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
      percent <- (correct/total) * 100
      all.eas[all.eas$sample_id == samples[i],]$percent_correct <- percent
      new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],
                    temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], 
                    temp.dat$`Admixture % SAS`[1], temp.dat$Model[1],percent )
      eas.new[nrow(eas.new) + 1,]<-new.list
    }
    
    #### BENCHMARK #### 
    df.eas.names <- paste0(eas.dir,"model_",eas.num,"/model_",eas.num, eas.anc,"_results.tsv") 
    eas.myfiles <- lapply(df.eas.names, readr::read_tsv)
    eas.bm <- do.call(rbind,eas.myfiles)
    
    all.eas.bm <- merge(eas.bm, samples.ancestry, by.x="sample_id", by.y="patient")
    all.eas.bm$Model <- "EasBM"
    all.eas.bm$percent_correct <- NA
    all.eas.bm$consensus_ancestry[all.eas.bm$consensus_ancestry %in% admix] <- "admixed"
    eas.new.bm <- data.frame(matrix(ncol = 9, nrow = 0))
    
    
    colnames(eas.new.bm) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                              "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
    
    samples <- unique(all.eas.bm$sample_id)
    for(i in 1:length(samples)){
      temp.dat <- all.eas.bm[all.eas.bm$sample_id == samples[i],]
      total <- nrow(temp.dat)
      correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
      percent <- (correct/total) * 100
      all.eas.bm[all.eas.bm$sample_id == samples[i],]$percent_correct <- percent
      new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],
                    temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], 
                    temp.dat$`Admixture % SAS`[1], temp.dat$Model[1],percent )
      eas.new.bm[nrow(eas.new.bm) + 1,]<-new.list
    }
    eas.new.bm$consensus_ancestry <- paste0(eas.new.bm$consensus_ancestry, "BM")
    eas.new$consensus_ancestry <- paste0(eas.new$consensus_ancestry, "PF")
    to.return <- rbind(eas.new, eas.new.bm)
    
  }else if(ancestry == "eur"){#####################################EUR#######################################################################
    pf.eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/eur/")
    eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/eur/")
    ## read in results for every ancestry  that eur model was tested on##
    pf.df.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num, eur.anc,"_results.tsv") 
    pf.eur.myfiles <- lapply(pf.df.eur.names, readr::read_tsv)
    ## this gives every models prediction of every sample tested no ## 
    pf.eur <- do.call(rbind,pf.eur.myfiles) 
    
    ## add ancestry percentages information to samples ##
    all.eur <- merge(pf.eur, samples.ancestry, by.x="sample_id", by.y="patient")
    all.eur$Model <- "PF.Eur"
    all.eur$percent_correct <- NA
    admix <- c("afr_admix","eas_admix","eur_admix", "admix")
    all.eur$consensus_ancestry[all.eur$consensus_ancestry %in% admix] <- "admixed"
    ## new dataframe for only keeping one instance of each sample ## 
    eur.new <- data.frame(matrix(ncol = 9, nrow = 0))
    
    colnames(eur.new) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                           "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
    ## for each sample, find % correctly idenfied by the model and add to dataframes ## 
    samples <- unique(all.eur$sample_id)
    for(i in 1:length(samples)){
      temp.dat <- all.eur[all.eur$sample_id == samples[i],] #get all predictions of sample
      total <- nrow(temp.dat)
      correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
      percent <- (correct/total) * 100
      all.eur[all.eur$sample_id == samples[i],]$percent_correct <- percent
      new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],
                    temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], 
                    temp.dat$`Admixture % SAS`[1], temp.dat$Model[1],percent )
      eur.new[nrow(eur.new) + 1,]<-new.list
    }
    
  
  #### BENCHMARK ### 
  df.eur.names <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num, eur.anc,"_results.tsv") 
  eur.myfiles <- lapply(df.eur.names, readr::read_tsv)
  eur.bm <- do.call(rbind,eur.myfiles)
  
  all.eur.bm <- merge(eur.bm, samples.ancestry, by.x="sample_id", by.y="patient")
  all.eur.bm$Model <- "Eur"
  all.eur.bm$percent_correct <- NA
  all.eur.bm$consensus_ancestry[all.eur.bm$consensus_ancestry %in% admix] <- "admixed"
  eur.new.bm <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(eur.new.bm) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                            "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")

  samples <- unique(all.eur.bm$sample_id)
  for(i in 1:length(samples)){
    temp.dat <- all.eur.bm[all.eur.bm$sample_id == samples[i],]
    total <- nrow(temp.dat)
    correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
    percent <- (correct/total) * 100
    all.eur.bm[all.eur.bm$sample_id == samples[i],]$percent_correct <- percent
    new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],temp.dat$`Admixture % AMR`[1],
                  temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], temp.dat$`Admixture % SAS`[1], 
                  temp.dat$Model[1],percent)
    eur.new.bm[nrow(eur.new.bm) + 1,]<-new.list
  }
  eur.new.bm$consensus_ancestry <- paste0(eur.new.bm$consensus_ancestry, "BM")
  eur.new$consensus_ancestry <- paste0(eur.new$consensus_ancestry, "PF")
  to.return <- rbind(eur.new.bm, eur.new)

  }else if(ancestry == "mixed"){#################################MIXED #######################################
    pf.mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/mixed/")
    mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/mixed/")
    ## read in results for every ancestry  that mixed model was tested on##
    pf.df.mixed.names <- paste0(pf.mixed.dir,"model_",mixed.num,"/model_",mixed.num, mixed.anc,"_results.tsv") 
    pf.mixed.myfiles <- lapply(pf.df.mixed.names, readr::read_tsv)
    ## this gives every models prediction of every sample tested no ## 
    pf.mixed <- do.call(rbind,pf.mixed.myfiles) 
    
    ## add ancestry percentages information to samples ##
    all.mixed <- merge(pf.mixed, samples.ancestry, by.x="sample_id", by.y="patient")
    all.mixed$Model <- "PF.Eur"
    all.mixed$percent_correct <- NA
    admix <- c("afr_admix","eas_admix","eur_admix", "admix")
    all.mixed$consensus_ancestry[all.mixed$consensus_ancestry %in% admix] <- "admixed"
    ## new dataframe for only keeping one instance of each sample ## 
    mixed.new <- data.frame(matrix(ncol = 9, nrow = 0))
    
    colnames(mixed.new) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                             "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
    ## for each sample, find % correctly idenfied by the model and add to dataframes ## 
    samples <- unique(all.mixed$sample_id)
    for(i in 1:length(samples)){
      temp.dat <- all.mixed[all.mixed$sample_id == samples[i],] #get all predictions of sample
      total <- nrow(temp.dat)
      correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
      percent <- (correct/total) * 100
      all.mixed[all.mixed$sample_id == samples[i],]$percent_correct <- percent
      new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],
                    temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], 
                    temp.dat$`Admixture % SAS`[1], temp.dat$Model[1],percent )
      mixed.new[nrow(mixed.new) + 1,]<-new.list
    }
    
    
    #### BENCHMARK ### 
    df.mixed.names <- paste0(mixed.dir,"model_",mixed.num,"/model_",mixed.num, mixed.anc,"_results.tsv") 
    mixed.myfiles <- lapply(df.mixed.names, readr::read_tsv)
    mixed.bm <- do.call(rbind,mixed.myfiles)
    
    all.mixed.bm <- merge(mixed.bm, samples.ancestry, by.x="sample_id", by.y="patient")
    all.mixed.bm$Model <- "Eur"
    all.mixed.bm$percent_correct <- NA
    all.mixed.bm$consensus_ancestry[all.mixed.bm$consensus_ancestry %in% admix] <- "admixed"
    mixed.new.bm <- data.frame(matrix(ncol = 9, nrow = 0))
    colnames(mixed.new.bm) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                                "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
    
    samples <- unique(all.mixed.bm$sample_id)
    for(i in 1:length(samples)){
      temp.dat <- all.mixed.bm[all.mixed.bm$sample_id == samples[i],]
      total <- nrow(temp.dat)
      correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
      percent <- (correct/total) * 100
      all.mixed.bm[all.mixed.bm$sample_id == samples[i],]$percent_correct <- percent
      new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],temp.dat$`Admixture % AMR`[1],
                    temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], temp.dat$`Admixture % SAS`[1], 
                    temp.dat$Model[1],percent)
      mixed.new.bm[nrow(mixed.new.bm) + 1,]<-new.list
    }
    mixed.new.bm$consensus_ancestry <- paste0(mixed.new.bm$consensus_ancestry, "BM")
    mixed.new$consensus_ancestry <- paste0(mixed.new$consensus_ancestry, "PF")
    to.return <- rbind(mixed.new.bm, mixed.new)
    
  }else{
    print("ENTER VALID ANCESTRY")
    
  }
}  


## validation set 
run <- function(disease, cur.dir, ancestry){

  cancer.type <- "BRCA"
  subtype1 <- "Basal"
  subtype2 <- "Luminal"
  eur.num <- list(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4), 
                  rep(10,4), rep(11,4), rep(12,4), rep(13,4), rep(14,4), rep(15,4), rep(16,4), rep(17,4))
  eur.num <- unlist(eur.num)
  afr.num <- list(rep(1,4), rep(2,4))
  afr.num <- unlist(afr.num)
  eas.num <- list(rep(1,3))
  eas.num <- unlist(eas.num)
  admixed.num <- list(rep(1,3))
  admixed.num <- unlist(admixed.num)
  mixed.num <- list(rep(1,4),rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4))
  mixed.num <- unlist(mixed.num)
  eur.anc <- list("eur", "afr", "eas", "admixed")#17 models 
  afr.anc <- list("eur", "afr", "eas", "admixed")#2 models 
  eas.anc <- list("eur", "afr", "admixed") #1 model, not enought eas to run eas model on eas samples 
  admixed.anc <- list("eur", "afr", "eas") #1 model 
  mixed.anc <- list("eur", "afr", "eas", "admixed") # 6 models 
    
  cur.dir <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinalVALIDATION/"
  pf.eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/eur/")
  eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/eur/")
  ## read in results for every ancestry  that eur model was tested on##
  pf.df.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num, eur.anc,"_results.tsv") 
  pf.eur.myfiles <- lapply(pf.df.eur.names, readr::read_tsv)
  ## this gives every models prediction of every sample tested no ## 
  pf.eur <- do.call(rbind,pf.eur.myfiles) 
  
  ## add ancestry percentages information to samples ##
  all.eur <- merge(pf.eur, samples.ancestry, by.x="sample_id", by.y="patient")
  all.eur$Model <- "PF.Eur"
  all.eur$percent_correct <- NA
  admix <- c("afr_admix","eas_admix","eur_admix", "admix")
  all.eur$consensus_ancestry[all.eur$consensus_ancestry %in% admix] <- "admixed"
  ## new dataframe for only keeping one instance of each sample ## 
  eur.new <- data.frame(matrix(ncol = 9, nrow = 0))
  
  colnames(eur.new) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                         "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
  ## for each sample, find % correctly idenfied by the model and add to dataframes ## 
  samples <- unique(all.eur$sample_id)
  for(i in 1:length(samples)){
    temp.dat <- all.eur[all.eur$sample_id == samples[i],] #get all predictions of sample
    total <- nrow(temp.dat)
    correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
    percent <- (correct/total) * 100
    all.eur[all.eur$sample_id == samples[i],]$percent_correct <- percent
    new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],
                  temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], 
                  temp.dat$`Admixture % SAS`[1], temp.dat$Model[1],percent )
    eur.new[nrow(eur.new) + 1,]<-new.list
  }
  
  
  #### BENCHMARK ### 
  df.eur.names <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num, eur.anc,"_results.tsv") 
  eur.myfiles <- lapply(df.eur.names, readr::read_tsv)
  eur.bm <- do.call(rbind,eur.myfiles)
  
  all.eur.bm <- merge(eur.bm, samples.ancestry, by.x="sample_id", by.y="patient")
  all.eur.bm$Model <- "Eur"
  all.eur.bm$percent_correct <- NA
  admix <- c("afr_admix","eas_admix","eur_admix", "admix")
  all.eur.bm$consensus_ancestry[all.eur.bm$consensus_ancestry %in% admix] <- "admixed"
  eur.new.bm <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(eur.new.bm) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                            "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
  
  samples <- unique(all.eur.bm$sample_id)
  for(i in 1:length(samples)){
    temp.dat <- all.eur.bm[all.eur.bm$sample_id == samples[i],]
    total <- nrow(temp.dat)
    correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
    percent <- (correct/total) * 100
    all.eur.bm[all.eur.bm$sample_id == samples[i],]$percent_correct <- percent
    new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],temp.dat$`Admixture % AMR`[1],
                  temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], temp.dat$`Admixture % SAS`[1], 
                  temp.dat$Model[1],percent)
    eur.new.bm[nrow(eur.new.bm) + 1,]<-new.list
  }
  eur.new.bm$consensus_ancestry <- paste0(eur.new.bm$consensus_ancestry, "BM")
  eur.new$consensus_ancestry <- paste0(eur.new$consensus_ancestry, "PF")
  to.return <- rbind(eur.new.bm, eur.new)
}



e+scale_fill_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#238B45", "#41AB5D","#74C476","#CE1256","#CA4136","#67000D", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#825283"))  +
scale_color_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#238B45", "#41AB5D","#74C476","#CE1256","#CA4136","#67000D", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#825283")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
