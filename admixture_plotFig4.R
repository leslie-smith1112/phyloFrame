### - FIGURE 4 B-C - ### 

######## ADMIXTURE FIGURE ########

source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breast.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/thyroid.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/uterine.R")

require(readr)  # for read_csv()
require(dplyr) #data frame handling 
library(tidyr) 
library(stringi) #string matching 
library(ggplot2)
####################################  DEFINE DISEASE AND MODEL NUMBERS FOR EACH ANCESTRY #################################### 
##TODO did this manually - change to read from directory## 
set.seed(123)
disease <- "breast"

## PHYLOFRAME
####################################  read in ancestry samples from different batches (each ancestry has different # of batches) #################################### 
# list of ancestries defines the ancestries the model was able to be tested on (EX: eas could not be tested on itself in BRCA)
if (disease == "breast"){
  cur.dir <- "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal"
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
  eur.anc <- list("eur", "afr", "eas", "admixed")
  afr.anc <- list("eur", "afr", "eas", "admixed")
  eas.anc <- list("eur", "afr", "admixed")
}else if(disease == "thyroid"){
  cur.dir <- "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal"
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
  eur.anc <- list("eur", "afr", "eas", "admixed")
  afr.anc <- list("eur", "eas", "admixed")
  eas.anc <- list("eur", "afr", "eas", "admixed")
}else if(disease == "uterine"){
  cur.dir <- "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal"
  cancer.type <- "UCEC"
  subtype1 <- "Endometrioid"
  subtype2 <- "Serous"
  ## UTERINE
  eur.num <- list(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4), 
                  rep(10,4), rep(11,4), rep(12,4))
  eur.num <- unlist(eur.num)
  afr.num <- list(rep(1,4), rep(2,4))
  afr.num <- unlist(afr.num)
  #not ennough samples for uterine east asians for model 
  #uterine 
  eur.anc <- list("eur", "afr", "eas", "admixed")
  afr.anc <- list("eur", "afr", "eas", "admixed")
}else{
  print("Please enter a valid disease, currently they are: 1. breast 2. thyroid 3. uterine")
  continue <- 0
}

####################################  READ IN ANCESTRY INFORMATION #################################### 
estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == cancer.type,] 
samples.ancestry <- samples.ancestry[,-4:-8]
samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
#################################### DEFINE PATHS AND READ IN FILES ####################################
#PHYLOFRAME
pf.eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/eur/")
pf.afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/afr/")
pf.eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/eas/")
#BENCHMARK 
eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/eur/")
afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/afr/")
eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/eas/")

# make disease capital for easy graphing title 
if(disease == "thyroid"){
  disease <- "Thyroid"
}else if(disease == "breast"){
  disease <- "Breast"
}else if(disease == "uterine"){
  disease <- "Uterine"
}else {
  
}
######################################################################## PHYLOFRAME ########################################################################

#################################### EUR MODEL  ####################################
## read in results for every ancestry  that eur model was tested on##
pf.df.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num, eur.anc,"_results.tsv") 
pf.eur.myfiles <- lapply(pf.df.eur.names, readr::read_tsv)
## this gives every models prediction of every sample tested no ## 
pf.eur <- do.call(rbind,pf.eur.myfiles) 

## add ancestry percentages information to samples ##
all.eur <- merge(pf.eur, samples.ancestry, by.x="sample_id", by.y="patient")
all.eur$Model <- "PF.Eur"
all.eur$percent_correct <- NA

## new dataframe for only keeping one instance of each sample and their ancestry %s ## 
eur.new <- data.frame(matrix(ncol = 9, nrow = 0))

colnames(eur.new) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                       "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
all.eur$consensus_ancestry <- substr(all.eur$consensus_ancestry, 1, 3) # get rid of admix ending if there is one 
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

eur.new$percent_correct <- as.numeric(eur.new$percent_correct)
eur.new$`Admixture % EAS` <- as.numeric(eur.new$`Admixture % EAS`)
eur.new$`Admixture % AFR` <- as.numeric(eur.new$`Admixture % AFR`)
eur.new$`Admixture % EUR` <- as.numeric(eur.new$`Admixture % EUR`)
eur.new$`Admixture % SAS` <- as.numeric(eur.new$`Admixture % SAS`)
eur.new$`Admixture % AMR` <- as.numeric(eur.new$`Admixture % AMR`)

colnames.dat <- c("sample_id", "consensus_ancestry", "admixture", "model", 
                  "percent_correct", "ancestry_percent")
## seperate samples in to rows, so thet have one row for each ancestry % ## 
eur.afr <- eur.new %>% dplyr::select(-`Admixture % AMR`, -`Admixture % EAS`,
                                     -`Admixture % EUR`,-`Admixture % SAS`)
eur.afr$ancestry_slice <- "afr"
colnames(eur.afr) <- colnames.dat
## percent eas
eur.eas <- eur.new %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                     -`Admixture % EUR`,-`Admixture % SAS`)
eur.eas$ancestry_slice <- "eas"
colnames(eur.eas) <- colnames.dat
## percent eur
eur.eur <- eur.new %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                     -`Admixture % EAS`,-`Admixture % SAS`)
eur.eur$ancestry_slice <- "eur"
colnames(eur.eur) <- colnames.dat
## percent sas
eur.sas <- eur.new %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                     -`Admixture % EAS`,-`Admixture % EUR`)
eur.sas$ancestry_slice <- "sas"
colnames(eur.sas) <- colnames.dat
##percent amr
eur.amr <- eur.new %>% dplyr::select(-`Admixture % SAS`, -`Admixture % AFR`,
                                     -`Admixture % EAS`,-`Admixture % EUR`)
eur.amr$ancestry_slice <- "amr"
colnames(eur.amr) <- colnames.dat

pf.all.eur <- rbind(eur.afr, eur.eas, eur.eur, eur.sas, eur.amr)
library(tidyr)
pf.all.eur <- pf.all.eur %>% drop_na(admixture)


######################################################################## BENCHMARK ########################################################################

#################################### EURO MODEL #################################### 
df.eur.names <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num, eur.anc,"_results.tsv") 
eur.myfiles <- lapply(df.eur.names, readr::read_tsv)
eur.bm <- do.call(rbind,eur.myfiles)

all.eur.bm <- merge(eur.bm, samples.ancestry, by.x="sample_id", by.y="patient")
all.eur.bm$Model <- "Eur"
all.eur.bm$percent_correct <- NA

eur.new.bm <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(eur.new.bm) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                          "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
all.eur.bm$consensus_ancestry <- substr(all.eur.bm$consensus_ancestry, 1, 3)

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

eur.new.bm$percent_correct <- as.numeric(eur.new.bm$percent_correct)
eur.new.bm$`Admixture % EAS` <- as.numeric(eur.new.bm$`Admixture % EAS`)
eur.new.bm$`Admixture % AFR` <- as.numeric(eur.new.bm$`Admixture % AFR`)
eur.new.bm$`Admixture % EUR` <- as.numeric(eur.new.bm$`Admixture % EUR`)
eur.new.bm$`Admixture % SAS` <- as.numeric(eur.new.bm$`Admixture % SAS`)
eur.new.bm$`Admixture % AMR` <- as.numeric(eur.new.bm$`Admixture % AMR`)

colnames.dat <- c("sample_id", "consensus_ancestry", "admixture", "model", 
                  "percent_correct", "ancestry_percent")

###### GET ROW FOR EACH ANCESTRY  #####
eur.afr <- eur.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % EAS`,
                                        -`Admixture % EUR`,-`Admixture % SAS`)
eur.afr$ancestry_slice <- "afr"
colnames(eur.afr) <- colnames.dat
## percent eas
eur.eas <- eur.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                        -`Admixture % EUR`,-`Admixture % SAS`)
eur.eas$ancestry_slice <- "eas"
colnames(eur.eas) <- colnames.dat
## percent eur
eur.eur <- eur.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                        -`Admixture % EAS`,-`Admixture % SAS`)
eur.eur$ancestry_slice <- "eur"
colnames(eur.eur) <- colnames.dat
## percent sas
eur.sas <- eur.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                        -`Admixture % EAS`,-`Admixture % EUR`)
eur.sas$ancestry_slice <- "sas"
colnames(eur.sas) <- colnames.dat
##percent amr
eur.amr <- eur.new.bm %>% dplyr::select(-`Admixture % SAS`, -`Admixture % AFR`,
                                        -`Admixture % EAS`,-`Admixture % EUR`)
eur.amr$ancestry_slice <- "amr"
colnames(eur.amr) <- colnames.dat
all.eur.bm <- rbind(eur.afr, eur.eas, eur.eur, eur.sas, eur.amr)
all.eur.bm <- all.eur.bm %>% drop_na(admixture)

###################### add for only african and european percentages ######################

afr.dat.pf <- pf.all.eur[pf.all.eur$consensus_ancestry == "afr",]
afr.dat.bm <- all.eur.bm[all.eur.bm$consensus_ancestry == "afr",]
afr.dat.pf$model <- "phyloFrame"
afr.dat.bm$model <- "benchmark"

eur.dat.pf <- pf.all.eur[pf.all.eur$consensus_ancestry == "eur",]
eur.dat.bm <- all.eur.bm[all.eur.bm$consensus_ancestry == "eur",]
eur.dat.pf$model <- "phyloFrame"
eur.dat.bm$model <- "benchmark"

eas.dat.pf <- pf.all.eur[pf.all.eur$consensus_ancestry == "eas",]
eas.dat.bm <- all.eur.bm[all.eur.bm$consensus_ancestry == "eas",]
eas.dat.pf$model <- "phyloFrame"
eas.dat.bm$model <- "benchmark"

afr.dat  <- rbind(afr.dat.pf, afr.dat.bm)
eur.dat <- rbind(eur.dat.pf, eur.dat.bm)
eas.dat <- rbind(eas.dat.pf, eas.dat.bm)

afr.dat <- afr.dat[afr.dat$admixture >= .50,]
eur.dat <- eur.dat[eur.dat$admixture >= .50,]
eas.dat <- eas.dat[eas.dat$admixture >= .50,]

## NOTE: figures are flipped in paper for easy interpretation.
ggplot(data= eas.dat, aes(x=admixture, y=percent_correct, group=model, color = model, fill = model))+  geom_smooth(method = loess, level = 0.95) +geom_point() + scale_y_continuous(limits=c(0, 100)) +
  scale_x_continuous(limits=c(0.5, 1.00)) +
  geom_hline(yintercept=50,linetype=2) + theme_minimal() +
  scale_color_manual(values = c(
    "phyloFrame" = "#496849",
    "benchmark" = "#63b34c"
  )) +
  scale_fill_manual(values = c("phyloFrame" = "#496849",
                               "benchmark" = "#63b34c"))


ggplot(data= eur.dat, aes(x=admixture, y=percent_correct, group=model, color = model, fill = model))+  geom_smooth(method = loess, level = 0.95) +geom_point() + scale_y_continuous(limits=c(0, 100)) +
  scale_x_continuous(limits=c(0.5, 1.00)) +
  geom_hline(yintercept=50,linetype=2) + theme_minimal() +
  scale_color_manual(values = c(
    "phyloFrame" = "#CA4136",
    "benchmark" = "#f5a2a5"
  )) +
  scale_fill_manual(values = c("phyloFrame" = "#CA4136",
                               "benchmark" = "#f5a2a5"))


ggplot(data= afr.dat, aes(x=admixture, y=percent_correct, group=model, color = model, fill = model))+  geom_smooth(method = loess, level = 0.95) +geom_point() + scale_y_continuous(limits=c(0, 100)) +
  scale_x_continuous(limits=c(0.5, 1.00)) +
  geom_hline(yintercept=50,linetype=2) + theme_minimal() +
  scale_color_manual(values = c("phyloFrame" = "#1F619E",
                                "benchmark" = "#9fccfa")) +
  scale_fill_manual(values = c("phyloFrame" = "#1F619E",
                               "benchmark" = "#9fccfa"))

##################################################
