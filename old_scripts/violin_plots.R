### SUPP FIG 5 ###

## function reads in model performance for each ancestry and plot 

run <- function(disease, dir, ancestry.sample, ancestry.model, mixture,metric){
  temp.list <- list.files(paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/samples/",ancestry.sample,"/"))
  num <- length(temp.list)
  num
  mod.num <- (1:num)
  df.path <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",dir,"/model_runs/",ancestry.model)
  
  all <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(all) <- c("model.num", "ancestry", "benchmark")
 
  df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"afr_metrics.tsv") #CHANGE HERE
  myfiles <- lapply(df.names, readr::read_tsv)

  # - for benchmark model - #
  afr.mod <- do.call("rbind", myfiles)
  afr.mod <- afr.mod[afr.mod$.metric == metric,]
  
  model.num <- paste0("model_",mod.num)
  ancestry <- rep("afr", num) 
  new.afr.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                           "benchmark" = afr.mod$.estimate)
  all <- rbind(all, new.afr.df)
  
  ###ADMIXED 
  df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"admixed_metrics.tsv") #CHANGE HERE
  myfiles <- lapply(df.names, readr::read_tsv)

  #- for benchmark model - #
  admixed.mod <- do.call("rbind", myfiles)
  admixed.mod <- admixed.mod[admixed.mod$.metric == metric,]

  model.num <- paste0("model_",mod.num)
  ancestry <- rep("admixed", num) 
  new.admixed.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                               "benchmark" = admixed.mod$.estimate)
  all <- rbind(all, new.admixed.df)

  ####EAS 
  df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eas_metrics.tsv") #CHANGE HERE

  myfiles <- lapply(df.names, readr::read_tsv)
  #- for benchmark model - #
  eas.mod <- do.call("rbind", myfiles)
  eas.mod <- eas.mod[eas.mod$.metric == metric,]
  
  model.num <- paste0("model_",mod.num)
  ancestry <- rep("eas", num)  
  new.eas.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                           "benchmark" = eas.mod$.estimate)
  all <- rbind(all, new.eas.df)
  
  #### EUR 
  df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eur_metrics.tsv") #CHANGE HERE

  myfiles <- lapply(df.names, readr::read_tsv)

  #- for benchmark model - #
  eur.mod <- do.call("rbind", myfiles)
  eur.mod <- eur.mod[eur.mod$.metric == metric,]

  model.num <- paste0("model_",mod.num)
  ancestry <- rep("eur", num)  
  new.eur.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                           "benchmark" = eur.mod$.estimate)
  all <- rbind(all, new.eur.df)

 return(all)
  
}  
all <- all[all$ancestry != "mixed"]

all <- run("thyroid", "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eur","eur",0.8, "recall")

all <-run("breast", "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eur","eur",0.8)

all <- run("uterine", "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal", "eur","eur",0.8, "roc_auc")


library(ggbeeswarm)
p <- ggplot(all, aes(x=ancestry, y=benchmark, fill = ancestry, color = ancestry)) + geom_hline(yintercept=mean(all$benchmark)) +
  geom_beeswarm() + scale_fill_manual(values = c("#FDC652","#1F619E",
                                                        "#496849",
                                                        "#CA4136",
                                                        "#8D7260")) + scale_color_manual(values = c("#FDC652","#1F619E",
                                                        "#496849",
                                                        "#CA4136",
                                                        "#8D7260"))+ theme_minimal() +
                                                          theme(axis.text=element_text(size=23),axis.title=element_text(size=3,face="bold")) + ylim(0,1)

p


#### Leslie please make a violin plot showing admixed vs not admixed (95%+ one ancestry) 
#and calculate a t-test to show difference between benchmark and phyloframe. 
#Do this across all EUR models, grouped together.

#rephrase: for each eur model across ancestries, use violin plot to show the performance of admixed vs not admixed. 
cancer.type <- c("BRCA","UCEC","THCA")
####################################  READ IN ANCESTRY INFORMATION #################################### 
estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type %in% cancer.type,] 
samples.ancestry <- samples.ancestry[,-4:-8]
samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
#################################### DEFINE PATHS AND READ IN FILES ####################################
# set directories 
pf.thyroid <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/ph08TESTAF001_1_network02_051_neigh2_top30VarFinal/model_runs/phyloFrame/eur/")
bm.thyroid<- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/ph08TESTAF001_1_network02_051_neigh2_top30VarFinal/model_runs/eur/")

pf.breast <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal/model_runs/phyloFrame/eur/")
bm.breast<- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal/model_runs/eur/")

pf.uterine <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal/model_runs/phyloFrame/eur/")
bm.uterine<- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal/model_runs/eur/")

#
ancestries <- c("afr", "admixed", "eas","eur")
### define model numbers for each ancestry ### 
head(samples.ancestry)
cancer.type <- "BRCA"
subtype1 <- "Basal"
subtype2 <- "Luminal"
eur.breast <- list(rep(1,5), rep(2,5), rep(3,5), rep(4,5), rep(5,5), rep(6,5), rep(7,5), rep(8,5), rep(9,5), 
                rep(10,5), rep(11,5), rep(12,5), rep(13,5), rep(14,5), rep(15,5), rep(16,5), rep(17,5))
eur.breast <- unlist(eur.breast)

cancer.type <- "THCA"
subtype1 <- "M0"
subtype2 <- "MX"
eur.thyroid <- list(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4), 
                    rep(10,4), rep(11,4), rep(12,4), rep(13,4), rep(14,4), rep(14,4), rep(16,4), rep(17,4), rep(18,4), rep(19,4), rep(20,4), rep(21,4), rep(22,4), rep(23,4))
eur.thyroid <- unlist(eur.thyroid)

cancer.type <- "UCEC"
subtype1 <- "Endometrioid"
subtype2 <- "Serous"
## UTERINE
eur.uterine <- list(rep(1,5), rep(2,5), rep(3,5), rep(4,5), rep(5,5), rep(6,5), rep(7,5), rep(8,5), rep(9,5), 
                rep(10,5), rep(11,5), rep(12,5))
eur.uterine <- unlist(eur.uterine)

#function to read in all eur files for each ancestru for a disease
#returns data frame of all phyloframe runs and another of all benchmark runs. 
read.myfiles <- function(bm.dir, pf.dir, number,ancestry.list){
  pf <- paste0(pf.dir,"model_",number,"/model_",number, ancestry.list,"_results.tsv") 
  pf.files <- lapply(pf, readr::read_tsv)
  
  bm <- paste0(bm.dir,"model_",number,"/model_",number, ancestry.list,"_results.tsv") 
  bm.files <- lapply(bm, readr::read_tsv)
  
  pf.all <- do.call(rbind,pf.files)
  bm.all <- do.call(rbind,bm.files)#this returns all predictions from all models of this disease
  
  
  files.list <- list("benchmark" = bm.all, "phyloFrame" = pf.all)
  return(files.list)
}

performance <- function(disease.dat, benchmark){
  return.dat <- data.frame(matrix(ncol = 7, nrow = 0))
  # get the samples in the disease
  new.samples <- samples.ancestry[samples.ancestry$patient %in% disease.dat$sample_id,]
  new.samples$max <- pmax(new.samples$`Admixture % AFR`, new.samples$`Admixture % AMR`, new.samples$`Admixture % EAS`,
                          new.samples$`Admixture % EUR`, new.samples$`Admixture % SAS`)
  new.samples$admixed <- "temp"
  new.samples$admixed[new.samples$max >= 0.95] <- "not_admixed"
  new.samples$admixed[new.samples$max < 0.95] <- "admixed"
  #new.samples <- new.samples[new.samples$admixed != "temp",] ## come back to this 
  # now we have the admixed info, calculate the % correct 
  
  merged <- merge(new.samples, disease.dat, by.x = "patient", by.y = "sample_id")
  merged$percent_correct <- 0
  samples <- unique(merged$patient)
  for(i in 1:length(samples)){
    
    temp.dat <- merged[merged$patient == samples[i],] #get all predictions of sample
    total <- nrow(temp.dat) # get total predictions 
    correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
    percent <- (correct/total) * 100
    temp.dat$percent_correct[temp.dat$patient == samples[i]] <- percent
    temp.dat$model <- "temp"
    if(benchmark){
      temp.dat$model <- "benchmark"
    }else{
      temp.dat$model <- "phyloFrame"
    }
    new.list <- c(temp.dat$patient[1], temp.dat$tumor_type[1],temp.dat$max[1],temp.dat$admixed[1],temp.dat$percent_correct[1], temp.dat$model[1], temp.dat$consensus_ancestry[1])
    #return.dat <- rbind(return.dat, new.list)
    return.dat[nrow(return.dat) + 1,]<-new.list
  }
  colnames(return.dat) <- c("sample_id","disease","max_admix","admixed","percent_correct","model","consensus_ancestry")
  # return.dat$admixed[return.dat$consensus_ancestry == "admix"] <- "admixed"
  # return.dat$admixed[return.dat$admixed == "temp"] <- "not_admixed"
  return.dat <- return.dat[!(return.dat$admixed == "temp"),]
  return(return.dat)
}

thyroid.files <-read.myfiles(bm.thyroid, pf.thyroid, eur.thyroid, ancestries)
thy.bm <- thyroid.files$benchmark
thy.pf <- thyroid.files$phyloFrame
thy.dat <- performance(thy.bm, TRUE)
thy.dat.pf <- performance(thy.pf, FALSE)
dim(thy.dat)

library(ggplot2)
p <- ggplot(thy.dat, aes(factor(admixed, levels = c("admixed", "not_admixed")), percent_correct)) +  geom_violin() 
p

ggplot(thy.dat, aes(x = admixed, y = percent_correct, fill = admixed)) +
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none")












#, fill = admixed, color = admixed)) + 
  #geom_violin() 
# + scale_fill_manual(values = c("#FDC652","#1F619E",
#                                                         "#496849",
#                                                         "#CA4136",
#                                                         "#8D7260")) + scale_color_manual(values = c("#FDC652","#1F619E",
#                                                         "#496849",
#                                                         "#CA4136",
#                                                         "#8D7260"))+ theme_minimal() 

p


p <- ggplot(mtcars, aes(factor(cyl), mpg))
p + geom_violin()
yy <- thy.dat
yy$admixed[yy$consensus_ancestry == "admix"] <- "admixed"
yy$admixed[yy$admixed == "temp"] <- "not_admixed"


ty <- thy.dat[is.na(thy.dat$max_admix),]



breast.files <-read.myfiles(bm.breast, pf.breast, eur.breast, ancestries)
bre.bm <- breast.files$benchmark
bre.pf <- breast.files$phyloFrame
colnames(bre.bm) <- c("sample_id", "subtype",".pred_class", "pred_sub1", "pred_sub2")
colnames(bre.pf) <- c("sample_id", "subtype",".pred_class", "pred_sub1", "pred_sub2")

uterine.files <-read.myfiles(bm.uterine, pf.uterine, eur.uterine, ancestries)
ute.bm <- uterine.files$benchmark
ute.pf <- uterine.files$phyloFrame
colnames(ute.bm) <- c("sample_id", "subtype",".pred_class", "pred_sub1", "pred_sub2")
colnames(ute.pf) <- c("sample_id", "subtype",".pred_class", "pred_sub1", "pred_sub2")
dim(ute.pf)

all.pf <- rbind(thy.pf, bre.pf, ute.pf)
all.bm <- rbind(thy.bm, bre.bm, ute.bm)




## stack all diseases 


## add ancestry percentages information to samples ##
all.eur <- merge(pf.eur, samples.ancestry, by.x="sample_id", by.y="patient")
all.eur$Model <- "PF.Eur"
all.eur$percent_correct <- NA

## new dataframe for only keeping one instance of each sample ## 
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




# make plot
p <- ggplot(all, aes(x=ancestry, y=phyloFrame, fill = ancestry, color = ancestry)) + 
  geom_violin() + scale_fill_manual(values = c("#FDC652","#1F619E",
                                                        "#496849",
                                                        "#CA4136",
                                                        "#8D7260")) + scale_color_manual(values = c("#FDC652","#1F619E",
                                                        "#496849",
                                                        "#CA4136",
                                                        "#8D7260"))+ theme_minimal() 
p



























