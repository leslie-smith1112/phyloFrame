######## BENCHMARK european performance plots ########
#### NOTE - IF YOU WANT TO RUN AFR AND EAS YOU HAVE TO MAKE SURE THEY MATCH EUR. 
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breast.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/thyroid.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/uterine.R")

require(readr)  # for read_csv()
require(dplyr) #data frame handling 
library(tidyr)
library(stringi)
library(ggplot2)
####################################  DEFINE DISEASE AND MODEL NUMBERS FOR EACH ANCESTRY #################################### 
##TODO did this manually - change to read from directory## 
set.seed(123)
disease <- "thyroid"

####################################  read in ancestry samples from different batches (each ancestry has different # of batches) #################################### 
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
  eur.anc <- list("eur", "afr", "eas", "admixed")
  afr.anc <- list("eur", "afr", "eas", "admixed")
  eas.anc <- list("eur", "afr", "admixed")
}else if(disease == "thyroid"){
  cancer.type <- "THCA"
  subtype1 <- "M0"
  subtype2 <- "MX"
  eur.num <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
  eur.num <- unlist(eur.num)
  afr.num <- list(rep(1,3))
  afr.num <- unlist(afr.num)
  eas.num <- list(rep(1,4),rep(2,4), rep(3,4))
  eas.num <- unlist(eas.num)
  eur.anc <- list("eur", "afr", "eas", "admixed")
  afr.anc <- list("eur", "eas", "admixed")
  eas.anc <- list("eur", "afr", "eas", "admixed")
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
  #not ennough samples for uterine east asians for model 
  #uterine 
  eur.anc <- c("eur", "afr", "eas", "admixed")
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

## do once for phyloFrame then once for benchmark (below)
## PHYLOFRAME
#cur.dir <- "breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal" #DEFINE PATH TO FILES 
cur.dir <- "ph08TESTAF001_1_network02_051_neigh2_top30VarFinal"
#cur.dir <- "uterineph08TESTAF001_1_network02_051_neigh2_top30VarFinal"
#pf.eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/phyloFrame/eur/")

eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",cur.dir,"/model_runs/eur/")

# make disease capital for easy graphing title 


######################################################################## BENCHMARK ########################################################################

#################################### EURO MODEL #################################### 
# read in files 
eur.afr <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num,"afr_results.tsv") 
eur.eur <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num,"eur_results.tsv") 
eur.eas<- paste0(eur.dir,"model_",eur.num,"/model_",eur.num, "eas_results.tsv") 
eur.admixed <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num,"admixed_results.tsv") 

eur.afr.files <- lapply(eur.afr, readr::read_tsv)
eur.eur.files <- lapply(eur.eur, readr::read_tsv)
eur.eas.files <- lapply(eur.eas, readr::read_tsv)
eur.admixed.files <- lapply(eur.admixed, readr::read_tsv)

eur.afr.bm <- do.call(rbind,eur.afr.files)
eur.eur.bm <- do.call(rbind,eur.eur.files)
eur.eas.bm <- do.call(rbind,eur.eas.files)
eur.admixed.bm <- do.call(rbind,eur.admixed.files)

eur.afr.bm$ancestry <- "afr"
eur.eur.bm$ancestry <- "eur"
eur.eas.bm$ancestry <- "eas"
eur.admixed.bm$ancestry <- "admixed"

samples2 <- function()
{### all samples plotted twice ### 
  afr.dat <- data.frame(matrix(ncol = 3, nrow = 0))
  #taking the 
  
  for(i in 1:nrow(eur.afr.bm)){
    current <- eur.afr.bm[i,]
    m0val <- current$.pred_M0
    mxval <- current$.pred_MX
    
    minimal <- min(m0val, mxval)
    maxval <- max(m0val, mxval)
    
    temp <- c(current$sample_id,"afr_min",minimal)
    temp2 <- c(current$sample_id,"afr_max",maxval)
    
    afr.dat <- rbind(afr.dat, temp, temp2)
    
  }
  colnames(afr.dat) <- c("sample_id", "ancestry","value")
  
  eur.dat <- data.frame(matrix(ncol = 3, nrow = 0))
  
  for(i in 1:nrow(eur.eur.bm)){
    current <- eur.eur.bm[i,]
    m0val <- current$.pred_M0
    mxval <- current$.pred_MX
    
    minimal <- min(m0val, mxval)
    maxval <- max(m0val, mxval)
    
    temp <- c(current$sample_id,"eur_min",minimal)
    temp2 <- c(current$sample_id,"eur_max",maxval)
    
    eur.dat <- rbind(eur.dat, temp, temp2)
    
  }
  colnames(eur.dat) <- c("sample_id", "ancestry","value")
  
  eas.dat <- data.frame(matrix(ncol = 3, nrow = 0))
  
  for(i in 1:nrow(eur.eas.bm)){
    current <- eur.eas.bm[i,]
    m0val <- current$.pred_M0
    mxval <- current$.pred_MX
    
    minimal <- min(m0val, mxval)
    maxval <- max(m0val, mxval)
    
    temp <- c(current$sample_id,"eas_min",minimal)
    temp2 <- c(current$sample_id,"eas_max",maxval)
    
    eas.dat <- rbind(eas.dat, temp, temp2)
    
  }
  colnames(eas.dat) <- c("sample_id", "ancestry","value")
  dim(eas.dat)
  
  admixed.dat <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(admixed.dat) <- c("sample_id", "ancestry","value")
  for(i in 1:nrow(eur.admixed.bm)){
    current <- eur.admixed.bm[i,]
    m0val <- current$.pred_M0
    mxval <- current$.pred_MX
    
    minimal <- min(m0val, mxval)
    maxval <- max(m0val, mxval)
    
    temp <- c(current$sample_id,"admixed_min",minimal)
    temp2 <- c(current$sample_id,"admixed_max",maxval)
    
    admixed.dat <- rbind(admixed.dat, temp, temp2)
    
  }
  dim(admixed.dat)
  colnames(admixed.dat) <- c("sample_id", "ancestry","value")
  new.all <- rbind(afr.dat, eur.dat, admixed.dat, eas.dat)
  write.table(new.all, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/Eursample_max_min_prediction_supplement.tsv", sep = '/t', col.names = TRUE, row.names = FALSE)
  p <- ggplot(new.all, aes(x = ancestry, y = value, fill = ancestry, color = ancestry)) +
    geom_point()+ scale_fill_manual(values = c("#f0f71e","#FDC652","#1eb2f7","#1F619E","#76c20d",
                                                        "#496849","#f5a2a5","#f0f71e","#FDC652","#1eb2f7","#1F619E","#76c20d",
                                                        "#496849","#f5a2a5",
                                                        "#CA4136"))+ theme_minimal() +
                                                          theme(axis.text=element_text(size=23),axis.title=element_text(size=3,face="bold"))
  
  p
}

####### same as above but taking the average of each sample ######

all.anc <- rbind(eur.admixed.bm, eur.eur.bm, eur.eas.bm, eur.afr.bm)
samples <- unique(all.anc$sample_id) #get unique samples 
length(samples)
new.anc <- data.frame(matrix(ncol = 4, nrow = 0))
for(i in 1:nrow(all.anc)){
  temp.dat <- all.anc[all.anc$sample_id == samples[i],] #get all instances of samples 
  temp.dat$min <- pmin(temp.dat$.pred_M0, temp.dat$.pred_MX) #get the min of the row
  temp.dat$max <- pmax(temp.dat$.pred_M0, temp.dat$.pred_MX) # get max of row
  temp.dat$minavg <- mean(temp.dat$min) # take min and max mean seperatley to get one value for each per sample
  temp.dat$maxavg <- mean(temp.dat$max)
  new.row <- c(samples[i], temp.dat$minavg[1], temp.dat$maxavg[1], temp.dat$ancestry[1])
  new.anc <- rbind(new.anc, new.row)
}
colnames(new.anc) <- c("sample_id","minavg","maxavg","ancestry")
tt <- reshape2::melt(new.anc, id = c("sample_id","ancestry"))
head(tt)
p <- ggplot(tt, aes(x = ancestry, y = value, fill = ancestry, color = ancestry)) +
  geom_point()+ scale_fill_manual(values = c("#f0f71e","#FDC652","#1eb2f7","#1F619E","#76c20d",
                                                      "#496849","#f5a2a5","#f0f71e","#FDC652","#1eb2f7","#1F619E","#76c20d",
                                                      "#496849","#f5a2a5",
                                                      "#CA4136"))+ theme_minimal() +
                                                        theme(axis.text=element_text(size=23),axis.title=element_text(size=3,face="bold"))

p








p <- ggplot(all, aes(x=ancestry, y=benchmark, fill = ancestry, color = ancestry)) + geom_hline(yintercept=mean(all$benchmark)) +
  geom_beeswarm() + scale_fill_manual(values = c("#FDC652","#1F619E",
                                                          "#496849",
                                                          "#CA4136",
                                                          "#8D7260")) + scale_color_manual(values = c("#FDC652","#1F619E",
                                                          "#496849",
                                                          "#CA4136",
                                                          "#8D7260"))+ theme_minimal() +
                                                            theme(axis.text=element_text(size=23),axis.title=element_text(size=3,face="bold")) + ylim(0,1)



geom_hline(yintercept = mean(Indcatotvalue), color="blue")




samples <- unique(eur.afr.bm$sample_id) #
temp.dat <- eur.afr.bm[eur.afr.bm$sample_id == samples[i],]
temp.dat$min <- pmin(temp.dat$.pred_M0, temp.dat$.pred_MX)
temp.dat$max <- pmax(temp.dat$.pred_M0, temp.dat$.pred_MX)
temp.dat$minavg <- mean(temp.dat$min)
temp.dat$maxavg <- mean(temp.dat$max)

to.add <- c(temp.dat$sample_id[1], "admixed_min",temp.dat$minavg[1])
to.add <- c(temp.dat$sample_id[1], "admixed_max",temp.dat$maxavg[1])

total <- nrow(temp.dat)
correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
percent <- (correct/total) * 100
all.eas.bm[all.eas.bm$sample_id == samples[i],]$percent_correct <- percent
new.list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1],
              temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1], 
              temp.dat$`Admixture % SAS`[1], temp.dat$Model[1],percent )
eas.new.bm[nrow(eas.new.bm) + 1,]<-new.list















head(eur.bm)
eur.bm$max <- 
max(eur.bm[1,])


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

# eur.new.bm$max  <- pmax(eur.new.bm$`Admixture % AFR`, eur.new.bm$`Admixture % AMR`, eur.new.bm$`Admixture % EAS`,eur.new.bm$`Admixture % EUR`, eur.new.bm$`Admixture % SAS`)
# eur.new.bm$ancestry_percent <- eur.new.bm$consensus_ancestry

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
bm.eur.edit <- all.eur.bm[all.eur.bm$ancestry_percent == "afr" | all.eur.bm$ancestry_percent == "eur",]
###################### end ######################


ggplot(data=all.eur.bm , aes(x=admixture, y=percent_correct, group=ancestry_percent, color = ancestry_percent, fill = ancestry_percent))+  geom_smooth(method = loess, level = 0.3) +
  labs(title=paste0(disease," Benchmark Eur Model Performance on Admixture"),x="Admixture Percentage", y = "Model Identification") +
  scale_color_manual(values = c("afr" = "#1F619E",
                                "sas"="#E29901",
                                "eas"="#496849",
                                "eur" = "#CA4136",
                                "amr" = "#6C1BA9")) +
  scale_fill_manual(values = c("afr" = "#1F619E",
                               "sas"="#E29901",
                               "eas"="#496849",
                               "eur" = "#CA4136",
                               "amr" = "#6C1BA9")) +
  scale_y_continuous(limits=c(0, 100)) + geom_hline(yintercept=50,linetype=2) + theme_minimal()+
  theme(axis.text=element_text(size=23),axis.title=element_text(size=3,face="bold"))


#write.table(all.eur.bm,"/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/figures/admixture_BM_EUR_uterine.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


#################################### AFR MODEL ####################################
df.afr.names <- paste0(afr.dir,"model_",afr.num,"/model_",afr.num, afr.anc,"_results.tsv") 
afr.myfiles <- lapply(df.afr.names, readr::read_tsv)
afr.bm <- do.call(rbind,afr.myfiles)

all.afr.bm <- merge(afr.bm, samples.ancestry, by.x="sample_id", by.y="patient")
all.afr.bm$Model <- "Afr"
all.afr.bm$percent_correct <- NA

afr.new.bm <- data.frame(matrix(ncol = 9, nrow = 0))

colnames(afr.new.bm) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                          "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
all.afr.bm$consensus_ancestry <- substr(all.afr.bm$consensus_ancestry, 1, 3)
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
afr.new.bm$percent_correct <- as.numeric(afr.new.bm$percent_correct)
afr.new.bm$`Admixture % EAS` <- as.numeric(afr.new.bm$`Admixture % EAS`)
afr.new.bm$`Admixture % AFR` <- as.numeric(afr.new.bm$`Admixture % AFR`)
afr.new.bm$`Admixture % EUR` <- as.numeric(afr.new.bm$`Admixture % EUR`)
afr.new.bm$`Admixture % SAS` <- as.numeric(afr.new.bm$`Admixture % SAS`)
afr.new.bm$`Admixture % AMR` <- as.numeric(afr.new.bm$`Admixture % AMR`)

colnames.dat <- c("sample_id", "consensus_ancestry", "admixture", "model", 
                  "percent_correct", "ancestry_percent")

###### GET ROW FOR EACH ANCESTRY  #####
afr.afr <- afr.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % EAS`,
                                        -`Admixture % EUR`,-`Admixture % SAS`)
afr.afr$ancestry_slice <- "afr"
colnames(afr.afr) <- colnames.dat
## percent eas
afr.eas <- afr.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                        -`Admixture % EUR`,-`Admixture % SAS`)
afr.eas$ancestry_slice <- "eas"
colnames(afr.eas) <- colnames.dat
## percent afr
afr.eur <- afr.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                        -`Admixture % EAS`,-`Admixture % SAS`)
afr.eur$ancestry_slice <- "eur"
colnames(afr.eur) <- colnames.dat
## percent sas
afr.sas <- afr.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                        -`Admixture % EAS`,-`Admixture % EUR`)
afr.sas$ancestry_slice <- "sas"
colnames(afr.sas) <- colnames.dat
##percent amr
afr.amr <- afr.new.bm %>% dplyr::select(-`Admixture % SAS`, -`Admixture % AFR`,
                                        -`Admixture % EAS`,-`Admixture % EUR`)
afr.amr$ancestry_slice <- "amr"
colnames(afr.amr) <- colnames.dat

all.afr.bm <- rbind(afr.afr, afr.eas, afr.eur, afr.sas, afr.amr)
###################### add for only african and european percentages ######################
bm.afr.edit <- all.afr.bm[all.afr.bm$ancestry_percent == "afr" | all.afr.bm$ancestry_percent == "eur",]
###################### end ######################
ggplot(data=all.afr.bm, aes(x=admixture, y=percent_correct, group=ancestry_percent, color = ancestry_percent, fill = ancestry_percent))+ geom_smooth(method = loess, level = 0.8) +
  labs(title=paste0(disease," Benchmark Afr Model Performance on Admixture"),x="Admixture Percentage", y = "Model Identification") +
  scale_color_manual(values = c("afr" = "#1F619E",
                                "sas"="#E29901",
                                "eas"="#496849",
                                "eur" = "#CA4136",
                                "amr" = "#6C1BA9")) +
  scale_fill_manual(values = c("afr" = "#1F619E",
                               "sas"="#E29901",
                               "eas"="#496849",
                               "eur" = "#CA4136",
                               "amr" = "#6C1BA9")) +
  scale_y_continuous(limits=c(0, 100)) + geom_hline(yintercept=50,linetype=2) + theme_minimal()+
  theme(axis.text=element_text(size=23),axis.title=element_text(size=3,face="bold"))

#################################### EAS MODEL #################################### 
df.eas.names <- paste0(eas.dir,"model_",eas.num,"/model_",eas.num, eas.anc,"_results.tsv") 
eas.myfiles <- lapply(df.eas.names, readr::read_tsv)
eas.bm <- do.call(rbind,eas.myfiles)

all.eas.bm <- merge(eas.bm, samples.ancestry, by.x="sample_id", by.y="patient")
all.eas.bm$Model <- "Eas"
all.eas.bm$percent_correct <- NA

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

eas.new.bm$percent_correct <- as.numeric(eas.new.bm$percent_correct)
eas.new.bm$`Admixture % EAS` <- as.numeric(eas.new.bm$`Admixture % EAS`)
eas.new.bm$`Admixture % AFR` <- as.numeric(eas.new.bm$`Admixture % AFR`)
eas.new.bm$`Admixture % EUR` <- as.numeric(eas.new.bm$`Admixture % EUR`)
eas.new.bm$`Admixture % SAS` <- as.numeric(eas.new.bm$`Admixture % SAS`)
eas.new.bm$`Admixture % AMR` <- as.numeric(eas.new.bm$`Admixture % AMR`)

colnames.dat <- c("sample_id", "consensus_ancestry", "admixture", "model", 
                  "percent_correct", "ancestry_percent")

###### GET ROW FOR EACH ANCESTRY  #####
eas.afr <- eas.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % EAS`,
                                        -`Admixture % EUR`,-`Admixture % SAS`)
eas.afr$ancestry_slice <- "afr"
colnames(eas.afr) <- colnames.dat
## percent eas
eas.eas <- eas.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                        -`Admixture % EUR`,-`Admixture % SAS`)
eas.eas$ancestry_slice <- "eas"
colnames(eas.eas) <- colnames.dat
## percent eas
eas.eur <- eas.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                        -`Admixture % EAS`,-`Admixture % SAS`)
eas.eur$ancestry_slice <- "eur"
colnames(eas.eur) <- colnames.dat
## percent sas
eas.sas <- eas.new.bm %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,
                                        -`Admixture % EAS`,-`Admixture % EUR`)
eas.sas$ancestry_slice <- "sas"
colnames(eas.sas) <- colnames.dat
##percent amr
eas.amr <- eas.new.bm %>% dplyr::select(-`Admixture % SAS`, -`Admixture % AFR`,
                                        -`Admixture % EAS`,-`Admixture % EUR`)
eas.amr$ancestry_slice <- "amr"
colnames(eas.amr) <- colnames.dat

all.eas.bm <- rbind(eas.afr, eas.eas, eas.eur, eas.sas, eas.amr)

###################### add for only african and european percentages ######################
bm.eas.edit <- all.eas.bm[all.eas.bm$ancestry_percent == "afr" | all.eas.bm$ancestry_percent == "eur",] #"lm", "glm", "gam", "loess"
###################### end ######################
ggplot(data=all.eas.bm , aes(x=admixture, y=percent_correct, group=ancestry_percent, color = ancestry_percent, fill = ancestry_percent))+stat_smooth(method = loess, level = 0.8)+#,fullrange = TRUE,level = 0.2) +
  labs(title=paste0(disease," Benchmark Eas Model Performance on Admixture"),x="Admixture Percentage", y = "Model Identification") +
  scale_color_manual(values = c("afr" = "#1F619E",
                                "sas"="#E29901",
                                "eas"="#496849",
                                "eur" = "#CA4136",
                                "amr" = "#6C1BA9")) +
  scale_fill_manual(values = c("afr" = "#1F619E",
                               "sas"="#E29901",
                               "eas"="#496849",
                               "eur" = "#CA4136",
                               "amr" = "#6C1BA9")) +
  scale_y_continuous(limits=c(0, 100)) + geom_hline(yintercept=50,linetype=2) + theme_minimal()+
  theme(axis.text=element_text(size=23),axis.title=element_text(size=3,face="bold"))

# breaks <- seq(0, 1, by = 0.1)
# temp.binning <- all.eas.bm %>% mutate(admix_bin = cut(admixture, breaks=breaks))
# table(temp.binning$admix_bin)