
require(readr)  # for read_csv()
require(dplyr) #data frame handling 
library(tidyr)
library(stringi) #string mnutations 
library(tidymodels)
library(workflows) #tidy models package for bundling model specs 
library(parsnip) #for modeling 
library(workflows) #put model in workflow 
library(rpart.plot)
library(vip)
library(igraph)
set.seed(4831)
#phyloFrame(network, exomeAF, dz.signature, consensus_genes, expression, samples, clinical, out, "breast")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/network_grid_search/new_exomeAF_nosex_ordering.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/network_grid_search/expression_elasticnet.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breast.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/network_grid_search/phyloFrame_main.R")

exomeAF <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/mean_enhancedAF_exome.tsv", col_names = TRUE)
consensus_genes <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/consensus_genes.txt", col_names = TRUE)
consensus_genes <- consensus_genes$`Gene Symbol`

## 1. CALL INIT METHHOD FOR GIVEN CANCER ##
initialize.disease <- breast.init() ######CHANGE PER CANCER TYPE

expr.mat <- initialize.disease$expr
clinical <- initialize.disease$clin
samples.ancestry <- initialize.disease$samples.anc
network <- initialize.disease$net
dz.signature <- initialize.disease$dz_sig

### phyloframe base ###
network.genes <- get.network(network, dz.signature, 2, 0.75)
humanbase_network <- c(network.genes, consensus.genes) # could add PAM50 back here 
top.anc.dat <- order_frequencies(humanbase_network, exomeAF)
## -- add genes proven as cancer drivers -- ##
top.anc.dat <- c(top.anc.dat, consensus.genes, "subtype")


expr <- trim.expr.matrix(expr.mat, NULL, NULL)
expression <- add.subtype.breast(expr, clinical) ######CHANGE PER CANCER TYPE 
dir <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/"

batch.count <- create.batches(dir, "BRCA", expression, "Basal", "Luminal") ######CHANGE PER CANCER TYPE
subtype1 <- "Basal"
subtype2 <- "Luminal"

afr.batch <- batch.count$afr
eas.batch <- batch.count$eas
eur.batch <- batch.count$eur
admixed.batch <- batch.count$admixed
mixed.batch <- batch.count$mixed

##samples for testing after training ## 
eur.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur",]$patient
afr.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr",]$patient
eas.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eas",]$patient
admix.samples.eas <- samples.ancestry[samples.ancestry$consensus_ancestry == "eas_admix",]$patient
admix.samples.afr <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr_admix",]$patient
admix.samples.eur <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur_admix",]$patient
admix.samples.gen <- samples.ancestry[samples.ancestry$consensus_ancestry == "admix",]$patient
admixed.samples <- c(admix.samples.eas, admix.samples.afr, admix.samples.eur, admix.samples.gen)
mixed.samples <- c(eur.samples, afr.samples, eas.samples, admix.samples)

#2. get ancestry samples to run for testing 
eur.expr <- trim.expr.matrix_OG(expression, eur.samples)
afr.expr <- trim.expr.matrix_OG(expression, afr.samples)
eas.expr <- trim.expr.matrix_OG(expression, eas.samples)
admixed.expr <- trim.expr.matrix_OG(expression, admixed.samples)
mixed.expr <- trim.expr.matrix_OG(expression, mixed.samples)

# set ancestry in directories 
eur.dir <- paste0(dir,"eur/")
afr.dir <- paste0(dir,"afr/")
eas.dir <- paste0(dir, "eas/")
admixed.dir <- paste0(dir, "admixed/")
mixed.dir <- paste0(dir, "mixed_ancestry/")

# set ancestry out directories 
afr.out <- "afr_samples"
eur.out <- "eur_samples"
eas.out <- "eas_samples"
admixed.out <- "admixed_samples"
mixed.out <- "mixed_samples"
pf.out <- "phyloFrame"

#dir <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/model_runs/phyloFrame/mixed_matrix/"

##########################################################################################################################################################################################

#################################### EUROPEAN MODEL RUN ON ALL ANCESTRIES ######################################
df.path <- paste0(dir,"samples/eur/") #dir to sample batches
dir <- paste0(dir,"model_runs/eur/") # dir to output
pf.dir <- paste0(dir,"model_runs/phyloFrame/eur/")
for(i in 1:eur.batch){
  dir.create(paste0(dir,"model_",i))
  new.dir <- paste0(dir,"model_",i)
  dir.create(paste0(pf.dir,"model_",i))
  new.pf.dir <- paste0(pf.dir,"model_",i)
  
  sample.num <- readr::read_tsv(paste0(df.path, "samples",i,".tsv"), col_names = TRUE)
  #for each batch exclude batch sampels from test set 
  eur.train <- eur.expr[rownames(eur.expr) %in% sample.num$x,]
  eur.test <- eur.expr[!(rownames(eur.expr) %in% sample.num$x),]
  eur.train$subtype <- as.factor(eur.train$subtype)
  if(nrow(eur.test) == 0){
    test.samples <- 0
  }else{
    test.samples <- 1
    eur.test$subtype <- as.factor(eur.test$subtype)
  }
  # - get samples for this batch - #
  #eur.dat <- eur.clin[rownames(eur.clin) %in% sample.num$X1,]
  out_f <- paste0("model_",i)
  eur_model <- elasticnet.run(eur.train, new.dir, out_f)
  pf.eur.model <- phyloFrame(top.anc.dat, eur.train, new.dir, out_f)
  afr.out <- paste0("model_",i,"afr")
  model.metrics(eur_model, afr.expr, new.dir, afr.out, subtype1, subtype2)
  model.metrics(pf.eur.model, afr.expr, pf.new.dir, afr.out, subtype1, subtype2)
  eas.out <- paste0("model_",i,"eas")
  model.metrics(eur_model, eas.expr, new.dir, eas.out, subtype1, subtype2)
  model.metrics(pf.eur.model, eas.expr, pf.new.dir, eas.out, subtype1, subtype2)
  if(test.samples == 1){ #make sure we actually have samples to test on 
    eur.out <- paste0("model_",i,"eur")
    model.metrics(eur_model, eur.test, new.dir, eur.out, subtype1, subtype2)
    model.metrics(pf.eur.model, eur.test, pf.new.dir, eur.out, subtype1, subtype2)
  }
  admixed.out <- paste0("model_",i,"admixed")
  model.metrics(eur_model, admixed.expr, new.dir, admixed.out, subtype1, subtype2)
  model.metrics(pf.eur.model, admixed.expr, pf.new.dir, admixed.out, subtype1, subtype2)
  mixed.out <- paste0("model_",i,"mixed")
  model.metrics(eur_model, mixed.expr, new.dir, mixed.out, subtype1, subtype2)
  model.metrics(pf.eur.model, mixed.expr, pf.new.dir, mixed.out, subtype1, subtype2)
}

#################################### AFRICAN MODEL RUN ON ALL ANCESTRIES ######################################
df.path <- paste0(dir,"samples/afr/") #dir to sample batches
dir <- paste0(dir,"model_runs/afr/") # dir to output
pf.dir <- paste0(dir,"model_runs/phyloFrame/afr/")
for(i in 1:afr.batch){
  dir.create(paste0(dir,"model_",i))
  new.dir <- paste0(dir,"model_",i)
  dir.create(paste0(pf.dir,"model_",i))
  new.pf.dir <- paste0(pf.dir,"model_",i)
  sample.num <- readr::read_tsv(paste0(df.path, "samples",i,".tsv"), col_names = TRUE)
  #for each batch exclude batch sampels from test set 
  afr.train <- afr.expr[rownames(afr.expr) %in% sample.num$x,]
  afr.test <- afr.expr[!(rownames(afr.expr) %in% sample.num$x),]
  afr.train$subtype <- as.factor(afr.train$subtype)
  if(nrow(afr.test) == 0){
    test.samples <- 0
  }else{
    test.samples <- 1
    afr.test$subtype <- as.factor(afr.test$subtype)
  }
  # - get samples for this batch - #
  #eur.dat <- eur.clin[rownames(eur.clin) %in% sample.num$X1,]
  out_f <- paste0("model_",i)
  afr_model <- elasticnet.run(afr.train, new.dir, out_f)
  pf.afr.model <- phyloFrame(top.anc.dat, afr.train, new.dir, out_f)
  if(test.samples == 1){
    afr.out <- paste0("model_",i,"afr")
    model.metrics(afr_model, afr.test, new.dir, afr.out, subtype1, subtype2)
    model.metrics(pf.afr.model, afr.test, pf.new.dir, afr.out, subtype1, subtype2)
  }
  eas.out <- paste0("model_",i,"eas")
  model.metrics(afr_model, eas.expr, new.dir, eas.out, subtype1, subtype2)
  model.metrics(pf.afr.model, eas.expr, pf.new.dir, eas.out, subtype1, subtype2)
  eur.out <- paste0("model_",i,"eur")
  model.metrics(afr_model, eur.expr, new.dir, eur.out, subtype1, subtype2)
  model.metrics(pf.afr.model, eur.expr, pf.new.dir, eur.out, subtype1, subtype2)
  admixed.out <- paste0("model_",i,"admixed")
  model.metrics(afr_model, admixed.expr, new.dir, admixed.out, subtype1, subtype2)
  model.metrics(pf.afr.model, admixed.expr, pf.new.dir, admixed.out, subtype1, subtype2)
  mixed.out <- paste0("model_",i,"mixed")
  model.metrics(afr_model, mixed.expr, new.dir, mixed.out, subtype1, subtype2)
  model.metrics(pf.afr.model, mixed.expr, pf.new.dir, mixed.out, subtype1, subtype2)
}

#################################### EAS MODEL RUN ON ALL ANCESTRIES ######################################
df.path <- paste0(dir,"samples/eas/") #dir to sample batches
dir <- paste0(dir,"model_runs/eas/") # dir to output
pf.dir <- paste0(dir,"model_runs/phyloFrame/eas/")
for(i in 1:eas.batch){
  dir.create(paste0(dir,"model_",i))
  new.dir <- paste0(dir,"model_",i)
  dir.create(paste0(pf.dir,"model_",i))
  new.pf.dir <- paste0(pf.dir,"model_",i)
  sample.num <- readr::read_tsv(paste0(df.path, "samples",i,".tsv"), col_names = TRUE)
  #for each batch exclude batch sampels from test set 
  eas.train <- eas.expr[rownames(eas.expr) %in% sample.num$x,]
  eas.test <- eas.expr[!(rownames(eas.expr) %in% sample.num$x),]
  eas.train$subtype <- as.factor(eas.train$subtype)
  if(nrow(eas.test) == 0){
    test.samples <- 0
  }else{
    test.samples <- 1
    eas.test$subtype <- as.factor(eas.test$subtype)
  }
  
  # - get samples for this batch - #
  #eur.dat <- eur.clin[rownames(eur.clin) %in% sample.num$X1,]
  out_f <- paste0("model_",i)
  eas.model <- elasticnet.run(eas.train, new.dir, out_f)
  pf.eas.model <- phyloFrame(top.anc.dat, eas.train, new.dir, out_f)
  afr.out <- paste0("model_",i,"afr")
  model.metrics(eas.model, afr.expr, new.dir, afr.out, subtype1, subtype2)
  model.metrics(pf.eas.model, afr.expr, pf.new.dir, afr.out, subtype1, subtype2)
  if(test.samples == 1){
    eas.out <- paste0("model_",i,"eas")
    model.metrics(eas.model, eas.test, new.dir, eas.out, subtype1, subtype2)
    model.metrics(pf.eas.model, eas.test, pf.new.dir, eas.out, subtype1, subtype2)
  }
  eur.out <- paste0("model_",i,"eur")
  model.metrics(eas.model, eur.expr, new.dir, eur.out, subtype1, subtype2)
  model.metrics(pf.eas.model, eur.expr, pf.new.dir, eur.out, subtype1, subtype2)
  admixed.out <- paste0("model_",i,"admixed")
  model.metrics(eas.model, admixed.expr, new.dir, admixed.out, subtype1, subtype2)
  model.metrics(pf.eas.model, admixed.expr, pf.new.dir, admixed.out, subtype1, subtype2)
  mixed.out <- paste0("model_",i,"mixed")
  model.metrics(eas.model, mixed.expr, new.dir, mixed.out, subtype1, subtype2)
  model.metrics(pf.eas.model, mixed.expr, pf.new.dir, mixed.out, subtype1, subtype2)
}


#################################### ADMIXED MODEL RUN ON ALL ANCESTRIES ######################################
df.path <- paste0(dir,"samples/admixed/") #dir to sample batches
dir <- paste0(dir,"model_runs/admixed/") # dir to output
pf.dir <- paste0(dir,"model_runs/phyloFrame/admixed/")
for(i in 1:admixed.batch){
  dir.create(paste0(dir,"model_",i))
  new.dir <- paste0(dir,"model_",i)
  dir.create(paste0(pf.dir,"model_",i))
  new.pf.dir <- paste0(pf.dir,"model_",i)
  sample.num <- readr::read_tsv(paste0(df.path, "samples",i,".tsv"), col_names = TRUE)
  #for each batch exclude batch sampels from test set 
  
  admixed.train <- admixed.expr[rownames(admixed.expr) %in% sample.num$x,]
  admixed.test <- admixed.expr[!(rownames(admixed.expr) %in% sample.num$x),]
  admixed.train$subtype <- as.factor(admixed.train$subtype)
  if(nrow(admixed.test) == 0){
    test.samples <- 0
  }else{
    test.samples <- 1
    admixed.test$subtype <- as.factor(admixed.test$subtype)
  }
  # - get samples for this batch - #
  #eur.dat <- eur.clin[rownames(eur.clin) %in% sample.num$X1,]
  out_f <- paste0("model_",i)
  admixed.model <- elasticnet.run(admixed.train, new.dir, out_f)
  pf.admixed.model <- phyloFrame(top.anc.dat, admixed.train, new.dir, out_f)
  afr.out <- paste0("model_",i,"afr")
  model.metrics(admixed.model, afr.expr, new.dir, afr.out, subtype1, subtype2)
  model.metrics(pf.admixed.model, afr.expr, pf.new.dir, afr.out, subtype1, subtype2)
  eas.out <- paste0("model_",i,"eas")
  model.metrics(admixed.model, eas.expr, new.dir, eas.out, subtype1, subtype2)
  model.metrics(pf.admixed.model, eas.expr, pf.new.dir, eas.out, subtype1, subtype2)
  eur.out <- paste0("model_",i,"eur")
  model.metrics(admixed.model, eur.expr, new.dir, eur.out, subtype1, subtype2)
  model.metrics(pf.admixed.model, eur.expr, pf.new.dir, eur.out, subtype1, subtype2)
  if(test.samples == 1){
    admixed.out <- paste0("model_",i,"admixed")
    model.metrics(admixed.model, admixed.test, new.dir, admixed.out, subtype1, subtype2)
    model.metrics(pf.admixed.model, admixed.test, pf.new.dir, admixed.out, subtype1, subtype2)
  }
  mixed.out <- paste0("model_",i,"mixed")
  model.metrics(admixed.model, mixed.expr, new.dir, mixed.out, subtype1, subtype2)
  model.metrics(pf.admixed.model, mixed.expr, pf.new.dir, mixed.out, subtype1, subtype2)
}

#################################### MIXED MODEL RUN ON ALL ANCESTRIES ######################################
df.path <- paste0(dir,"samples/mixed/") #dir to sample batches
dir <- paste0(dir,"model_runs/mixed/") # dir to output
pf.dir <- paste0(dir,"model_runs/phyloFrame/mixed/")
for(i in 1:mixed.batch){
  dir.create(paste0(dir,"model_",i))
  new.dir <- paste0(dir,"model_",i)
  dir.create(paste0(pf.dir,"model_",i))
  new.pf.dir <- paste0(pf.dir,"model_",i)
  sample.num <- readr::read_tsv(paste0(df.path, "samples",i,".tsv"), col_names = TRUE)
  #for each batch exclude batch sampels from test set 
  
  mixed.train <- mixed.expr[rownames(mixed.expr) %in% sample.num$x,]
  mixed.test <- mixed.expr[!(rownames(mixed.expr) %in% sample.num$x),]
  mixed.train$subtype <- as.factor(mixed.train$subtype)
  if(nrow(mixed.test) == 0){
    test.samples <- 0
  }else{
    test.samples <- 1
    mixed.test$subtype <- as.factor(mixed.test$subtype)
  }
  # - get samples for this batch - #
  #eur.dat <- eur.clin[rownames(eur.clin) %in% sample.num$X1,]
  out_f <- paste0("model_",i)
  mixed.model <- elasticnet.run(mixed.train, new.dir, out_f)
  pf.mixed.model <- phyloFrame(top.anc.dat, mixed.train, new.dir, out_f)
  afr.out <- paste0("model_",i,"afr")
  model.metrics(mixed.model, afr.expr, new.dir, afr.out, subtype1, subtype2)
  model.metrics(pf.admixed.model, afr.expr, pf.new.dir, afr.out, subtype1, subtype2)
  eas.out <- paste0("model_",i,"eas")
  model.metrics(mixed.model, eas.expr, new.dir, eas.out, subtype1, subtype2)
  model.metrics(pf.admixed.model, eas.expr, pf.new.dir, eas.out, subtype1, subtype2)
  eur.out <- paste0("model_",i,"eur")
  model.metrics(mixed.model, eur.expr, new.dir, eur.out, subtype1, subtype2)
  model.metrics(pf.admixed.model, eur.expr, pf.new.dir, eur.out, subtype1, subtype2)
  admixed.out <- paste0("model_",i,"admixed")
  model.metrics(mixed.model, admixed.expr, new.dir, admixed.out, subtype1, subtype2)
  model.metrics(pf.admixed.model, admixed.expr, pf.new.dir, admixed.out, subtype1, subtype2)
  if(test.samples == 1){
    mixed.out <- paste0("model_",i,"mixed")
    model.metrics(mixed.model, mixed.test, new.dir, mixed.out, subtype1, subtype2)
    model.metrics(pf.admixed.model, mixed.expr, pf.new.dir, mixed.out, subtype1, subtype2)
  }
  
}
