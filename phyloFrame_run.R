
library("optparse") # command line input 
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

source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/new_exomeAF_nosex_ordering.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/expression_elasticnet.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/phyloFrame_main.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/sample_selection.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/diseases/breast.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/diseases/thyroid.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/diseases/uterine.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/figure_scripts/scaled_importance_matrix.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/figure_scripts/get_ancestry_models.R")



### ARGUMENTS TO BE PASSED IN VIA COMMAND LINE 
# - tissue of disease 
# - penalty for model 

option_list <- list(
  make_option(c("-d", "--disease"), type="character", default=NULL, 
              help="Tissue diseases occurs in", metavar="character"),
  make_option(c("-m","--mixture"), type="numeric", default=0.5, 
              help="Mixture of LASSO and Ridge", metavar="numeric"),
  make_option(c("-p","--path"), type="character", default=NULL, 
              help="Model directory", metavar="character"),
  make_option(c("-V","--version"), type="numeric", default=1, 
              help="Model version", metavar="numeric"),
  make_option(c("-b", "--bench_penalty"),type="numeric", default=0.8, 
              help="Benchmark penalty for V1", metavar="numeric"),
  make_option(c("-g", "--variable_genes"),type="numeric", default=10000, 
              help="n variable genes", metavar="numeric")
 
); 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

disease <- opt$disease
en.mixture <- opt$mixture
model.path <- opt$path
mod.version <- opt$version
run.penalties <- opt$bench_penalty
variable.genes <- opt$variable_genes

print("############################### INPUT VARIABLES ###############################")
print(paste0("DISEASE: ", disease, " V2 AND BASE PF PENALTY: ", en.mixture, " MODEL PATH: ", model.path, " MODEL VERSION: ", mod.version, " V1 BENCHMARK PENALTY: ",run.penalties))

exomeAF <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/mean_enhancedAF_exome.tsv", col_names = TRUE)

## CALL INIT METHHOD FOR GIVEN CANCER ##
if (disease == "breast"){
  initialize.disease <- breast.init()
  cancer.type <- "BRCA"
  subtype1 <- "Basal"
  subtype2 <- "Luminal"
  continue <- 1
}else if(disease == "thyroid"){
  initialize.disease <- thyroid.init()
  cancer.type <- "THCA"
  subtype1 <- "M0"
  subtype2 <- "MX"
  continue <- 1
}else if(disease == "uterine"){
  initialize.disease <- uterine.init()
  cancer.type <- "UCEC"
  subtype1 <- "Endometrioid"
  subtype2 <- "Serous"
  continue <- 1
}else{
  print("Please enter a valid disease, currently they are: 1. breast 2. thyroid 3. uterine")
  continue <- 0
  stop()
}

if(continue == 1){
  expr.mat <- initialize.disease$expr
  clinical <- initialize.disease$clin
  samples.ancestry <- initialize.disease$samples.anc
  network <- initialize.disease$net
  
  ######################################################################
  #TEMPORARY ADD FOR THE VALIDAITON SET - TO RERUN VAL UNCOMMENT OUT 
  ######################################################################
  # validaiton set has new alias for some genes and excludes others from TCGA data 
  # this makes sure model only selects genes that are in the validation set
  # we just read this in to get the genes
  # expr.dat <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/AACR_2023_poster/GSE211167_tpm_matrix.txt")
  # -- replace gnees with alias in the validation set, then only keep common genes between the datasets 
  # names(expr.dat)[names(expr.dat) == "HID1"] <- "C17orf28"
  # names(expr.dat)[names(expr.dat) == "CYP2B7P"] <- "CYP2B7P1"
  # names(expr.dat)[names(expr.dat) == "MISP"] <- "C19orf21"
  # names(expr.dat)[names(expr.dat) == "BRINP3"] <- "FAM5C"
  # names(expr.dat)[names(expr.dat) == "CT83"] <- "CXorf61"
  # names(expr.dat)[names(expr.dat) == "BPIFB1"] <- "C20orf114"
  # names(expr.dat)[names(expr.dat) == "FDCSP"] <- "C4orf7"
  # names(expr.dat)[names(expr.dat) == "SOWAHA"] <- "ANKRD43"
  # names(expr.dat)[names(expr.dat) == "BRINP3"] <- "FAM5C"
  # names(expr.dat)[names(expr.dat) == "MS4A8"] <- "MS4A8B"
  # names(expr.dat)[names(expr.dat) == "NEURL1"] <- "NEURL"
  # names(expr.dat)[names(expr.dat) == "NSG2"] <- "HMP19"
  # names(expr.dat)[names(expr.dat) == "GPC1"] <- "GPC1-AS1"
  # 
  # allowed.genes <- colnames(expr.dat)
  
  ######################################################################
  
  ## trim matrix to delete NA columns, modifies samples and genes if provided in 2nd and 3rd param## 
  expr <- trim.expr.matrix(expr.mat, NULL, NULL) #for validation set change 2nd param to allowed.genes - then change back to NULL 
  
  # if genes have 2 sets of reads, take average 
  expr <- expr %>% 
    group_by(Hugo_Symbol) %>% 
    summarise_if(is.numeric, mean, na.rm = TRUE)

  if (disease == "breast"){
    expression <- add.subtype.breast(expr, clinical)

  }else if(disease == "thyroid"){
    expression <- add.subtype.thyroid(expr, clinical)
    
  }else if(disease == "uterine"){
    expression <- add.subtype.uterine(expr, clinical)
    
  }else{
    print("ERROR 1")
  }

  dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease, "/")######CHANGE PER CANCER TYPE
  
  ######################################################################
  #FOR CREATING NEW SAMPLE BATCHES IN THE RUN
  ######################################################################
 #  batch.dir <- paste0(dir, "samples/")
 #  batch.count <- create.batches(batch.dir, cancer.type, expression, subtype1, subtype2) ######CHANGE PER CANCER TYPE
  # afr.batch <- batch.count$afr
  # eas.batch <- batch.count$eas
  # eur.batch <- batch.count$eur
  # admixed.batch <- batch.count$admixed
  # mixed.batch <- batch.count$mixed
  # mixed.samples <- batch.count$mixed_samples
  #### 
  
  # for not making new samples each run. 
  afr.batch <- length(list.files(paste0(dir,"samples/afr/")))
  eas.batch <- length(list.files(paste0(dir,"samples/eas/")))
  eur.batch <- length(list.files(paste0(dir,"samples/eur/")))
  admixed.batch <- length(list.files(paste0(dir,"samples/admixed/")))
  mixed.batch <- length(list.files(paste0(dir,"samples/mixed_ancestry/")))
  ## to get all mixed sample IDs, this is to test on mixed samples, because mixed samples are interwoven 
  ## with all other ancestries, so we have to read in all batches ## 
  mixed <- list.files(paste0(dir,"samples/mixed_ancestry/"))
  mixed <- paste0(dir,"samples/mixed_ancestry/", mixed)
  all.mixed.samples <- lapply(mixed, readr::read_tsv)
  all.mixed <- do.call(rbind, all.mixed.samples)
  mixed.samples <- all.mixed$x
  
  ## seperate samples by ancestry ## 
  ## for passing samples to test models on ## 
  eur.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur",]$patient
  afr.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr",]$patient
  eas.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eas",]$patient
  admix.samples.eas <- samples.ancestry[samples.ancestry$consensus_ancestry == "eas_admix",]$patient
  admix.samples.afr <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr_admix",]$patient
  admix.samples.eur <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur_admix",]$patient
  admix.samples.gen <- samples.ancestry[samples.ancestry$consensus_ancestry == "admix",]$patient
  admixed.samples <- c(admix.samples.eas, admix.samples.afr, admix.samples.eur, admix.samples.gen)
  
  ## get ancestry samples expression matrices to run for testing 
  eur.expr <- trim.expr.matrix_OG(expression, eur.samples)
  eur.expr$subtype <- as.factor(eur.expr$subtype)
  afr.expr <- trim.expr.matrix_OG(expression, afr.samples)
  afr.expr$subtype <- as.factor(afr.expr$subtype)
  eas.expr <- trim.expr.matrix_OG(expression, eas.samples)
  eas.expr$subtype <- as.factor(eas.expr$subtype)
  admixed.expr <- trim.expr.matrix_OG(expression, admixed.samples)
  admixed.expr$subtype <- as.factor(admixed.expr$subtype)
  mixed.expr <- trim.expr.matrix_OG(expression, mixed.samples)
  mixed.expr$subtype <- as.factor(mixed.expr$subtype)

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

  #################################### EUROPEAN MODEL RUN ON ALL ANCESTRIES ######################################
  print("################################### STARTING EUROPEAN ANALYSIS ###################################")
  
  df.path.eur <- paste0(dir,"samples/eur/") #dir to sample batches
  dir.eur <- paste0(dir, model.path,"/model_runs/eur/") # dir to output
  pf.dir.eur <- paste0(dir, model.path,"/model_runs/phyloFrame/eur/")
  pf.base.eur <- paste0(dir, model.path,"/model_runs/phyloFrame/base/eur/")
  ancestry.dir <- paste0(dir, model.path,"/ancestry_genes/")
  
  for(i in 1:eur.batch){
    #create model specific dirs
    dir.create(paste0(dir.eur,"model_",i))
    new.dir.eur <- paste0(dir.eur,"model_",i)
    dir.create(paste0(pf.dir.eur,"model_",i))
    pf.new.dir.eur <- paste0(pf.dir.eur,"model_",i)
    dir.create(paste0(pf.base.eur,"model_",i))
    pf.base.new.eur <- paste0(pf.base.eur,"model_",i)
    #read in sample batch
    sample.num <- readr::read_tsv(paste0(df.path.eur, "samples",i,".tsv"), col_names = TRUE)
    # keep batch samples for training
    eur.train <- eur.expr[rownames(eur.expr) %in% sample.num$x,]
    table(eur.train$subtype)
    eur.train$subtype <- as.factor(eur.train$subtype)
    
    # put non-batch samples for testing 
    # make sure there were enough samples for a test set - otherwise mark 0 and skip testset for ancestry
    #EX admixed and afr rarely have enough models to train and have a seperate test set for themselves
    eur.test <- eur.expr[!(rownames(eur.expr) %in% sample.num$x),]
    if(nrow(eur.test) == 0){
      test.samples <- 0
    }else{
      test.samples <- 1
      eur.test$subtype <- as.factor(eur.test$subtype)
    }
    out_f <- paste0("model_",i)
    
    ## keep top 10000 genes for workflow
    base.eur.genes <- pf_top_varying_genes(eur.train, variable.genes)
    ## get the baseline disease signature - NOTE: phyloFrame function just runs elastic net on a subset of genes - the phyloFrame method is primary implemented in 
    # the function V2, where the network walk and allele frequency sorting occur. - TODO function is poorly named and needs to be be changed. 
    eur.base_pf <- phyloFrame(base.eur.genes, eur.train, pf.base.new.eur, out_f, en.mixture) 
    
    ## read back in the signature to use genes as start for network walk
    model.genes <- read.table(paste0(pf.base.new.eur, "/",out_f,"_all_sig.txt"))
    model.genes <- model.genes[-1,]
    colnames(model.genes) <- c("Variable", "Importance", "Sign")
    model.genes <- model.genes$Variable

    #get network genes for phyloframe and varying genes for benchmark for this batch
    model.input.genes <- get.genes.V2(network, model.genes, eur.train, exomeAF) 
    ancestry.genes <- model.input.genes$phyloFrame
    
    # for writing ancestry specific genes to file - just FYI not further used on model
    temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
    anc.genes.dat <- as.data.frame(temp.anc.genes)
  	gene.list <- base.eur.genes[-length(base.eur.genes)]
  	
  	# sample benchmark genes to match signature of phyloFrame 
  	index <- sample(1:length(gene.list), length(temp.anc.genes))  
  	benchmark.genes <- gene.list[index]
  	benchmark.genes <- c(benchmark.genes, model.genes,"subtype")
  	
    #write ancestry specific genes to file
    write.table(anc.genes.dat, file = paste0(ancestry.dir, "/eur_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
    eur_model <- phyloFrame(benchmark.genes, eur.train, new.dir.eur, out_f, run.penalties) #BENCHMARK - penalty here will always be 0
    pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")
    pf.eur.model <- phyloFrame(pf.pass.genes, eur.train, pf.new.dir.eur, out_f, run.penalties) #PHYLOFRAME
    
    ### run model on other ancestries ###
    afr.out <- paste0("model_",i,"afr")
    model.metrics(eur_model, afr.expr, new.dir.eur, afr.out, subtype1, subtype2)
    model.metrics(pf.eur.model, afr.expr, pf.new.dir.eur, afr.out, subtype1, subtype2)
    eas.out <- paste0("model_",i,"eas")
    model.metrics(eur_model, eas.expr, new.dir.eur, eas.out, subtype1, subtype2)
    model.metrics(pf.eur.model, eas.expr, pf.new.dir.eur, eas.out, subtype1, subtype2)
    if(test.samples == 1){ #make sure we actually have samples to test on for ancestry
      eur.out <- paste0("model_",i,"eur")
      model.metrics(eur_model, eur.test, new.dir.eur, eur.out, subtype1, subtype2)
      model.metrics(pf.eur.model, eur.test, pf.new.dir.eur, eur.out, subtype1, subtype2)
    }
    admixed.out <- paste0("model_",i,"admixed")
    model.metrics(eur_model, admixed.expr, new.dir.eur, admixed.out, subtype1, subtype2)
    model.metrics(pf.eur.model, admixed.expr, pf.new.dir.eur, admixed.out, subtype1, subtype2)
    mixed.out <- paste0("model_",i,"mixed")
    temp.mixed.expr <- mixed.expr[!(rownames(mixed.expr) %in% sample.num$x),]
    model.metrics(eur_model, temp.mixed.expr, new.dir.eur, mixed.out, subtype1, subtype2)
    model.metrics(pf.eur.model, temp.mixed.expr, pf.new.dir.eur, mixed.out, subtype1, subtype2)
    
  }
  
  #################################### AFRICAN MODEL RUN ON ALL ANCESTRIES ######################################
  print("################################### STARTING AFRICAN ANALYSIS ###################################")
  
  df.path.afr <- paste0(dir,"samples/afr/") #dir to sample batches
  dir.afr <- paste0(dir, model.path,"/model_runs/afr/") # dir to output
  pf.dir.afr <- paste0(dir, model.path,"/model_runs/phyloFrame/afr/")
  pf.base.afr <- paste0(dir, model.path, "/model_runs/phyloFrame/base/afr/")
  for(i in 1:afr.batch){
    dir.create(paste0(dir.afr,"model_",i))
    new.dir.afr <- paste0(dir.afr,"model_",i)
    dir.create(paste0(pf.dir.afr,"model_",i))
    pf.new.dir.afr <- paste0(pf.dir.afr,"model_",i)
    dir.create(paste0(pf.base.afr,"model_",i))
    pf.base.new.afr <- paste0(pf.base.afr,"model_",i)
    
    sample.num <- readr::read_tsv(paste0(df.path.afr, "samples",i,".tsv"), col_names = TRUE)
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
    #only use top 10000 most variable genes 
    base.afr.genes <- pf_top_varying_genes(afr.train, variable.genes)
    out_f <- paste0("model_",i)
    afr.base_pf <- phyloFrame(base.afr.genes, afr.train, pf.base.new.afr, out_f, en.mixture)
    ## read back in the signature to use genes as start for network walk
    model.genes <- read.table(paste0(pf.base.new.afr, "/",out_f,"_all_sig.txt"))
    model.genes <- model.genes[-1,]
    colnames(model.genes) <- c("Variable", "Importance", "Sign")
    model.genes <- model.genes$Variable
    ## get network and ordered genes for phyloFrame and most varying genes of same number for benchmark
    model.input.genes <- get.genes.V2(network, model.genes, afr.train, exomeAF)
    ancestry.genes <- model.input.genes$phyloFrame 
    temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
    anc.genes.dat <- as.data.frame(temp.anc.genes)
    gene.list <- base.afr.genes[-length(base.afr.genes)]
    index <- sample(1:length(gene.list), length(temp.anc.genes)) 
    benchmark.genes <- gene.list[index]
    benchmark.genes <- c(benchmark.genes, model.genes,"subtype")
    #this is just for knowin which genes got added from ancestry - not used further in the model 
    write.table(anc.genes.dat, file = paste0(ancestry.dir, "/afr_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
    #penalties for models should always be 0
    afr_model <- phyloFrame(benchmark.genes, afr.train, new.dir.afr, out_f, run.penalties) #BENCHMARK
    pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")
    pf.afr.model <- phyloFrame(pf.pass.genes, afr.train, pf.new.dir.afr, out_f, run.penalties) #PHYLOFRAME 
   
    ### run model on other ancestries ### 
    if(test.samples == 1){
      afr.out <- paste0("model_",i,"afr")
      model.metrics(afr_model, afr.test, new.dir.afr, afr.out, subtype1, subtype2)
      model.metrics(pf.afr.model, afr.test, pf.new.dir.afr, afr.out, subtype1, subtype2)
    }
    eas.out <- paste0("model_",i,"eas")
    model.metrics(afr_model, eas.expr, new.dir.afr, eas.out, subtype1, subtype2)
    model.metrics(pf.afr.model, eas.expr, pf.new.dir.afr, eas.out, subtype1, subtype2)
    eur.out <- paste0("model_",i,"eur")
    model.metrics(afr_model, eur.expr, new.dir.afr, eur.out, subtype1, subtype2)
    model.metrics(pf.afr.model, eur.expr, pf.new.dir.afr, eur.out, subtype1, subtype2)
    admixed.out <- paste0("model_",i,"admixed")
    model.metrics(afr_model, admixed.expr, new.dir.afr, admixed.out, subtype1, subtype2)
    model.metrics(pf.afr.model, admixed.expr, pf.new.dir.afr, admixed.out, subtype1, subtype2)
    mixed.out <- paste0("model_",i,"mixed")
    temp.mixed.expr <- mixed.expr[!(rownames(mixed.expr) %in% sample.num$x),]
    model.metrics(afr_model, temp.mixed.expr, new.dir.afr, mixed.out, subtype1, subtype2)
    model.metrics(pf.afr.model, temp.mixed.expr, pf.new.dir.afr, mixed.out, subtype1, subtype2)
  }
  
  #################################### EAS MODEL RUN ON ALL ANCESTRIES ######################################
  print("################################### STARTING EAST ASIAN ANALYSIS ###################################")
  ## note: uterine does not have enought samples to create model for east asian 
  if(disease != "uterine"){
    df.path.eas <- paste0(dir,"samples/eas/") #dir to sample batches
    dir.eas <- paste0(dir, model.path,"/model_runs/eas/") # dir to output
    pf.dir.eas <- paste0(dir, model.path,"/model_runs/phyloFrame/eas/")
    pf.base.eas <- paste0(dir, model.path,"/model_runs/phyloFrame/base/eas/")
    for(i in 1:eas.batch){
      dir.create(paste0(dir.eas,"model_",i))
      new.dir.eas <- paste0(dir.eas,"model_",i)
      dir.create(paste0(pf.dir.eas,"model_",i))
      pf.new.dir.eas <- paste0(pf.dir.eas,"model_",i)
      dir.create(paste0(pf.base.eas,"model_",i))
      pf.base.new.eas <- paste0(pf.base.eas,"model_",i)
      sample.num <- readr::read_tsv(paste0(df.path.eas, "samples",i,".tsv"), col_names = TRUE)
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
      out_f <- paste0("model_",i)
      # - get samples for this batch - #
      base.eas.genes <- pf_top_varying_genes(eas.train,variable.genes)
      eas.base_pf <- phyloFrame(base.eas.genes, eas.train, pf.base.new.eas, out_f, en.mixture) 
      ## read back in the signature to use genes as start for network walk
      model.genes <- read.table(paste0(pf.base.new.eas, "/",out_f,"_all_sig.txt"))
      model.genes <- model.genes[-1,]
      colnames(model.genes) <- c("Variable", "Importance", "Sign")
      model.genes <- model.genes$Variable
      model.input.genes <- get.genes.V2(network, model.genes, eas.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
      ancestry.genes <- model.input.genes$phyloFrame
      temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
      anc.genes.dat <- as.data.frame(temp.anc.genes)
      gene.list <- base.eas.genes[-length(base.eas.genes)]
      # sample genes to get the same amount of genes for benchmark sig
      index <- sample(1:length(gene.list), length(temp.anc.genes)) 
      benchmark.genes <- gene.list[index]
      benchmark.genes <- c(benchmark.genes, model.genes,"subtype")
      write.table(anc.genes.dat, file = paste0(ancestry.dir, "/eas_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
      eas.model <- phyloFrame(benchmark.genes, eas.train, new.dir.eas, out_f, run.penalties) #BENCHMARK
      pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")
      pf.eas.model <- phyloFrame(pf.pass.genes, eas.train, pf.new.dir.eas, out_f, run.penalties) #PHYLOFRAME
     
      ### run model on other ancestries ### 
      afr.out <- paste0("model_",i,"afr")
      model.metrics(eas.model, afr.expr, new.dir.eas, afr.out, subtype1, subtype2)
      model.metrics(pf.eas.model, afr.expr, pf.new.dir.eas, afr.out, subtype1, subtype2)
      if(test.samples == 1){
        eas.out <- paste0("model_",i,"eas")
        model.metrics(eas.model, eas.test, new.dir.eas, eas.out, subtype1, subtype2)
        model.metrics(pf.eas.model, eas.test, pf.new.dir.eas, eas.out, subtype1, subtype2)
      }
      eur.out <- paste0("model_",i,"eur")
      model.metrics(eas.model, eur.expr, new.dir.eas, eur.out, subtype1, subtype2)
      model.metrics(pf.eas.model, eur.expr, pf.new.dir.eas, eur.out, subtype1, subtype2)
      admixed.out <- paste0("model_",i,"admixed")
      model.metrics(eas.model, admixed.expr, new.dir.eas, admixed.out, subtype1, subtype2)
      model.metrics(pf.eas.model, admixed.expr, pf.new.dir.eas, admixed.out, subtype1, subtype2)
      mixed.out <- paste0("model_",i,"mixed")
      temp.mixed.expr <- mixed.expr[!(rownames(mixed.expr) %in% sample.num$x),]
      model.metrics(eas.model, temp.mixed.expr, new.dir.eas, mixed.out, subtype1, subtype2)
      model.metrics(pf.eas.model, temp.mixed.expr, pf.new.dir.eas, mixed.out, subtype1, subtype2)
    }
  }
  
  
  #################################### ADMIXED MODEL RUN ON ALL ANCESTRIES ######################################
  print("################################### STARTING ADMIXED ANALYSIS ###################################")
  
  df.path.adm <- paste0(dir,"samples/admixed/") #dir to sample batches
  dir.adm <- paste0(dir,model.path,"/model_runs/admixed/") # dir to output
  pf.dir.adm <- paste0(dir, model.path,"/model_runs/phyloFrame/admixed/")
  pf.base.adm <- paste0(dir, model.path,"/model_runs/phyloFrame/base/admixed/")
  for(i in 1:admixed.batch){
    dir.create(paste0(dir.adm,"model_",i))
    new.dir.adm <- paste0(dir.adm,"model_",i)
    dir.create(paste0(pf.dir.adm,"model_",i))
    pf.new.dir.adm <- paste0(pf.dir.adm,"model_",i)
    dir.create(paste0(pf.base.adm,"model_",i))
    pf.base.new.adm <- paste0(pf.base.adm,"model_",i)
    
    sample.num <- readr::read_tsv(paste0(df.path.adm, "samples",i,".tsv"), col_names = TRUE)
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
    out_f <- paste0("model_",i)
    base.admixed.genes <- pf_top_varying_genes(admixed.train, variable.genes)
    admixed.base_pf <- phyloFrame(base.admixed.genes, admixed.train, pf.base.new.adm, out_f, en.mixture)  
    ## read back in the signature to use genes as start for network walk
    model.genes <- read.table(paste0(pf.base.new.adm, "/",out_f,"_all_sig.txt"))
    model.genes <- model.genes[-1,]
    colnames(model.genes) <- c("Variable", "Importance", "Sign")
    model.genes <- model.genes$Variable
    #get network genes for phyloframe and varying genes for benchmark for this batch
    model.input.genes <- get.genes.V2(network, model.genes,  admixed.train, exomeAF) 
    ancestry.genes <- model.input.genes$phyloFrame # these get passed into phyloframe with no penalty
    temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
    anc.genes.dat <- as.data.frame(temp.anc.genes)
    gene.list <- base.admixed.genes[-length(base.admixed.genes)]
    index <- sample(1:length(gene.list), length(temp.anc.genes)) # need to add to each ancestry.
    benchmark.genes <- gene.list[index]
    benchmark.genes <- c(benchmark.genes, model.genes,"subtype")
    write.table(anc.genes.dat, file = paste0(ancestry.dir, "/admixed_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
    admixed.model <- phyloFrame(benchmark.genes, admixed.train, new.dir.adm, out_f, run.penalties) #BENCHMARK
    pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")
    pf.admixed.model <- phyloFrame(pf.pass.genes, admixed.train, pf.new.dir.adm, out_f, run.penalties) #PHYLOFRAME
    
    ### run model on other ancestries ## 
    afr.out <- paste0("model_",i,"afr")
    model.metrics(admixed.model, afr.expr, new.dir.adm, afr.out, subtype1, subtype2)
    model.metrics(pf.admixed.model, afr.expr, pf.new.dir.adm, afr.out, subtype1, subtype2)
    eas.out <- paste0("model_",i,"eas")
    model.metrics(admixed.model, eas.expr, new.dir.adm, eas.out, subtype1, subtype2)
    model.metrics(pf.admixed.model, eas.expr, pf.new.dir.adm, eas.out, subtype1, subtype2)
    eur.out <- paste0("model_",i,"eur")
    model.metrics(admixed.model, eur.expr, new.dir.adm, eur.out, subtype1, subtype2)
    model.metrics(pf.admixed.model, eur.expr, pf.new.dir.adm, eur.out, subtype1, subtype2)
    if(test.samples == 1){
      admixed.out <- paste0("model_",i,"admixed")
      model.metrics(admixed.model, admixed.test, new.dir.adm, admixed.out, subtype1, subtype2)
      model.metrics(pf.admixed.model, admixed.test, pf.new.dir.adm, admixed.out, subtype1, subtype2)
    }
    mixed.out <- paste0("model_",i,"mixed")
    temp.mixed.expr <- mixed.expr[!(rownames(mixed.expr) %in% sample.num$x),]
    model.metrics(admixed.model, temp.mixed.expr, new.dir.adm, mixed.out, subtype1, subtype2)
    model.metrics(pf.admixed.model, temp.mixed.expr, pf.new.dir.adm, mixed.out, subtype1, subtype2)
  }
  
  #################################### MIXED MODEL RUN ON ALL ANCESTRIES ######################################
  print("################################### STARTING MIXED ANALYSIS ###################################")
  
  df.path.mix <- paste0(dir,"samples/mixed_ancestry/") #dir to sample batches
  dir.mix <- paste0(dir,model.path,"/model_runs/mixed/") # dir to output
  pf.dir.mix <- paste0(dir,model.path,"/model_runs/phyloFrame/mixed/")
  pf.base.mix <- paste0(dir, model.path,"/model_runs/phyloFrame/base/mixed/")
  for(i in 1:mixed.batch){
    dir.create(paste0(dir.mix,"model_",i))
    new.dir.mix <- paste0(dir.mix,"model_",i)
    dir.create(paste0(pf.dir.mix,"model_",i))
    pf.new.dir.mix <- paste0(pf.dir.mix,"model_",i)
    dir.create(paste0(pf.base.mix,"model_",i))
    pf.base.new.mix <- paste0(pf.base.mix,"model_",i)
    sample.num <- readr::read_tsv(paste0(df.path.mix, "samples",i,".tsv"), col_names = TRUE)
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
    out_f <- paste0("model_",i)
    base.mixed.genes <- pf_top_varying_genes(mixed.train, variable.genes)
    mixed.base_pf <- phyloFrame(base.mixed.genes, mixed.train, pf.base.new.mix, out_f, en.mixture) ## PF BASELINE SIG - this gives back a model - we dont actually need it 
    ## read back in the signature to use genes as start for network walk
    model.genes <- read.table(paste0(pf.base.new.mix, "/",out_f,"_all_sig.txt"))
    model.genes <- model.genes[-1,]
    colnames(model.genes) <- c("Variable", "Importance", "Sign")
    model.genes <- model.genes$Variable
    ## get network and ordered genes for phyloFrame and most varying genes of same number for benchmark
    model.input.genes <- get.genes.V2(network, model.genes, mixed.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
    ancestry.genes <- model.input.genes$phyloFrame
    temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
    anc.genes.dat <- as.data.frame(temp.anc.genes)
    gene.list <- base.mixed.genes[-length(base.mixed.genes)]
    index <- sample(1:length(gene.list), length(temp.anc.genes)) # need to add to each ancestry.
    benchmark.genes <- gene.list[index]
    benchmark.genes <- c(benchmark.genes, model.genes,"subtype")
    
    write.table(anc.genes.dat, file = paste0(ancestry.dir, "/mixed_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
    mixed.model <- phyloFrame(benchmark.genes, mixed.train, new.dir.mix, out_f,run.penalties) #BENCHMARK
    pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")
    pf.mixed.model <- phyloFrame(pf.pass.genes, mixed.train, pf.new.dir.mix, out_f, run.penalties) #PHYLOFRAME
    
    ## run models on other ancestries ## 
    afr.out <- paste0("model_",i,"afr")
    temp.afr.expr <- afr.expr[!(rownames(afr.expr) %in% sample.num$x),]
    model.metrics(mixed.model, temp.afr.expr, new.dir.mix, afr.out, subtype1, subtype2)
    model.metrics(pf.mixed.model, temp.afr.expr, pf.new.dir.mix, afr.out, subtype1, subtype2)
    eas.out <- paste0("model_",i,"eas")
    temp.eas.expr <- eas.expr[!(rownames(eas.expr) %in% sample.num$x),]
    model.metrics(mixed.model, temp.eas.expr, new.dir.mix, eas.out, subtype1, subtype2)
    model.metrics(pf.mixed.model, temp.eas.expr, pf.new.dir.mix, eas.out, subtype1, subtype2)
    eur.out <- paste0("model_",i,"eur")
    temp.eur.expr <- eur.expr[!(rownames(eur.expr) %in% sample.num$x),]
    model.metrics(mixed.model, temp.eur.expr, new.dir.mix, eur.out, subtype1, subtype2)
    model.metrics(pf.mixed.model, temp.eur.expr, pf.new.dir.mix, eur.out, subtype1, subtype2)
    admixed.out <- paste0("model_",i,"admixed")
    temp.admixed.expr <- admixed.expr[!(rownames(admixed.expr) %in% sample.num$x),]
    model.metrics(mixed.model, temp.admixed.expr, new.dir.mix, admixed.out, subtype1, subtype2)
    model.metrics(pf.mixed.model, temp.admixed.expr, pf.new.dir.mix, admixed.out, subtype1, subtype2)
    if(test.samples == 1){
      mixed.out <- paste0("model_",i,"mixed")
      model.metrics(mixed.model, mixed.test, new.dir.mix, mixed.out, subtype1, subtype2)
      model.metrics(pf.mixed.model, mixed.test, pf.new.dir.mix, mixed.out, subtype1, subtype2)
    }
  }
}


##### CREATE SCATTER PLOTS AND IMPORTANCE HEATMAP#### 
run(disease, model.path, "admixed", "admixed", en.mixture)
run(disease, model.path, "mixed_ancestry", "mixed", en.mixture)
run(disease, model.path, "eur", "eur", en.mixture)
run(disease, model.path, "eas", "eas", en.mixture)
run(disease, model.path, "afr", "afr", en.mixture)


