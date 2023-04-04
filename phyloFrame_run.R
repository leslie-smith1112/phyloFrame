
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
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/differential_expression_GSEA_analysis.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/scaled_importance_matrix.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/get_ancestry_models.R")



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
  # make_option(c("-b", "--bench_penalty"),type="numeric", default=0.8, 
  #           help="Benchmark penalty for V1", metavar="numeric")),
  
  # ,
  # make_option(c("-n","--node"), type="numeric", default=2, 
  #             help="Neighbors to include in network", metavar="character"),
  # make_option(c("-e","--edge"), type="numeric", default=0.75, 
  #             help="Edge weight of genes to include in network]", metavar="character")
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

## 1. CALL INIT METHHOD FOR GIVEN CANCER ##

 ### TO DO FOR EAACH DISEASE ADD EDGE AND NODE VALUES 
if (disease == "breast"){
  initialize.disease <- breast.init()
  cancer.type <- "BRCA"
  subtype1 <- "Basal"
  subtype2 <- "Luminal"
  node <- 2 #1 # defined by grid search 
  edge <- 0.65#0.85 # defined by grid search
  continue <- 1
}else if(disease == "thyroid"){
  initialize.disease <- thyroid.init()
  cancer.type <- "THCA"
  subtype1 <- "M0"
  subtype2 <- "MX"
  node <- 2 ##2#1 #3
  edge <- 0.75#0.75# 0.65 #0.85
  continue <- 1
}else if(disease == "uterine"){
  initialize.disease <- uterine.init()
  cancer.type <- "UCEC"
  subtype1 <- "Endometrioid"
  subtype2 <- "Serous"
  node <- 2 #3
  edge <-  0.65  # 0.5 0.65
  continue <- 1
}else{
  print("Please enter a valid disease, currently they are: 1. breast 2. thyroid 3. uterine")
  continue <- 0
}

if(continue == 1){
  expr.mat <- initialize.disease$expr
  clinical <- initialize.disease$clin
  samples.ancestry <- initialize.disease$samples.anc
  network <- initialize.disease$net
  
  ## trim matrix to delete NA columns ## 
  expr <- trim.expr.matrix(expr.mat, NULL, NULL)
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
  

  ## differential expression analysis ##  - CHANGING HERE TO SEE WHAT HAPPEND IF GIVE BOTH VARIANCE 
  # de.expr <- column_to_rownames(expr, "Hugo_Symbol")  #DE needs samples in columns and genes in rows
  # de.expr <- de.expr[,colnames(de.expr) %in% rownames(expression)]
  # de.meta <- data.frame(rownames(expression), expression$subtype)
  # colnames(de.meta) <- c("sample_id", "subtype")
  # de.meta$subtype <- as.factor(de.meta$subtype)
  # de.expr <- de.expr %>%
  #   dplyr::select(de.meta$sample_id)
  # all.equal(colnames(de.expr), de.meta$sample_id)
  # de.genes <- compute_DE(de.expr, de.meta)
  # de.sig <- de.genes[de.genes$threshold == TRUE,]
  # dim(de.sig)
  # dz.genes <- de.sig$Gene

  dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease, "/")######CHANGE PER CANCER TYPE
  
  #### FOR CREATING NEW SAMPLE BATCHES IN THE RUN #### 
 #  batch.dir <- paste0(dir, "samples/")
 #  batch.count <- create.batches(batch.dir, cancer.type, expression, subtype1, subtype2) ######CHANGE PER CANCER TYPE
  # afr.batch <- batch.count$afr
  # eas.batch <- batch.count$eas
  # eur.batch <- batch.count$eur
  # admixed.batch <- batch.count$admixed
  # mixed.batch <- batch.count$mixed
  # mixed.samples <- batch.count$mixed_samples
  
  # for not making new samples each run. 
  afr.batch <- length(list.files(paste0(dir,"samples/afr/")))
  eas.batch <- length(list.files(paste0(dir,"samples/eas/")))
  eur.batch <- length(list.files(paste0(dir,"samples/eur/")))
  admixed.batch <- length(list.files(paste0(dir,"samples/admixed/")))
  mixed.batch <- length(list.files(paste0(dir,"samples/mixed_ancestry/")))
  ## to get all mixed sample IDs, this is to test on mixed samples, becasue mixed samples are interwovem with all othe ancestries ## 
  mixed <- list.files(paste0(dir,"samples/mixed_ancestry/"))
  mixed <- paste0(dir,"samples/mixed_ancestry/", mixed)
  all.mixed.samples <- lapply(mixed, readr::read_tsv)
  all.mixed <- do.call(rbind, all.mixed.samples)
  mixed.samples <- all.mixed$x
  
  ##seperate samples by ancestry ## 
  ## for passing samples to test models on ## 
  eur.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur",]$patient
  afr.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr",]$patient
  eas.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eas",]$patient
  admix.samples.eas <- samples.ancestry[samples.ancestry$consensus_ancestry == "eas_admix",]$patient
  admix.samples.afr <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr_admix",]$patient
  admix.samples.eur <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur_admix",]$patient
  admix.samples.gen <- samples.ancestry[samples.ancestry$consensus_ancestry == "admix",]$patient
  admixed.samples <- c(admix.samples.eas, admix.samples.afr, admix.samples.eur, admix.samples.gen)
  
  #2. get ancestry samples to run for testing 
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
    dir.create(paste0(dir.eur,"model_",i))
    new.dir.eur <- paste0(dir.eur,"model_",i)
    dir.create(paste0(pf.dir.eur,"model_",i))
    pf.new.dir.eur <- paste0(pf.dir.eur,"model_",i)
    dir.create(paste0(pf.base.eur,"model_",i))
    pf.base.new.eur <- paste0(pf.base.eur,"model_",i)
    
    sample.num <- readr::read_tsv(paste0(df.path.eur, "samples",i,".tsv"), col_names = TRUE)
    #for each batch exclude batch sampels from test set 
    eur.train <- eur.expr[rownames(eur.expr) %in% sample.num$x,]
    table(eur.train$subtype)
    
    # make sure there were enough samples for a test set - otherwise skip test for ancestry
    eur.test <- eur.expr[!(rownames(eur.expr) %in% sample.num$x),]
    eur.train$subtype <- as.factor(eur.train$subtype)
    if(nrow(eur.test) == 0){
      test.samples <- 0
    }else{
      test.samples <- 1
      eur.test$subtype <- as.factor(eur.test$subtype)
    }
    out_f <- paste0("model_",i)
    
    ## get base genes for phyloFrame
    base.eur.genes <- pf_top_varying_genes(eur.train, variable.genes)
    eur.base_pf <- phyloFrame(base.eur.genes, eur.train, pf.base.new.eur, out_f, en.mixture) ## PF BASELINE SIG - this gives back a model - we dont actually need it 
    
    ## read back in the signature to use genes as start for network walk
    model.genes <- read.table(paste0(pf.base.new.eur, "/",out_f,"_all_sig.txt"))
    model.genes <- model.genes[-1,]
    colnames(model.genes) <- c("Variable", "Importance", "Sign")
    model.genes <- model.genes$Variable
    if(mod.version == 2){
      model.input.genes <- get.genes.V2(network, model.genes, node, edge, eur.train, exomeAF) # CHANGED EXPRESSION -> EUR.TRAIN#get network genes for phyloframe and varying genes for benchmark for this batch
      ancestry.genes <- model.input.genes$phyloFrame
      temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
      #temp.anc.genes <- ancestry.genes
      anc.genes.dat <- as.data.frame(temp.anc.genes)
    	variable.genes <- length(anc.genes.dat$temp.anc.genes)
    	benchmark.genes <- pf_top_varying_genes(eur.train,length(temp.anc.genes))## ADDED
    	benchmark.genes <- c(benchmark.genes, model.genes)

      write.table(anc.genes.dat, file = paste0(ancestry.dir, "/eur_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
      eur_model <- phyloFrame(NULL, eur.train, new.dir.eur, out_f, run.penalties) #BENCHMARK
	## ADDED
	    pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")#ADDED
      # find max value in expression matrix 
      rescaled.eur.train <- ancestry.rescale(eur.train, temp.anc.genes)
      #keep.genes <- c(base.eur.genes,model.genes, temp.anc.genes, "subtype")
      #rescaled.eur.train <- rescaled.eur.train[,colnames(rescaled.eur.train) %in% keep.genes]
      #pf.in.genes <- c(pf.in.genes, "subtype")
      pf.eur.model <- phyloFrame(NULL, rescaled.eur.train, pf.new.dir.eur, out_f, run.penalties) #PHYLOFRAME
    
    }else{
      model.input.genes <- get.genes.V1(network, model.genes, node, edge, eur.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
      temp.train <- eur.train[,colnames(eur.train) %in% model.input.genes]
      top.var.anc <- pf_top_varying_genes(temp.train, 150)
      top.var.anc <- c(top.var.anc, model.genes, "subtype")
      #rescaled.eur.train <- ancestry.rescale(eur.train, c(anc.model.genes, model.genes))
      #eur.base.anc.pf <- phyloFrame(model.input.genes, eur.train, pf.base.new.eur, paste0(out_f,"_ancestry"), en.mixture) #
      
      ## read back in the signature to use genes as start for network walk
      # anc.model.genes <- read.table(paste0(pf.base.new.eur, "/",paste0(out_f,"_ancestry"),"_all_sig.txt"))
      # anc.model.genes <- anc.model.genes[-1,]
      # colnames(anc.model.genes) <- c("Variable", "Importance", "Sign")
      # anc.model.genes <- anc.model.genes$Variable
      
      
      #rescaled.eur.train <- ancestry.rescale(eur.train, c(anc.model.genes, model.genes))
      #num.genes <- length(anc.model.genes)
      #model.input.genes <- c(model.input.genes, model.genes, "subtype")
      
      expr.variance <- apply(eur.train, 2, var)
      var.ordered <- expr.variance[order(expr.variance, decreasing = TRUE)]
      genes.keep <- var.ordered[1:(length(top.var.anc ))] 
      n.variable.genes<- c(names(genes.keep), model.genes, "subtype")
      #n.variable.genes <- c(model.genes, n.variable.genes, "subtype")
    
      # ancestry.genes <- model.input.genes$phyloFrame # these get passed into phyloframe with no penalty
      # temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
      # anc.genes.dat <- as.data.frame(temp.anc.genes)
      # write.table(anc.genes.dat, file = paste0(ancestry.dir, "/eur_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
      # benchmark.genes <- model.input.genes$benchmark
      eur_model <- phyloFrame(n.variable.genes, eur.train, new.dir.eur, out_f, 0) #BENCHMARK
      pf.eur.model <- phyloFrame(top.var.anc ,eur.train, pf.new.dir.eur, out_f, 0) #PHYLOFRAME # anc.model.genes == top.var.anc
    }
    
    ### run model on other ancestries ###
    afr.out <- paste0("model_",i,"afr")
    model.metrics(eur_model, afr.expr, new.dir.eur, afr.out, subtype1, subtype2)
    model.metrics(pf.eur.model, afr.expr, pf.new.dir.eur, afr.out, subtype1, subtype2)
    eas.out <- paste0("model_",i,"eas")
    model.metrics(eur_model, eas.expr, new.dir.eur, eas.out, subtype1, subtype2)
    model.metrics(pf.eur.model, eas.expr, pf.new.dir.eur, eas.out, subtype1, subtype2)
    if(test.samples == 1){ #make sure we actually have samples to test on 
      eur.out <- paste0("model_",i,"eur")
      model.metrics(eur_model, eur.test, new.dir.eur, eur.out, subtype1, subtype2)
      model.metrics(pf.eur.model, eur.test, pf.new.dir.eur, eur.out, subtype1, subtype2)
    }
    admixed.out <- paste0("model_",i,"admixed")
    model.metrics(eur_model, admixed.expr, new.dir.eur, admixed.out, subtype1, subtype2)
    model.metrics(pf.eur.model, admixed.expr, pf.new.dir.eur, admixed.out, subtype1, subtype2)
    mixed.out <- paste0("model_",i,"mixed")
    temp.mixed.expr <- mixed.expr[!(rownames(mixed.expr) %in% sample.num),]
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
    base.afr.genes <- pf_top_varying_genes(afr.train, variable.genes)
    out_f <- paste0("model_",i)
    afr.base_pf <- phyloFrame(base.afr.genes, afr.train, pf.base.new.afr, out_f, en.mixture) ## PF BASELINE SIG - this gives back a model - we dont actually need it 
    ## read back in the signature to use genes as start for network walk
    model.genes <- read.table(paste0(pf.base.new.afr, "/",out_f,"_all_sig.txt"))
    model.genes <- model.genes[-1,]
    colnames(model.genes) <- c("Variable", "Importance", "Sign")
    model.genes <- model.genes$Variable
    ## get network and ordered genes for phyloFrame and most varying genes of same number for benchmark
    
    if(mod.version == 2){
      model.input.genes <- get.genes.V2(network, model.genes, node, edge, afr.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
      ancestry.genes <- model.input.genes$phyloFrame # these get passed into phyloframe with no penalty
      temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
      anc.genes.dat <- as.data.frame(temp.anc.genes)
      variable.genes <- length(anc.genes.dat$temp.anc.genes)
      benchmark.genes <- pf_top_varying_genes(afr.train,length(temp.anc.genes))## ADDED
      benchmark.genes <- c(benchmark.genes, model.genes)
      
      write.table(anc.genes.dat, file = paste0(ancestry.dir, "/afr_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
      afr_model <- phyloFrame(NULL, afr.train, new.dir.afr, out_f, run.penalties) #BENCHMARK
      rescaled.afr.train <- ancestry.rescale(afr.train, temp.anc.genes)
      pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")
      pf.afr.model <- phyloFrame(NULL, rescaled.afr.train, pf.new.dir.afr, out_f, run.penalties) #PHYLOFRAME - all 20000 get run 
    }else{
      model.input.genes <- get.genes.V1(network, model.genes, node, edge, afr.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
      ancestry.genes <- model.input.genes$phyloFrame # these get passed into phyloframe with no penalty
      temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
      anc.genes.dat <- as.data.frame(temp.anc.genes)
      write.table(anc.genes.dat, file = paste0(ancestry.dir, "/afr_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
      benchmark.genes <- model.input.genes$benchmark
      afr_model <- phyloFrame(benchmark.genes, afr.train, new.dir.afr, out_f, 0) #BENCHMARK
      length(ancestry.genes)
      pf.afr.model <- phyloFrame(ancestry.genes, afr.train, pf.new.dir.afr, out_f, 0) #PHYLOFRAME
      
    }
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
    temp.mixed.expr <- mixed.expr[!(rownames(mixed.expr) %in% sample.num),]
    model.metrics(afr_model, temp.mixed.expr, new.dir.afr, mixed.out, subtype1, subtype2)
    model.metrics(pf.afr.model, temp.mixed.expr, pf.new.dir.afr, mixed.out, subtype1, subtype2)
  }
  
  #################################### EAS MODEL RUN ON ALL ANCESTRIES ######################################
  print("################################### STARTING EAST ASIAN ANALYSIS ###################################")
  
  ## COMMENTED OUT BECAUSE UTERINE DOESNT HAVE ENOUGH EAST ASIAN TO TRAIN A SAMPLE 
  
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
      #eur.dat <- eur.clin[rownames(eur.clin) %in% sample.num$X1,]
      base.eas.genes <- pf_top_varying_genes(eas.train,variable.genes)
      #eas.model <- phyloFrame(n.variable.genes, eas.train, new.dir.eas, out_f, en.mixture) # BENCHMARK
      #pf.eas.model <- phyloFrame(top.anc.dat, eas.train, pf.new.dir.eas, out_f,en.mixture) #PHYLOFRAME
      eas.base_pf <- phyloFrame(base.eas.genes, eas.train, pf.base.new.eas, out_f, en.mixture) ## PF BASELINE SIG - this gives back a model - we dont actually need it 
      ## read back in the signature to use genes as start for network walk
      model.genes <- read.table(paste0(pf.base.new.eas, "/",out_f,"_all_sig.txt"))
      model.genes <- model.genes[-1,]
      colnames(model.genes) <- c("Variable", "Importance", "Sign")
      model.genes <- model.genes$Variable
      if(mod.version == 2){
        model.input.genes <- get.genes.V2(network, model.genes, node, edge, eas.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
        ancestry.genes <- model.input.genes$phyloFrame
        temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
        anc.genes.dat <- as.data.frame(temp.anc.genes)
        variable.genes <- length(anc.genes.dat$temp.anc.genes)
        benchmark.genes <- pf_top_varying_genes(eas.train,length(temp.anc.genes))## ADDED
        benchmark.genes <- c(benchmark.genes, model.genes)
        
        write.table(anc.genes.dat, file = paste0(ancestry.dir, "/eas_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
        eas.model <- phyloFrame(NULL, eas.train, new.dir.eas, out_f, run.penalties) #BENCHMARK
        rescaled.eas.train <- ancestry.rescale(eas.train, temp.anc.genes)
        pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")
        pf.eas.model <- phyloFrame(NULL, rescaled.eas.train, pf.new.dir.eas, out_f, run.penalties) #PHYLOFRAME
      }else{
        model.input.genes <- get.genes.V1(network, model.genes, node, edge, eas.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
        ancestry.genes <- model.input.genes$phyloFrame # these get passed into phyloframe with no penalty
        temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
        anc.genes.dat <- as.data.frame(temp.anc.genes)
        write.table(anc.genes.dat, file = paste0(ancestry.dir, "/eas_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
        benchmark.genes <- model.input.genes$benchmark
        eas.model <- phyloFrame(benchmark.genes, eas.train, new.dir.eas, out_f, 0) #BENCHMARK
        pf.eas.model <- phyloFrame(ancestry.genes, eas.train, pf.new.dir.eas, out_f, 0) #PHYLOFRAME
      }
      
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
      temp.mixed.expr <- mixed.expr[!(rownames(mixed.expr) %in% sample.num),]
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
    #eur.dat <- eur.clin[rownames(eur.clin) %in% sample.num$X1,]
    out_f <- paste0("model_",i)
    base.admixed.genes <- pf_top_varying_genes(admixed.train, variable.genes)
    #admixed.model <- phyloFrame(n.variable.genes, admixed.train, new.dir.adm, out_f,en.mixture) #BENCHMARK
    #pf.admixed.model <- phyloFrame(top.anc.dat, admixed.train, pf.new.dir.adm, out_f,en.mixture) #PHYLOFRAME
    admixed.base_pf <- phyloFrame(base.admixed.genes, admixed.train, pf.base.new.adm, out_f, en.mixture) ## PF BASELINE SIG - this gives back a model - we dont actually need it 
    ## read back in the signature to use genes as start for network walk
    model.genes <- read.table(paste0(pf.base.new.adm, "/",out_f,"_all_sig.txt"))
    model.genes <- model.genes[-1,]
    colnames(model.genes) <- c("Variable", "Importance", "Sign")
    model.genes <- model.genes$Variable
    if(mod.version == 2){
      model.input.genes <- get.genes.V2(network, model.genes, node, edge, admixed.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
      ancestry.genes <- model.input.genes$phyloFrame # these get passed into phyloframe with no penalty
      temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
      anc.genes.dat <- as.data.frame(temp.anc.genes)
      variable.genes <- length(anc.genes.dat$temp.anc.genes)
      benchmark.genes <- pf_top_varying_genes(admixed.train,length(temp.anc.genes))## ADDED
      benchmark.genes <- c(benchmark.genes, model.genes)
      
      write.table(anc.genes.dat, file = paste0(ancestry.dir, "/admixed_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
      admixed.model <- phyloFrame(NULL, admixed.train, new.dir.adm, out_f, run.penalties) #BENCHMARK
      rescaled.admixed.train <- ancestry.rescale(admixed.train, temp.anc.genes)
      pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")
      pf.admixed.model <- phyloFrame(NULL, rescaled.admixed.train, pf.new.dir.adm, out_f, run.penalties) #PHYLOFRAME
    }else{
      model.input.genes <- get.genes.V1(network, model.genes, node, edge, admixed.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
      ancestry.genes <- model.input.genes$phyloFrame # these get passed into phyloframe with no penalty
      temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
      anc.genes.dat <- as.data.frame(temp.anc.genes)
     
      
      write.table(anc.genes.dat, file = paste0(ancestry.dir, "/admixed_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
      benchmark.genes <- model.input.genes$benchmark
      admixed.model <- phyloFrame(benchmark.genes, admixed.train, new.dir.adm, out_f, 0 ) #BENCHMARK
      pf.admixed.model <- phyloFrame(ancestry.genes, admixed.train, pf.new.dir.adm, out_f, 0) #PHYLOFRAME
    }
    
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
    temp.mixed.expr <- mixed.expr[!(rownames(mixed.expr) %in% sample.num),]
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
    #eur.dat <- eur.clin[rownames(eur.clin) %in% sample.num$X1,]
    out_f <- paste0("model_",i)
    base.mixed.genes <- pf_top_varying_genes(mixed.train, variable.genes)
    mixed.base_pf <- phyloFrame(base.mixed.genes, mixed.train, pf.base.new.mix, out_f, en.mixture) ## PF BASELINE SIG - this gives back a model - we dont actually need it 
    ## read back in the signature to use genes as start for network walk
    model.genes <- read.table(paste0(pf.base.new.mix, "/",out_f,"_all_sig.txt"))
    model.genes <- model.genes[-1,]
    colnames(model.genes) <- c("Variable", "Importance", "Sign")
    model.genes <- model.genes$Variable
    ## get network and ordered genes for phyloFrame and most varying genes of same number for benchmark
    if(mod.version == 2){
      model.input.genes <- get.genes.V2(network, model.genes, node, edge, mixed.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
      ancestry.genes <- model.input.genes$phyloFrame
      temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)]
      anc.genes.dat <- as.data.frame(temp.anc.genes)
      variable.genes <- length(anc.genes.dat$temp.anc.genes)
      benchmark.genes <- pf_top_varying_genes(mixed.train,length(temp.anc.genes))## ADDED
      benchmark.genes <- c(benchmark.genes, model.genes)
      
      write.table(anc.genes.dat, file = paste0(ancestry.dir, "/mixed_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
      mixed.model <- phyloFrame(NULL, mixed.train, new.dir.mix, out_f,run.penalties) #BENCHMARK
      rescaled.mixed.train <- ancestry.rescale(mixed.train, temp.anc.genes)
      pf.pass.genes <- c(temp.anc.genes, model.genes, "subtype")#ADDED
      pf.mixed.model <- phyloFrame(NULL, rescaled.mixed.train, pf.new.dir.mix, out_f, run.penalties) #PHYLOFRAME
    }else{
      model.input.genes <- get.genes.V1(network, model.genes, node, edge, mixed.train, exomeAF) #get network genes for phyloframe and varying genes for benchmark for this batch
      ancestry.genes <- model.input.genes$phyloFrame # 
      temp.anc.genes <- ancestry.genes[!(ancestry.genes %in% model.genes)] #just for writing to file
      anc.genes.dat <- as.data.frame(temp.anc.genes)
      write.table(anc.genes.dat, file = paste0(ancestry.dir, "/mixed_",out_f,"_genes.txt" ), sep = "\t", col.names = TRUE, row.names = FALSE)
      benchmark.genes <- model.input.genes$benchmark
      mixed.model <- phyloFrame(benchmark.genes, mixed.train, new.dir.mix, out_f,0) #BENCHMARK
      pf.mixed.model <- phyloFrame(ancestry.genes, mixed.train, pf.new.dir.mix, out_f, 0) #PHYLOFRAME
    }
    
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
    temp.eur.expr <- eur.expr[!(rownames(eur.expr) %in% sample.num),]
    model.metrics(mixed.model, temp.eur.expr, new.dir.mix, eur.out, subtype1, subtype2)
    model.metrics(pf.mixed.model, temp.eur.expr, pf.new.dir.mix, eur.out, subtype1, subtype2)
    admixed.out <- paste0("model_",i,"admixed")
    temp.admixed.expr <- admixed.expr[!(rownames(admixed.expr) %in% sample.num),]
    model.metrics(mixed.model, temp.admixed.expr, new.dir.mix, admixed.out, subtype1, subtype2)
    model.metrics(pf.mixed.model, temp.admixed.expr, pf.new.dir.mix, admixed.out, subtype1, subtype2)
    if(test.samples == 1){
      mixed.out <- paste0("model_",i,"mixed")
      model.metrics(mixed.model, mixed.test, new.dir.mix, mixed.out, subtype1, subtype2)
      model.metrics(pf.mixed.model, mixed.test, pf.new.dir.mix, mixed.out, subtype1, subtype2)
    }
    
  }
}

##### CREATE PLOTS #### 
run(disease, model.path, "admixed", "admixed", en.mixture)
run(disease, model.path, "mixed_ancestry", "mixed", en.mixture)
run(disease, model.path, "eur", "eur", en.mixture)
run(disease, model.path, "eas", "eas", en.mixture)
run(disease, model.path, "afr", "afr", en.mixture)
#### SCALED IMPORTANCE HEATMAP ### - creates plots and saves them. 


