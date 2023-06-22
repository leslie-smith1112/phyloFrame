############ SIGNATURE SIGNATURE OVERLAP PLOT  ##############

library(ComplexHeatmap)
set.seed(123)
disease <- "breast"

#dir <- "/pf08_03thyroid_10000_V2_expr700_Netowrk_and_allele_change_ANCVAR_poster_version"
dir <- "/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal"
############ EUR ############
## define lists for number of models in each ancestry for help reading in files ## 
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

# disease <- "thyroid/pf4_thyroid_new_node2_5000V2_rescale500"

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
pf.admixed.names <- paste0(pf.admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt")
pf.admixed.myfiles <- lapply(pf.admixed.names, readr::read_tsv)
pf.mixed.names <- paste0(pf.mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt")
pf.mixed.myfiles <- lapply(pf.mixed.names, readr::read_tsv)

df.eur.names <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
eur.myfiles <- lapply(df.eur.names, readr::read_tsv)
df.afr.names <- paste0(afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
afr.myfiles <- lapply(df.afr.names, readr::read_tsv)
df.eas.names <- paste0(eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt") 
eas.myfiles <- lapply(df.eas.names, readr::read_tsv)
df.admixed.names <- paste0(admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt")
admixed.myfiles <- lapply(df.admixed.names, readr::read_tsv)
df.mixed.names <- paste0(mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt")
mixed.myfiles <- lapply(df.mixed.names, readr::read_tsv)


### SIGNAURE SIGNATURE CORRELATION ## 
# ------------------- EUR ------------------- # 
# data frame for phyloFrame and benchmark signtature percentages - 20 models total (17 eur, 2 afr, 1 eas)
eur.df <- data.frame(matrix(ncol = all.num,nrow = 0))
pf.eur.df <- data.frame(matrix(ncol = all.num,nrow = 0))
temp <- c()
pf.temp <- c()
for(i in 1:eur.count){

  for(j in 1:eur.count){
    ### check the intersection of model i with model j
    if(i != j){ #excluse overlap with itself 
      the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable))
      ### check how many genes total from the 2 models
      gene.union <- union(eur.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable)
      denom <- length(gene.union)
      overlap <- (the.int/denom) * 100 # get in percentage
      temp <- c(temp,overlap)
      ## phyloframe ##
      pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable))
      pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable)
      pf.denom <- length(pf.gene.union)
      pf.overlap <- (pf.the.int/pf.denom) * 100
      pf.temp <- c(pf.temp, pf.overlap)

    }
    print(length(pf.temp))
    print(length(temp))
  #
  }
  for(j in 1:afr.count){
    ### checking model1 of phyloframe against all other models of benchamrk
      the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable))
      gene.union <- union(eur.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable)
      denom <- length(gene.union)
      overlap <- (the.int/denom) * 100
      temp<- c(temp, overlap) # add 17 because there are 17 eur models already in dataframe
      ## phyloframe ##
      pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable))
      pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable)
      pf.denom <- length(pf.gene.union)
      pf.overlap <- (pf.the.int/pf.denom) * 100
      pf.temp <- c(pf.temp, pf.overlap) # add 17 because there are 17 eur models already in dataframe
  }

  print(length(pf.temp))
  print(length(temp))
  for(j in 1:eas.count){

    ### checking model1 of phyloframe against all other models of benchamrk
    the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable))
    gene.union <- union(eur.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp <- c(temp,overlap)
    ## phyloframe ##
    pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp <- c(pf.temp,pf.overlap)

  }
  for(j in 1:admixed.count){
    ### checking model1 of phyloframe against all other models of benchamrk
    the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,admixed.myfiles[j][[1]]$Variable))
    gene.union <- union(eur.myfiles[i][[1]]$Variable,admixed.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp <- c(temp,overlap)
    ## phyloframe ##
    pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.admixed.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.admixed.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp <- c(pf.temp,pf.overlap)
  }
  print(length(pf.temp))
  print(length(temp))
  

  print(length(pf.temp))
  print(length(temp))
  for(j in 1:mixed.count){
    ### checking model1 of phyloframe against all other models of benchamrk
    the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,mixed.myfiles[j][[1]]$Variable))
    gene.union <- union(eur.myfiles[i][[1]]$Variable,mixed.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp <- c(temp,overlap)
    ## phyloframe ## 
    pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.mixed.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.mixed.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp <- c(pf.temp,pf.overlap)

  }
  print(length(pf.temp))
  print(length(temp))
}
head(pf.temp)
mean(pf.temp)
mean(temp)

