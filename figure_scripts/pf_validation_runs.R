## FIGURE 5 B-C ## 

#### RUN MODELS ON VALIDATION SETS ##### 

### pseudp line up ### 
# seperate data by african subtype
#run models on each one to see if we get accurate preductiion ofr botoh phyloFrame nad benchmark. 
#Read in all models 
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/expression_elasticnet.R")

full.expr <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/AACR_2023_poster/GSE211167_tpm_matrix.txt")
full.meta <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/AACR_2023_poster/TNBC_African_validation_set.tsv")
full.expr[1:5,1:5]

## get the ids we need from ID columns to match with expression matrix 
full.meta$Label <-  substr(full.meta$ID, 17,24)
#race: self-reported race: African American self-reported race: African/Ethiopian  self-reported race: African/Ghanaian
table(full.meta$race) # cut race to get only specific African Race
full.meta$race[1:7] <- substr(full.meta$race[1:7], 21,36)
full.meta$race[8:14] <- substr(full.meta$race[8:14], 21,37)
full.meta$race[15:22] <- substr(full.meta$race[15:22], 21,36)
full.meta$race[23:26] <- substr(full.meta$race[23:26], 21,37)


full.expr[1:5,1:5]
expr.dat <- as.data.frame(full.expr)
samp.names <- as.character(full.expr$Label)

rownames(expr.dat) <- samp.names
expr.dat <- expr.dat[,-1]
expr.dat[1:5,1:5]
expr.dat <- log2(expr.dat + 1) 
expr.dat$subtype <- "Basal" # all subtypes are Basal 

# seperate out African expression matrices for testing 
G.meta <- full.meta$Label[full.meta$race == "African/Ghanaian"]
Ghanaian.expr <- expr.dat[rownames(expr.dat) %in% G.meta,]
Ghanaian.expr$subtype <- as.factor(Ghanaian.expr$subtype)
A.meta <- full.meta$Label[full.meta$race == "African American"]
American.expr <- expr.dat[rownames(expr.dat) %in% A.meta,]
E.meta <- full.meta$Label[full.meta$race == "African/Ethiopian"]
Ethiopian.expr <- expr.dat[rownames(expr.dat) %in% E.meta,]

#Read in all models and test 

#define dirs 
dir <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinalVALIDATION/"
out.dir <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinalVALIDATION/validation_set/"

alias <- list("C17orf28" = "HID1", 
              "CYP2B7P1" = "CYP2B7P", 
              "C19orf21" = "MISP", 
              "FAM5C" = "BRINP3", 
              "CXorf61" = "CT83", 
              "C20orf114" = "BPIFB1",
              "C4orf7" = "FDCSP", 
              "ANKRD43" = "SOWAHA", 
              "FAM5C" = "BRINP3",
              "MS4A8B" = "MS4A8",
              "NEURL" = "NEURL1")

expr.dat <- column_to_rownames(full.expr, "Label")
expr.dat[1:5,1:5]
expr.dat$subtype <- "Basal"
expr.dat$subtype <- as.factor(expr.dat$subtype)

### #replace genes that were alias for the trained data #### 
names(expr.dat)[names(expr.dat) == "HID1"] <- "C17orf28"
names(expr.dat)[names(expr.dat) == "CYP2B7P"] <- "CYP2B7P1"
names(expr.dat)[names(expr.dat) == "MISP"] <- "C19orf21"
names(expr.dat)[names(expr.dat) == "BRINP3"] <- "FAM5C"
names(expr.dat)[names(expr.dat) == "CT83"] <- "CXorf61"
names(expr.dat)[names(expr.dat) == "BPIFB1"] <- "C20orf114"
names(expr.dat)[names(expr.dat) == "FDCSP"] <- "C4orf7"
names(expr.dat)[names(expr.dat) == "SOWAHA"] <- "ANKRD43"
names(expr.dat)[names(expr.dat) == "BRINP3"] <- "FAM5C"
names(expr.dat)[names(expr.dat) == "MS4A8"] <- "MS4A8B"
names(expr.dat)[names(expr.dat) == "NEURL1"] <- "NEURL"
names(expr.dat)[names(expr.dat) == "NSG2"] <- "HMP19"
names(expr.dat)[names(expr.dat) == "GPC1"] <- "GPC1-AS1"


# function to test all ancestry model on validation set, recall is calcualted and written to file
#there is a seperate function of model.metrics() defined in expression_elasticnet for the validation set
get.models <- function(dir,ancestry){
  pf.mod <- length(list.files(paste0(dir,"model_runs/phyloFrame/",ancestry,"/")))
  mod.num <- 1:pf.mod
  pf.mod.paths <- paste0(dir,"model_runs/phyloFrame/",ancestry,"/model_",mod.num,"/model_",mod.num,"_EN_model.rds")
  pf.models <- lapply(pf.mod.paths, read_rds)
  
  mod.paths <- paste0(dir,"model_runs/",ancestry,"/model_",mod.num,"/model_",mod.num,"_EN_model.rds")
  models <- lapply(mod.paths, read_rds)
  
  E.expr <- expr.dat[rownames(expr.dat) %in% E.meta,]
  A.expr <- expr.dat[rownames(expr.dat) %in% A.meta,]
  G.expr <- expr.dat[rownames(expr.dat) %in% G.meta,]
  
  for(i in 1:pf.mod){
    curr.mod <- pf.models[i][[1]]
    model.metrics(curr.mod, E.expr, out.dir, paste0("PF",ancestry,"_model_",i,"Ethiopian"), "Luminal","Basal")
    model.metrics(curr.mod, A.expr, out.dir, paste0("PF",ancestry,"_model_",i, "AfricanAmerican"), "Luminal","Basal")
    model.metrics(curr.mod, G.expr, out.dir, paste0("PF",ancestry,"_model_",i, "Ghanian"), "Luminal","Basal")
  }
  for(i in 1:pf.mod){
    curr.mod <- models[i][[1]]
    model.metrics(curr.mod, E.expr, out.dir, paste0("BM",ancestry,"_model_",i,"Ethiopian"), "Luminal","Basal")
    model.metrics(curr.mod, A.expr, out.dir, paste0("BM",ancestry,"_model_",i,"AfricanAmerican"), "Luminal","Basal")
    model.metrics(curr.mod, G.expr, out.dir, paste0("BM",ancestry,"_model_",i,"Ghanian"), "Luminal","Basal")
    
  }
  return(models)
}

get.models(dir, "eur")
get.models(dir,"afr")
get.models(dir,"eas")
get.models(dir,"admixed")
get.models(dir,"mixed")

#E.meta A.meta G.meta

######################## VALIDATION PLOTS  ######################## 

## for making the plots from validation set ## 

## file outlines: ## 
# BMadmixed_model_1AfricanAmerican_results.tsv
# BMadmixed_model_1Ethiopian_results.tsv
# BMadmixed_model_1Ghanian_results.tsv
# 
# PFadmixed_model_1AfricanAmerican_results.tsv
# PFadmixed_model_1Ethiopian_results.tsv
# PFadmixed_model_1Ghanian_results.tsv
#define model numbers
admixed <- 1
afr <- 1:2
eas <- 1
eur <- 1:17
mixed <- 1:6
ancestries <- c("admixed","afr","eas","eur","mixed")
dir <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinalVALIDATION/validation_set/"
afr.ancestry <- "AfricanAmerican"

bm.admixed <- paste0(dir,"BMadmixed_model_",admixed,afr.ancestry,"_results.tsv")
bm.admixed.files <- readr::read_tsv(bm.admixed)
#bm.admixed.files <- bm.admixed.files[1][[1]]
bm.admixed.files$ancestry <- "admixed"
colnames(bm.admixed.files) <- c("mode", "value", "ancestry")

bm.afr <- paste0(dir,"BMafr_model_",afr,afr.ancestry,"_results.tsv")
bm.afr.files <- lapply(bm.afr, readr::read_tsv)
bm.afr.files <- do.call(rbind,bm.afr.files)
bm.afr.files$ancestry <- "afr"
colnames(bm.afr.files) <- c("mode", "value","ancestry")

bm.eas <- paste0(dir,"BMeas_model_",eas,afr.ancestry,"_results.tsv")
bm.eas.files <- readr::read_tsv(bm.eas)
#bm.eas.files <- bm.eas.files[1][[1]]
bm.eas.files$ancestry <- "eas"
colnames(bm.eas.files) <- c("mode", "value", "ancestry")

bm.eur <- paste0(dir,"BMeur_model_",eur,afr.ancestry,"_results.tsv")
bm.eur.files <- readr::read_tsv(bm.eur)
#bm.eur.files <- bm.eur.files[1][[1]]
bm.eur.files$ancestry <- "eur"
colnames(bm.eur.files) <- c("mode", "value", "ancestry")

bm.mixed <- paste0(dir,"BMmixed_model_",mixed,afr.ancestry,"_results.tsv")
bm.mixed.files <- readr::read_tsv(bm.mixed)
#bm.mixed.files <- bm.mixed.files[1][[1]]
bm.mixed.files$ancestry <- "mixed"
colnames(bm.mixed.files) <- c("mode", "value", "ancestry")

BM.all <- rbind(bm.admixed.files, bm.afr.files, bm.eas.files, bm.eur.files, bm.mixed.files)
BM.all <- BM.all[BM.all$mode == "Recall",]
head(BM.all)
colnames(BM.all) <- c("mode","Benchmark","Model")

BM.all$model <- c("admixed1", "afr1", "afr2", "eas1",paste0("eur",1:17), paste0("mixed", 1:6))

pf.admixed <- paste0(dir,"PFadmixed_model_",admixed,afr.ancestry,"_results.tsv")
pf.admixed.files <- lapply(pf.admixed, readr::read_tsv)
pf.admixed.files <- do.call(rbind,pf.admixed.files)
pf.admixed.files$ancestry <- "admixed"
colnames(pf.admixed.files) <- c("mode", "value","ancestry")

pf.afr <- paste0(dir,"PFafr_model_",afr,afr.ancestry,"_results.tsv")
pf.afr.files <- lapply(pf.afr, readr::read_tsv)
pf.afr.files <- do.call(rbind,pf.afr.files)
pf.afr.files$ancestry <- "afr"
colnames(pf.afr.files) <- c("mode", "value","ancestry")

pf.eas <- paste0(dir,"PFeas_model_",eas,afr.ancestry,"_results.tsv")
pf.eas.files <- lapply(pf.eas, readr::read_tsv)
pf.eas.files <- do.call(rbind,pf.eas.files)
pf.eas.files$ancestry <- "eas"
colnames(pf.eas.files) <- c("mode", "value","ancestry")

pf.eur <- paste0(dir,"PFeur_model_",eur,afr.ancestry,"_results.tsv")
pf.eur.files <- lapply(pf.eur, readr::read_tsv)
pf.eur.files <- do.call(rbind,pf.eur.files)
pf.eur.files$ancestry <- "eur"
colnames(pf.eur.files) <- c("mode", "value","ancestry")

pf.mixed <- paste0(dir,"PFmixed_model_",mixed,afr.ancestry,"_results.tsv")
pf.mixed.files <- lapply(pf.mixed, readr::read_tsv)
pf.mixed.files <- do.call(rbind,pf.mixed.files)
pf.mixed.files$ancestry <- "mixed"
colnames(pf.mixed.files) <- c("mode", "value","ancestry")

PF.all <- rbind(pf.admixed.files, pf.afr.files, pf.eas.files, pf.eur.files, pf.mixed.files)
PF.all <- PF.all[PF.all$mode == "Recall",]
colnames(PF.all) <- c("mode", "PhyloFrame","ancestry")
PF.all$model <- c("admixed1", "afr1", "afr2", "eas1",paste0("eur",1:17), paste0("mixed", 1:6))
head(PF.all)

allsies <- cbind(PF.all, BM.all)
head(allsies)
allsies$testset <- afr.ancestry
allsies <- allsies[,-4]
to.add <- allsies %>% dplyr::select(ancestry, PhyloFrame,Benchmark, testset)
head(to.add)

mean(PF.all$PhyloFrame) #Ethiopian avg = 0.7777778
mean(BM.all$benchmark) #Ethiopian avg = 0.6969697
median(PF.all$PhyloFrame)
median(BM.all$benchmark)
head(PF.all)

temp <- reshape2::melt(to.add[,c("Benchmark", "PhyloFrame","ancestry")])

temp$variable <- factor(temp$variable, levels = c("PhyloFrame","Benchmark"))

ggplot(temp, aes(x=variable, y = value, fill = ancestry,color = ancestry)) + geom_dotplot(binaxis='y', stackdir='center',  dotsize=1.5,method="dotdensity", stackgroups = T, binpositions="all")+
  scale_fill_manual(values = c("afr" = "#1F619E",
                                "admixed"="#FDC652",
                                "eas"="#496849",
                                "eur" = "#CA4136",
                                "mixed" = "#8D7260")) + theme_bw()+scale_color_manual(values = c("afr" = "#1F619E",
                                                                                                 "admixed"="#FDC652",
                                                                                                 "eas"="#496849",
                                                                                                 "eur" = "#CA4136",
                                                                                                 "mixed" = "#8D7260")) +
  theme(axis.text=element_text(size=23),axis.title=element_text(size=3,face="bold")) + ylim(0,1)
