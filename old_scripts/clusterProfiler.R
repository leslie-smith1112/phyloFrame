#############################################################################
#SCRIPT FOR GENE WEIGHT COMPARISON BETWEEN BENCHMARK AND PHYLOFRAME SIGNATURES. 
#############################################################################
library(pheatmap)
library(grid)
#################################################################################
library(clusterProfiler)
library("org.Hs.eg.db")
library(DESeq2)
library(msigdbr)
library(magrittr)
library(EnhancedVolcano)
library(apeglm)
library(ggplot2)

work.dir <- "pf08_03breast_10000_V2_expr700_Netowrk_and_allele_change_ANCVAR50_poster_version"
disease <- "breast"
# disease <- "thyroid" ## CHANGE HERE
# work.dir <- "pen1_version1_3" ## CHANGE HERE 
# get.matrix(disease, work.dir)
get.matrix <- function(disease, work.dir, plot.dir){
  dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/",disease,"/",work.dir,"/model_runs/")
  
  eur.num <- length(list.files(paste0(dir,"eur")))
  eur.num <- list(1:eur.num)
  eur.num <- unlist(eur.num)
  afr.num <- length(list.files(paste0(dir,"afr")))
  afr.num <- list(1:afr.num)
  afr.num <- unlist(afr.num)
  eas.num <- length(list.files(paste0(dir,"eas")))
  eas.num <- list(1:eas.num)
  eas.num <- unlist(eas.num)
  admixed.num <- length(list.files(paste0(dir,"admixed")))
  admixed.num <- list(1:admixed.num)
  admixed.num <- unlist(admixed.num)
  mixed.num <- length(list.files(paste0(dir,"mixed")))
  mixed.num <- list(1:mixed.num)
  mixed.num <- unlist(mixed.num)
  #############################################################################
  
  eur.count <- length(eur.num)
  afr.count <- length(afr.num)
  eas.count <- length(eas.num)
  admixed.count <- length(admixed.num)
  mixed.count <- length(mixed.num)
  
  eur.model.names <- paste0("eur_model_",eur.num)
  afr.model.names <- paste0("afr_model_",afr.num)
  eas.model.names <- paste0("eas_model_",eas.num)
  admixed.names <- paste0("admixed_",admixed.num)  
  mixed.names <- paste0("mixed_model_", mixed.num)
  
  ## read in signatures from phyloframe and bechmark ancestry models 
  pf.eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/eur/")
  pf.afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/afr/")
  pf.admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/admixed/")
  pf.eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/eas/")
  pf.mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/mixed/")
  
  eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/eur/")
  afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/afr/")
  admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/admixed/")
  eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/eas/")
  mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/mixed/")
  
  pf.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
  pf.eur.myfiles <- lapply(pf.eur.names, readr::read_tsv)
  pf.afr.names <- paste0(pf.afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
  pf.afr.myfiles <- lapply(pf.afr.names, readr::read_tsv)
  if(disease != "uterine"){
    pf.eas.names <- paste0(pf.eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt") 
    pf.eas.myfiles <- lapply(pf.eas.names, readr::read_tsv)
  }
  pf.admixed.names <- paste0(pf.admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt") 
  pf.admixed.myfiles <- lapply(pf.admixed.names, readr::read_tsv)
  pf.mixed.names <- paste0(pf.mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt") 
  pf.mixed.myfiles <- lapply(pf.mixed.names, readr::read_tsv)
  
  df.eur.names <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
  eur.myfiles <- lapply(df.eur.names, readr::read_tsv)
  df.afr.names <- paste0(afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
  afr.myfiles <- lapply(df.afr.names, readr::read_tsv)
  if(disease != "uterine"){
    df.eas.names <- paste0(eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt") 
    eas.myfiles <- lapply(df.eas.names, readr::read_tsv)
  }
  df.admixed.names <- paste0(admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt") 
  admixed.myfiles <- lapply(df.admixed.names, readr::read_tsv)
  df.mixed.names <- paste0(mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt") 
  mixed.myfiles <- lapply(df.mixed.names, readr::read_tsv)
  

  ## get list of all genes 
  
  
  add.model.names <- function(anc.list, ancestry, pf){
    dat <- data.frame(matrix(ncol = 5, nrow = 0))
    if(pf == TRUE){
      for (i in 1:length(anc.list)) {
        temp.dat <- anc.list[i]
        temp.dat <- as.data.frame(temp.dat)
        class(temp.dat)
        model.type <- rep(paste0("pf.",ancestry,".model.",i), nrow(temp.dat))
        temp.dat$model.type <- model.type
        temp.dat$scaled <- (temp.dat$Importance - min(temp.dat$Importance))/(max(temp.dat$Importance) - min(temp.dat$Importance))
        dat <- rbind(dat, temp.dat)
      }
      dat <- dat %>% dplyr::select(c(Variable, model.type, scaled))
    }else{
      for (i in 1:length(anc.list)) {
        temp.dat <- anc.list[i]
        temp.dat <- as.data.frame(temp.dat)
        class(temp.dat)
        model.type <- rep(paste0(ancestry,".model.",i), nrow(temp.dat))
        temp.dat$model.type <- model.type
        temp.dat$scaled <- (temp.dat$Importance - min(temp.dat$Importance))/(max(temp.dat$Importance) - min(temp.dat$Importance))
        dat <- rbind(dat, temp.dat)
        
      }
      dat <- dat %>% dplyr::select(c(Variable, model.type, scaled))
    }
    return(dat)
  }
  pf.eur.dat <- add.model.names(pf.eur.myfiles,"eur", TRUE)  ## repeat this for every ancestry
  pf.afr.dat <- add.model.names(pf.afr.myfiles,"afr", TRUE)
  if(disease != "uterine"){
    pf.eas.dat <- add.model.names(pf.eas.myfiles,"eas", TRUE)
  }
  pf.admixed.dat <- add.model.names(pf.admixed.myfiles,"admixed", TRUE)
  pf.mixed.dat <- add.model.names(pf.mixed.myfiles,"mixed", TRUE)
  
  eur.dat <- add.model.names(eur.myfiles,"eur", FALSE)  ## repeat this for every ancestry
  afr.dat <- add.model.names(afr.myfiles,"afr", FALSE)
  if(disease != "uterine"){
    eas.dat <- add.model.names(eas.myfiles,"eas", FALSE)
  }
  admixed.dat <- add.model.names(admixed.myfiles,"admixed", FALSE)
  mixed.dat <- add.model.names(mixed.myfiles,"mixed", FALSE)
  
  expression <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/brca_data_mrna_seq_v2_rsem.txt", col_names = TRUE)
  clinical <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/brca_tcga_pan_can_atlas_2018_clinical_data.tsv", col_names = TRUE)
  estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
  samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == "BRCA",]
  samples.ancestry <- samples.ancestry %>% dplyr::select(patient, consensus_ancestry)
  samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
  
  pf.eur.genes <- pf.eur.dat$Variable
  pf.afr.genes <- pf.afr.dat$Variable
  pf.eas.genes <- pf.eas.dat$Variable
  pf.admixed.genes <- pf.admixed.dat$Variable
  pf.mixed.genes <- pf.mixed.dat$Variable
  
  pf.all <- c(pf.eur.genes, pf.afr.genes,pf.eas.genes, pf.admixed.genes, pf.mixed.genes)
  
  eur.genes <- eur.dat$Variable
  afr.genes <- afr.dat$Variable
  eas.genes <- eas.dat$Variable
  admixed.genes <- admixed.dat$Variable
  mixed.genes <- mixed.dat$Variable
  all <- c(eur.genes, afr.genes, eas.genes, admixed.genes, mixed.genes)
  
  #####//function pass in expression matrix, ancestry ancestry genes and samples.ancestry and clinical #####
#  euro <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur",]
  #clusterProfiler <- function(expression, clinical, ancestry, samples.ancestry)
  clin.cut <- data.frame(clinical$`Sample ID`,clinical$Subtype)
  colnames(clin.cut) <- c( "sample_id","subtype")
  #eur.clin <- clin.cut[clin.cut$sample_id %in% colnames(eur.expr),]# make sure samples are in expression matrix
  
  ## keep only basal and luminal 
  temp <- clin.cut
  tcga.binomial <- temp[temp$subtype == "BRCA_Basal",] #160 patients
  tcga.binomial1 <- temp[temp$subtype == "BRCA_LumA",] #20 patients
  binomial <- rbind(tcga.binomial,tcga.binomial1)
  #binomial[1:5,1:5]
  tcga.binomial2 <- temp[temp$subtype == "BRCA_LumB",]
  #tcga.binomial2 <- na.omit(tcga.binomial2)
  binomial2 <- rbind(binomial,tcga.binomial2)
  binomial2 <- na.omit(binomial2)
  pattern1 <- "BRCA_LumA"
  pattern2 <- "BRCA_LumB"
  pattern3 <- "BRCA_Basal"
  subtype <- binomial2$subtype
  #### replace subtypes for logistic regression model####
  binomial2$subtype <- stri_replace_all_fixed(binomial2$subtype, pattern2,"Luminal")
  binomial2$subtype <- stri_replace_all_fixed(binomial2$subtype, pattern1,"Luminal")
  binomial2$subtype <- stri_replace_all_fixed(binomial2$subtype, pattern3,"Basal")
  binomial2$subtype <- as.factor(binomial2$subtype)
  meta.dat <- binomial2
  
  all.expr <- expression[expression$Hugo_Symbol %in% all, colnames(expression) %in% c(meta.dat$sample_id,"Hugo_Symbol")]
  en.expr <- column_to_rownames(pf.all.expr, "Hugo_Symbol")
  meta.dat <- meta.dat[match(colnames(en.expr), meta.dat$sample_id),]
  all.equal(colnames(en.expr), meta.dat$sample_id) ## CHANGE HERE
  # make sure subtype is a factor
  meta.dat$subtype <- as.factor(meta.dat$subtype) 
  rownames(meta.dat) <- NULL
  head(meta.dat)
  
  ## do differentual expression analysis 
  deseq_object <- compute_DE(en.expr, meta.dat)
  write.table(deseq_object,"/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/AACR_2023_poster/benchmark_allanc_de.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
  pf.deseq_object <- deseq_object
  
  pf.og.gene.list <- deseq_object$log2FoldChange 
  names(pf.og.gene.list) <- deseq_object$Gene
  gene.list <- na.omit(pf.og.gene.list)
  gene.list <- gene.list[order(gene.list, decreasing = TRUE)]
  organism <- "org.Hs.eg.db"
  keytypes(org.Hs.eg.db)
  gse <- gseGO(geneList=gene.list, 
               ont ="ALL", 
               keyType = "SYMBOL", 
               #nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
  write.table(gse,"/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/AACR_2023_poster/admixedall_BRCA_gse.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
  
  dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  ridgeplot(gse) + labs(x = "enrichment distribution")
}