


breast.init <- function(){
  expression <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/brca_data_mrna_seq_v2_rsem.txt", col_names = TRUE)
  clinical <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/brca_tcga_pan_can_atlas_2018_clinical_data.tsv", col_names = TRUE)
  estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
  samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == "BRCA",]
  samples.ancestry <- samples.ancestry %>% dplyr::select(patient, consensus_ancestry)
  samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
  samples.ancestry  <- samples.ancestry[samples.ancestry$patient %in% colnames(expression),]
  network <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/disease_networks/mammary_epithelium_symbol.tsv")
  network <- network %>% dplyr::select(-4)
  # dz.signature <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/breast_cancer_signature.txt", col_names=FALSE)
  # dz.signature <- dz.signature$X1
  print("Breast completed all read files.")
  breast.values <- list("expr" = expression, "clin" = clinical, "samples.anc" = samples.ancestry, "net" = network)
  return(breast.values)
}

add.subtype.breast <- function(expression, tcga.clinical){
  genes <- expression$Hugo_Symbol
  samples <- colnames(expression)
  samples <- samples[-1]
  head(samples)
  colnames(expression) <- NA
  expression <- expression[,-1] #get rid of genes - later added as row names 
  tran <- data.table::transpose(expression)
  tran <- log2(tran + 1)
  tran$sample_id <- samples
  ####################### -- unique to each cancer -- ####################### 
  #cut clinical data 
  clin.cut <- data.frame( tcga.clinical$`Sample ID`,tcga.clinical$Subtype)
  colnames(clin.cut) <- c( "sample_id","subtype")
  temp <- merge(x = tran, y = clin.cut, by = "sample_id")
  #keep.subtypes <- c("BRCA_Basal","BRCA_LumA","BRCA_LumB")
  #binomial2 <- temp[temp$subtype %in% keep.subtypes,]
  #binomial2 <- na.omit(binomial2)
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
  
  cols <- c("sample_id", genes, "subtype")
  colnames(binomial2) <- cols
  samples <- binomial2$sample_id
  rownames(binomial2) <- samples 
  binomial2 <- binomial2[,-1]
  return(binomial2)
}

