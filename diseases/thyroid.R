

thyroid.init<- function(){
  expression <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/thca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt", col_names = TRUE)
  clinical <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/thca_tcga_pan_can_atlas_2018_clinical_data.tsv")
  estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
  samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == "THCA",]
  samples.ancestry <- samples.ancestry %>% dplyr::select(patient, consensus_ancestry)
  samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
  #narrow samples down to those actually in the expression matrix
  samples.ancestry  <- samples.ancestry[samples.ancestry$patient %in% colnames(expression),]
  network <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/disease_networks/thyroid_gland_symbol.tsv")
  network <- network %>% dplyr::select(-4)
  # dz.signature <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/thyroid_cancer_signature.txt", col_names=FALSE)
  # dz.signature <- dz.signature$X1
  print("Thyroid completed all read files.")
  thyroid.values <- list("expr" = expression, "clin" = clinical, "samples.anc" = samples.ancestry, "net" = network) 
  return(thyroid.values)
}

add.subtype.thyroid <- function(expression, tcga.clinical){
  genes <- expression$Hugo_Symbol
  samples <- colnames(expression)
  samples <- samples[-1] #get rid of Hugo Symbol at the beginning of the list
  head(samples)
  colnames(expression) <- NA
  expression <- expression[,-1] #get rid of genes - later added as col names 
  tran <- data.table::transpose(expression)
  tran <- log2(tran + 1)
  tran$sample_id <- samples
  ####################### -- unique to each cancer -- ####################### 
  #cut clinical data 
  clin.cut <- data.frame( tcga.clinical$`Sample ID`,tcga.clinical$`American Joint Committee on Cancer Metastasis Stage Code`)
  colnames(clin.cut) <- c( "sample_id","subtype")
  #narrow down to subtypes we want 
  clin.cut <- clin.cut[clin.cut$subtype == "M0" | clin.cut$subtype == "MX",]
  
  temp <- merge(x = tran, y = clin.cut, by = "sample_id")
  binomial2 <- na.omit(temp) #get rid of nas if there are any
  
  binomial2$subtype <- as.factor(binomial2$subtype)
  
  cols <- c("sample_id", genes, "subtype")
  colnames(binomial2) <- cols
  samples <- binomial2$sample_id
  rownames(binomial2) <- samples 
  binomial2 <- binomial2[,-1]
  return(binomial2)
}

