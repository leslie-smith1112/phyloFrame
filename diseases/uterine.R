uterine.init<- function(){
  expression <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/ucec_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt", col_names = TRUE)
  clinical <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/ucec_tcga_pan_can_atlas_2018_clinical_data.tsv")
  estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
  samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == "UCEC",]
  samples.ancestry <- samples.ancestry %>% dplyr::select(patient, consensus_ancestry)
  samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
  #narrow samples down to those actually in the expression matrix
  samples.ancestry  <- samples.ancestry[samples.ancestry$patient %in% colnames(expression),]
  
  network <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/disease_networks/uterine_endometrium_symbol.tsv")
  network <- network %>% dplyr::select(-4)
  de.meta <- clinical %>% dplyr::select(`Sample ID`, Subtype)
  colnames(de.meta) <- c("sample_id", "subtype")
  # dz.signature <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/uterine_cancer_signature.txt", col_names=FALSE)
  # dz.signature <- dz.signature$X1
  print("Uterine completed reading in all needed files.")
  uterine.values <- list("expr" = expression, "clin" = clinical, "samples.anc" = samples.ancestry, "net" = network, "diff.meta" = de.meta)
  return(uterine.values)
}

add.subtype.uterine <- function(expression, tcga.clinical){
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
  clin.cut <- data.frame( tcga.clinical$`Sample ID`,tcga.clinical$`Tumor Type`)
  colnames(clin.cut) <- c( "sample_id","subtype")
  #narrow down to subtypes we want 
  clin.cut <- clin.cut[clin.cut$subtype == "Endometrioid Endometrial Adenocarcinoma" | clin.cut$subtype == "Serous Endometrial Adenocarcinoma",]
  pattern1 <- "Endometrioid Endometrial Adenocarcinoma"
  pattern2 <- "Serous Endometrial Adenocarcinoma"
  
  clin.cut$subtype <- stri_replace_all_fixed(clin.cut$subtype, pattern1,"Endometrioid")
  clin.cut$subtype <- stri_replace_all_fixed(clin.cut$subtype, pattern2,"Serous")
  
  temp <- merge(x = tran, y = clin.cut, by = "sample_id")
  temp[is.na(temp)] <- 0
  #binomial2 <- na.omit(temp) #get rid of nas if there are any
  binomial2 <- temp
  binomial2$subtype <- as.factor(binomial2$subtype)
  
  cols <- c("sample_id", genes, "subtype")
  colnames(binomial2) <- cols
  samples <- binomial2$sample_id
  rownames(binomial2) <- samples 
  binomial2 <- binomial2[,-1]
  return(binomial2)
}

