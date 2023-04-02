#1. randomly select samples for subsample runs through elasticnet model
set.seed(4831)

#dir <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/samples/"
create.batches <- function(directory, cancer.type, expression, subtype1, subtype2){
  print(cancer.type)
  #get samples and count of each ancestry in given cancer 
  estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
  samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == cancer.type,]
  samples.ancestry <- samples.ancestry %>% dplyr::select(patient, consensus_ancestry)
  #add -01 to match the sample numbers to the tcga data
  samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
  #narrow samples down to those in expression matrix
  samples.ancestry  <- samples.ancestry[samples.ancestry$patient %in% rownames(expression),]
  #get count of each ancestry in the samples 
  eur.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur",]$patient
  eur.count <- length(eur.samples)
  afr.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr",]$patient
  afr.count <- length(afr.samples)
  eas.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eas",]$patient
  eas.count <- length(eas.samples)
  admix.samples.eas <- samples.ancestry[samples.ancestry$consensus_ancestry == "eas_admix",]$patient
  admix.samples.afr <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr_admix",]$patient
  admix.samples.eur <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur_admix",]$patient
  admix.samples.gen <- samples.ancestry[samples.ancestry$consensus_ancestry == "admix",]$patient
  admix.samples <- c(admix.samples.eas, admix.samples.afr, admix.samples.eur, admix.samples.gen)
  admix.count <- length(admix.samples)
  
  ### downsampling european for the mixed samples ### 
  eur.values <- list(1:length(eur.samples))
  eur.values <- unlist(eur.values)
  
  ##  we keep the sme number of euro samples as the max number of samples from the other ancestries 
  max.of.samples <- max(eas.count, afr.count, admix.count) 
  ## randomly sample european samples to keep
  keep.index <- sample(eur.values, max.of.samples, replace = FALSE)
  # - find samples to keep - # 
  eur.shortened <- eur.samples[keep.index]
  # we keep all samples & all.count to make mixed sample batches later
  all.samples <- c(eur.shortened, afr.samples, eas.samples, admix.samples)
  all.count <- length(eur.shortened) + afr.count + eas.count + admix.count
  
  # find smaller sample of the ancestries - this is how we will create our batches 
  smallest.batch <- min(eur.count, afr.count, eas.count, admix.count)
  
  #Expression matrices for each subtype that we can pull rows from 
  ##doing this to ensure that each batch has some of each subtype
  temp.expr1 <- expression[expression$subtype == subtype1,]
  temp.expr2 <- expression[expression$subtype == subtype2,]
  
  #for each ancestry make batches to match smallest size - all other ancestries have same format
  ##############EUR################# 
  #find the number of batches needed 
  eur.num.batches <- floor(eur.count/smallest.batch)
  
  #get ancstry samples in each subtype 
  sub1.eur <- rownames(temp.expr1)[rownames(temp.expr1) %in% eur.samples]
  sub2.eur <- rownames(temp.expr2)[rownames(temp.expr2) %in% eur.samples]
  
  #find number of samples in each subtype
  sub1.len <- length(sub1.eur)
  sub2.len <- length(sub2.eur)
  #get how many samples will be in each batch for subtype 1 
  subtype.one.count <- floor(sub1.len/eur.num.batches) 
  #how many samples will be leftover - these samples will be added individually to batches until there is no more
  subtype.one.left <- sub1.len - (subtype.one.count * eur.num.batches)
  
  #get how many samples will be in each batch for subtype 2
  subtype.two.count <- floor(sub2.len/eur.num.batches)
  #how many samples will be leftover - these samples will be added individually to batches until there is no more
  subtype.two.left <- sub2.len - (subtype.two.count * eur.num.batches)
  
  #for each subtype set up to randomly select samples from each subtype to be in a batch
  #lists to hold row number for each sample
  sub1.values <- list(1:sub1.len)
  sub1.values <- unlist(sub1.values)
  sub2.values <- list(1:sub2.len)
  sub2.values <- unlist(sub2.values)
  #for as many batches as we need, randomly select rownumber for samples to be in batch, 
  #if there is a left over sample, assign it to current batch (by taking one extra sample)
  #randomly sample rownumbers from the master list
  #remove chosen samples from master list so they arent selected again
  #repeat with subtype 2
  #
  for(i in 1:eur.num.batches){
    #subtype one 
    if(subtype.one.left != 0){
      size.batch <- subtype.one.count + 1
      subtype.one.left <- subtype.one.left - 1
    }else{
      size.batch <- subtype.one.count
    }
    sample1.row <- sample(x  = sub1.values, size = size.batch, replace = FALSE)
    sub1.values <- sub1.values[!(sub1.values %in% sample1.row)]
    current.sub1 <- sub1.eur[sample1.row]
    
    #subtype 2 
    if(subtype.two.left != 0){
      size.batch <- subtype.two.count + 1
      subtype.two.left <- subtype.two.left - 1
    }else{
      size.batch <- subtype.two.count
    }
    
    sample2.row <- sample(x  = sub2.values, size = size.batch, replace = FALSE)
    sub2.values <- sub2.values[!(sub2.values %in% sample2.row)]
    current.sub2 <- sub2.eur[sample2.row]
    #combine subtype 1 and subtype 2 samples and write to file
    all.subs <- c(current.sub1, current.sub2)
    print(paste0(directory,"eur/samples",i,".tsv"))
    write.table(all.subs, paste0(directory,"eur/samples",i,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  
  ##############AFR#################
  afr.num.batches <- floor(afr.count/smallest.batch)
  #afr.samples.left <-  afr.count %% afr.num.batches
  #samples in each subtype 
  sub1.afr <- rownames(temp.expr1)[rownames(temp.expr1) %in% afr.samples]
  sub2.afr <- rownames(temp.expr2)[rownames(temp.expr2) %in% afr.samples]
  
  sub1.len <- length(sub1.afr)
  sub2.len <- length(sub2.afr)
  subtype.one.count <- floor(sub1.len/afr.num.batches) #how many samples will be in each batch 
  subtype.one.left <- sub1.len - (subtype.one.count * afr.num.batches)  #how many samples will be left 
  
  subtype.two.count <- floor(sub2.len/afr.num.batches)
  subtype.two.left <- sub2.len - (subtype.two.count * afr.num.batches)
  
  sub1.values <- list(1:sub1.len)
  sub1.values <- unlist(sub1.values)
  sub2.values <- list(1:sub2.len)
  sub2.values <- unlist(sub2.values)
  #afr.batches <- data.frame(matrix(ncol = afr.num.batches, nrow = 0))
  for(i in 1:afr.num.batches){
    #subtype one 
    if(subtype.one.left != 0){
      size.batch <- subtype.one.count + 1
      subtype.one.left <- subtype.one.left - 1
    }else{
      size.batch <- subtype.one.count
    }
    sample1.row <- sample(x  = sub1.values, size = size.batch, replace = FALSE)
    sub1.values <- sub1.values[!(sub1.values %in% sample1.row)]
    current.sub1 <- sub1.afr[sample1.row]
    
    #subtype 2 
    if(subtype.two.left != 0){
      size.batch <- subtype.two.count + 1
      subtype.two.left <- subtype.two.left - 1
    }else{
      size.batch <- subtype.two.count
    }
    
    sample2.row <- sample(x  = sub2.values, size = size.batch, replace = FALSE)
    sub2.values <- sub2.values[!(sub2.values %in% sample2.row)]
    current.sub2 <- sub2.afr[sample2.row]
    
    all.subs <- c(current.sub1, current.sub2)
    write.table(all.subs, paste0(directory,"afr/samples",i,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  ##############EAS#################
  eas.num.batches <- floor(eas.count/smallest.batch)
  #eas.samples.left <-  eas.count %% eas.num.batches
  #samples in each subtype 
  sub1.eas <- rownames(temp.expr1)[rownames(temp.expr1) %in% eas.samples]
  sub2.eas <- rownames(temp.expr2)[rownames(temp.expr2) %in% eas.samples]
  
  sub1.len <- length(sub1.eas)
  sub2.len <- length(sub2.eas) ##TODO 
  subtype.one.count <- floor(sub1.len/eas.num.batches) #how many samples will be in each batch  
  subtype.one.left <- sub1.len - (subtype.one.count * eas.num.batches)  #how many samples will be left 
  
  subtype.two.count <- floor(sub2.len/eas.num.batches)
  subtype.two.left <- sub2.len - (subtype.two.count * eas.num.batches) 
  
  sub1.values <- list(1:sub1.len)
  sub1.values <- unlist(sub1.values)
  sub2.values <- list(1:sub2.len)
  sub2.values <- unlist(sub2.values)
  #eas.batches <- data.frame(matrix(ncol = eas.num.batches, nrow = 0))
  for(i in 1:eas.num.batches){
    #subtype one 
    if(subtype.one.left != 0){
      size.batch <- subtype.one.count + 1 #TODO 
      subtype.one.left <- subtype.one.left - 1
    }else{
      size.batch <- subtype.one.count
    }
    sample1.row <- sample(x  = sub1.values, size = size.batch, replace = FALSE)
    sub1.values <- sub1.values[!(sub1.values %in% sample1.row)]
    current.sub1 <- sub1.eas[sample1.row]
    
    #subtype 2 
    if(subtype.two.left != 0){
      size.batch <- subtype.two.count + 1
      subtype.two.left <- subtype.two.left - 1
    }else{
      size.batch <- subtype.two.count
    }
    
    sample2.row <- sample(x  = sub2.values, size = size.batch, replace = FALSE)
    sub2.values <- sub2.values[!(sub2.values %in% sample2.row)]
    current.sub2 <- sub2.eas[sample2.row]
    
    all.subs <- c(current.sub1, current.sub2)
    write.table(all.subs, paste0(directory,"eas/samples",i,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  
  ##############ADMIXED#################
  #EDITED 1/6/23 ADD BACK IN ADMIX SAMPLES IN MIN SAMPLE CALCULATION
  ## EDITED 1/1/23 - admixed is typically the smallest sample size, we want a larger sample size so
  #we just make one batch of admixed samples 
  admix.num.batches <- floor(admix.count/smallest.batch) #EDITED
  #i <- 1 - COMMENTED OUT 
  #admix.num.batches <- 1 - COMMENTED OUT 
  #admix.samples.left <-  admix.count %% admix.num.batches #EDITED
  #samples in each subtype 
  sub1.admix <- rownames(temp.expr1)[rownames(temp.expr1) %in% admix.samples]
  sub2.admix <- rownames(temp.expr2)[rownames(temp.expr2) %in% admix.samples]
  ## #EDITED 1/6/23 ADD BACK IN ADMIX SAMPLES IN MIN SAMPLE CALCULATION
  #all.subs <- c(sub1.admix, sub2.admix) - COMMENTED OUT 
  
  sub1.len <- length(sub1.admix)
  sub2.len <- length(sub2.admix)
  subtype.one.count <- floor(sub1.len/admix.num.batches) #how many samples will be in each batch for subtype one 
  subtype.one.left <- sub1.len - (subtype.one.count * admix.num.batches)  #how many samples will be left

  subtype.two.count <- floor(sub2.len/admix.num.batches)
  subtype.two.left <- sub2.len - (subtype.two.count * admix.num.batches)

  sub1.values <- list(1:sub1.len)
  sub1.values <- unlist(sub1.values)
  sub2.values <- list(1:sub2.len)
  sub2.values <- unlist(sub2.values)
  #admix.batches <- data.frame(matrix(ncol = admix.num.batches, nrow = 0))
  for(i in 1:admix.num.batches){
    #subtype one
    if(subtype.one.left != 0){
      size.batch <- subtype.one.count + 1
      subtype.one.left <- subtype.one.left - 1 #TODO
    }else{
      size.batch <- subtype.one.count
    }
    sample1.row <- sample(x  = sub1.values, size = size.batch, replace = FALSE)
    sub1.values <- sub1.values[!(sub1.values %in% sample1.row)]
    current.sub1 <- sub1.admix[sample1.row]

    #subtype 2
    if(subtype.two.left != 0){
      size.batch <- subtype.two.count + 1
      subtype.two.left <- subtype.two.left - 1
    }else{
      size.batch <- subtype.two.count
    }

    sample2.row <- sample(x  = sub2.values, size = size.batch, replace = FALSE)
    sub2.values <- sub2.values[!(sub2.values %in% sample2.row)]
    current.sub2 <- sub2.admix[sample2.row]

    all.subs <- c(current.sub1, current.sub2)
    write.table(all.subs, paste0(directory,"admixed/samples",i,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  #write.table(all.subs, paste0(directory,"admixed/samples",i,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
  ##############MIX OF ALL ANCESTRIES#################
  all.num.batches <- floor(all.count/smallest.batch)
  #all.samples.left <-  all.count %% all.num.batches
  #samples in each subtype 
  sub1.all <- rownames(temp.expr1)[rownames(temp.expr1) %in% all.samples]
  sub2.all <- rownames(temp.expr2)[rownames(temp.expr2) %in% all.samples]
  
  sub1.len <- length(sub1.all)
  sub2.len <- length(sub2.all) #TODO 
  subtype.one.count <- floor(sub1.len/all.num.batches) #how many samples will be in each batch 
  subtype.one.left <- sub1.len - (subtype.one.count * all.num.batches)  #how many samples will be left 
  
  subtype.two.count <- floor(sub2.len/all.num.batches)
  subtype.two.left <- sub2.len - (subtype.two.count * all.num.batches)#TODO 
  
  sub1.values <- list(1:sub1.len)
  sub1.values <- unlist(sub1.values)
  sub2.values <- list(1:sub2.len)
  sub2.values <- unlist(sub2.values)
  #all.batches <- data.frame(matrix(ncol = all.num.batches, nrow = 0))
  for(i in 1:all.num.batches){
    #subtype one 
    if(subtype.one.left != 0){
      size.batch <- subtype.one.count + 1 #TODO 
      subtype.one.left <- subtype.one.left - 1
    }else{
      size.batch <- subtype.one.count
    }
    sample1.row <- sample(x  = sub1.values, size = size.batch, replace = FALSE)
    sub1.values <- sub1.values[!(sub1.values %in% sample1.row)]
    current.sub1 <- sub1.all[sample1.row]
    
    #subtype 2 
    if(subtype.two.left != 0){
      size.batch <- subtype.two.count + 1
      subtype.two.left <- subtype.two.left - 1
    }else{
      size.batch <- subtype.two.count
    }
    
    sample2.row <- sample(x  = sub2.values, size = size.batch, replace = FALSE)
    sub2.values <- sub2.values[!(sub2.values %in% sample2.row)]
    current.sub2 <- sub2.all[sample2.row]
    
    all.subs <- c(current.sub1, current.sub2)
    write.table(all.subs, paste0(directory,"mixed_ancestry/samples",i,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  # we pass back the sample id for all samples to ccount for the downsampling of european sampels 
  ancestry.counts <- list("afr" = afr.num.batches, "eas" = eas.num.batches, "eur" = eur.num.batches, "admixed" = admix.num.batches, "mixed" = all.num.batches, "mixed_samples" = all.samples)
  return(ancestry.counts)
}

# 
# #2. get vector of value from 1 - number of samples in an ancestry
# #these will be the indexes we pull from each expression matrix
# values <- list(1:665) #- for eur this was set to 665, afr was 90
# values <- unlist(values)
# 
# #4. from phyloFrame main read in sample.ancestry, used to seperate samples by ancestry ### 
# eur.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eur",]$patient
# afr.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "afr",]$patient
# eas.samples <- samples.ancestry[samples.ancestry$consensus_ancestry == "eas",]$patient
# 
# ## read in the expression matrices to be cut ## 
# #dat  <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/brca_data_mrna_seq_v2_rsem.txt")
# dat <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/thca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt", col_names = TRUE)
# the.mat <- trim.expr.matrix_OG(dat, eur.samples)
# the.mat <- add.subtype.thyroid(the.mat, clinical)
# samples <- rownames(the.mat)
# 
# ## select basal and luminal samples seperatey to make sure each batch has enough - ## 
# 
# basal <- the.mat[the.mat$subtype == "M0",]
# luminal <- the.mat[the.mat$subtype == "MX",]
# 
# ## for eur: total luminal samples: 558, total basal samples: 107 ##
# ### we choose 6 basal 12 times and 7 basal samples 5 times ## 
# ### we choose 32 samples 3 times and 33 samples 14 times ## 
# #### this is done manually, segments for each range of quantity we needed: 1-12, 13-14, 15-17 ## 
# bas.values <- list(1:37)
# bas.values <- unlist(bas.values)
# lum.values <- list(1:53)
# lum.values <- unlist(lum.values)
# for(i in 2:2){
#   
#   bas.sample <- sample(x  = bas.values, size = 19, replace = FALSE)
#   bas.values <- bas.values[!(bas.values %in% bas.sample)]
#   length(bas.values)
#   lum.sample <- sample(x  = lum.values, size = 27, replace = FALSE)
#   lum.values <- lum.values[!(lum.values %in% lum.sample)]
#   length(lum.values)
#   #keep rows we chose ## 
#   temp.basal <- basal[bas.sample,]
#   temp.luminal <- luminal[lum.sample,]
#   
#   all <- rbind(temp.basal, temp.luminal)
#   to.keep <- rownames(all)
#   write.table(to.keep,paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/samples/afr/samples",i,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
#   
# }
# 
# 
# ## old from when we were selecting row number, writing to file and then reading back in and getting sampled ##  
# # 
# 
# values <- list(1:352)
# values <- unlist(values)
# for(i in 11:19){
#   sample1 <- sample(x = values, size = 18,replace = FALSE)
#   values <- values[!(values %in% sample1)]
#   print(length(values))
#   write.table(sample1,paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/thyroid/samples/eur/samples",i,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# #   dat <- files[i][[1]] #get list of samples
#   kept <- samples[dat$x]
#   write.table(kept,paste0(df.path,"sampleID_",i,".tsv"), sep = "\t", row.names = FALSE, col.names = FALSE)
#}
# 
# #7. check samples to make sure they have some basal samples - this prints out the subtype distribution 
# for(i in 1:17){
#   sample.num <- readr::read_tsv(paste0(df.path, "sampleID_",i,".tsv"), col_names = FALSE)
#   dat <- eur.mat[rownames(eur.mat) %in% sample.num$X1,]
#   #dat <- trim.expr.matrix_OG(expression, sample.num$X1)
#   #dat <- add.subtype(dat, clinical)
#   print(paste0("sample: ", i, "subtype distribution: ")) 
#   print(table(dat$subtype))
#   dat[1:5,1:5]
#   print("#####################################################################")
# }


## old from when we tried up sampling ## 
# #1.5 - added, upsampling each ancestry. 
# eas.t <- rownames_to_column(eas.clin,"sample_id")
# eas.up <- upSample(eas.t, eas.clin$subtype)
# eas.up[1:5,1:5]
# write.table(eas.up,"/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/samples/upsampled_matrices/eas_expression.tsv",sep = "\t", col.names = TRUE, row.names = FALSE)
# 
# eur.t <- rownames_to_column(eur.clin,"sample_id")
# eur.up <- upSample(eur.t, eur.clin$subtype)
# eur.up[1:5,1:5]
# dim(eur.up)
# table(eur.up$subtype)
# write.table(eur.up,"/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/samples/upsampled_matrices/eur_expression.tsv",sep = "\t", col.names = TRUE, row.names = FALSE)
# 
# afr.t <- rownames_to_column(afr.clin,"sample_id")
# afr.up <- upSample(afr.t, afr.clin$subtype)
# afr.up[1:5,1:5]
# table(afr.up$subtype)
# write.table(afr.up,"/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/breast/samples/upsampled_matrices/afr_expression.tsv",sep = "\t", col.names = TRUE, row.names = FALSE)
# #### end of add froom upsampling ####













  
  
  
  
