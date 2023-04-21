########################################################################
### EXOME ALLELE FREQUENCY ORDERING ###
#######################################################################

require(tibble)
require(readr)  # for read_csv()
require(dplyr) #data frame handling 
require(tidyr)
require(stringi)
#source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/phyloFrame_main.R")

#######################################################################
### FUNCITON TO GET GENES IN NETWORK AND SORT EACH ANCESTRY IN DECREASING ORDER BY AF###
#######################################################################
pf_top_varying_genes <- function(expression.dat,top.genes){
  expr.variance <- apply(expression.dat, 2, var)
  var.ordered <- expr.variance[order(expr.variance, decreasing = TRUE)]
  genes.keep <- var.ordered[1:top.genes]
  genes.keep <- names(genes.keep)
  genes.keep <- c(genes.keep,"subtype")
  return(genes.keep)
}


order_frequencies <- function(include_genes, exome_file, expression){
  #### keep only include genes that are relevant in network ####
  no_sex <- exome_file[exome_file$gene %in% include_genes,] # genes from the human base network 
  keep.num <- 50
  upper.bound <- 1 #1#
  lower.bound <- 0.001 #0.5#
  
  mutationburden.upper <- 10000
  mutationburden.lower <- 0
  #### dataframe to hold the genes selected for each ancestry ####
  #dat <- list()
  
  ########################################################################
  ### For each ancestry, get range where right peak of enhanced frequencies
  ### is highest and keep only genes with snps in that range. Add genes to
  ### list holding all ancestry selected genes. 
  #######################################################################

  ########################################################################
  ### AFR ###
  #######################################################################
  #--- get enhanced genes for each ancestry ---#
 
  
  afr_snps <- no_sex[no_sex$afr > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_onf < 0 & no_sex$amr < 0 & no_sex$eas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 &
                        no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$eas_oea < 0 & no_sex$nfe_est < 0 & no_sex$nfe < 0 & no_sex$fin < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  
  afr.t <- table(afr_snps$gene)
  afr.t <- afr.t[order(afr.t, decreasing = TRUE)]
  afr.k <- afr.t[afr.t < mutationburden.upper & afr.t > mutationburden.lower]
  afr.nsnps <- afr_snps[afr_snps$gene %in% names(afr.k),] # only keep snps with bounded mutation burden 
  afr.snps.k <- afr.nsnps[(afr.nsnps$afr > lower.bound & afr.nsnps$afr < upper.bound),] # only keeps snps within enhanced allele frequency range
  
  #afr.o <- afr_snps[order(afr_snps$afr, decreasing = TRUE),]
  #afr.keep <- unique(na.omit(afr.snps.k$gene))
  afr.keep <- na.omit(afr.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% afr.keep]
  afr.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  #afr.keep <- afr.keep[1:keep.num] # keep top most enhanced genes 
  
  dat <- afr.keep
  

  
  
  ########################################################################
  ### EAS ###
  #######################################################################
  eas_snps <- no_sex[no_sex$eas > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_onf < 0 & no_sex$amr < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 &
                       no_sex$nfe_est < 0 & no_sex$nfe < 0 & no_sex$fin < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  eas.t <- table(eas_snps$gene)
  eas.t <- eas.t[order(eas.t, decreasing = TRUE)]
  eas.k <- eas.t[eas.t < mutationburden.upper & eas.t > mutationburden.lower]
  eas.nsnps <- eas_snps[eas_snps$gene %in% names(eas.k),]
  eas.snps.k <- eas.nsnps[(eas.nsnps$eas > lower.bound & eas.nsnps$eas < upper.bound),]
  
  #eas.o <- eas_snps[order(eas_snps$eas, decreasing = TRUE),]
  eas.keep <- na.omit(eas.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% eas.keep]
  eas.keep <- pf_top_varying_genes(temp.expr, keep.num)
  #eas.keep <- eas.keep[1:keep.num] # keep top most enhanced genes 

  dat <- c(dat,eas.keep)
  
  # ########################################################################
  # ### NFE ###
  # #######################################################################
  nfe_snps <- no_sex[no_sex$nfe > 0 & no_sex$sas < 0 & no_sex$amr < 0 & no_sex$afr < 0 & no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$eas_oea < 0 & no_sex$eas < 0 & 
                       no_sex$asj < 0 & no_sex$oth < 0 & no_sex$fin < 0, ]

  nfe.t <- table(nfe_snps$gene)
  nfe.t <- nfe.t[order(nfe.t, decreasing = TRUE)]
  nfe.k <- nfe.t[nfe.t < mutationburden.upper & nfe.t > mutationburden.lower]
  nfe.nsnps <- nfe_snps[nfe_snps$gene %in% names(nfe.k),]
  nfe.snps.k <- nfe.nsnps[(nfe.nsnps$nfe > lower.bound & nfe.nsnps$nfe < upper.bound),]
  
  #nfe.o <- nfe_snps[order(nfe_snps$nfe, decrnfeing = TRUE),]
  nfe.keep <- na.omit(nfe.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% nfe.keep]
  nfe.keep <- pf_top_varying_genes(temp.expr, keep.num)
  #nfe.keep <- nfe.keep[1:keep.num] # keep top most enhanced genes 
  
  dat <- c(dat,nfe.keep)
 
  # ########################################################################
  # ### AMR ###
  # #######################################################################
  amr_snps <- no_sex[no_sex$amr > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_onf < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 &
                      no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$eas_oea < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$fin < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
 
  amr.t <- table(amr_snps$gene)
  amr.t <- amr.t[order(amr.t, decreasing = TRUE)]
  amr.k <- amr.t[amr.t < mutationburden.upper & amr.t > mutationburden.lower]
  amr.nsnps <- amr_snps[amr_snps$gene %in% names(amr.k),]
  amr.snps.k <- amr.nsnps[(amr.nsnps$amr > lower.bound & amr.nsnps$amr < upper.bound),]
  
  #amr.o <- amr_snps[order(amr_snps$amr, decramring = TRUE),]
  
  amr.keep <- na.omit(amr.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% amr.keep]
  amr.keep <- pf_top_varying_genes(temp.expr, keep.num)
  #amr.keep <- amr.keep[1:keep.num] # keep top most enhanced genes 

  dat <- c(dat,amr.keep)
  
  
  # ########################################################################
  # ### FIN ###
  # #######################################################################
  fin_snps <- no_sex[no_sex$fin > 0 & no_sex$sas < 0 & no_sex$afr < 0 & no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$nfe_seu < 0 & no_sex$nfe < 0 & no_sex$nfe_bgr < 0 & no_sex$nfe_onf < 0 & no_sex$nfe_swe < 0 &
                       no_sex$nfe_nwe < 0 & no_sex$nfe_est < 0 & no_sex$eas_oea < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  fin.t <- table(fin_snps$gene)
  fin.t <- fin.t[order(fin.t, decreasing = TRUE)]
  fin.k <- fin.t[fin.t < mutationburden.upper & fin.t > mutationburden.lower]
  fin.nsnps <- fin_snps[fin_snps$gene %in% names(fin.k),]
  fin.snps.k <- fin.nsnps[(fin.nsnps$fin > lower.bound & fin.nsnps$fin < upper.bound),]
  
  #fin.o <- fin_snps[order(fin_snps$fin, decrfining = TRUE),]
  fin.keep <- na.omit(fin.snps.k$gene)# top most enhanced genes 
  temp.expr <- expression[,colnames(expression) %in% fin.keep]
  fin.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,fin.keep)

  # ########################################################################
  # ### EAS_OEA ###
  # #######################################################################
  eas_oea_snps <- no_sex[no_sex$eas_oea > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_onf < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 &
                       no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0  & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  eas_oea.t <- table(eas_oea_snps$gene)
  eas_oea.t <- eas_oea.t[order(eas_oea.t, decreasing = TRUE)]
  eas_oea.k <- eas_oea.t[eas_oea.t < mutationburden.upper & eas_oea.t > mutationburden.lower]
  eas_oea.nsnps <- eas_oea_snps[eas_oea_snps$gene %in% names(eas_oea.k),]
  eas_oea.snps.k <- eas_oea.nsnps[(eas_oea.nsnps$eas_oea > lower.bound & eas_oea.nsnps$eas_oea < upper.bound),]
  
  #eas_oea.o <- eas_oea_snps[order(eas_oea_snps$eas_oea, decreas_oeaing = TRUE),]
  eas_oea.keep <- na.omit(eas_oea.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% eas_oea.keep]
  eas_oea.keep <- pf_top_varying_genes(temp.expr, keep.num)
  dat <- c(dat,eas_oea.keep)
  
  # ########################################################################
  # ### NFE_NWE ###
  # #######################################################################
  nfe_nwe_snps <- no_sex[no_sex$nfe_nwe > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_onf < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  nfe_nwe.t <- table(nfe_nwe_snps$gene)
  nfe_nwe.t <- nfe_nwe.t[order(nfe_nwe.t, decreasing = TRUE)]
  nfe_nwe.k <- nfe_nwe.t[nfe_nwe.t < mutationburden.upper & nfe_nwe.t > mutationburden.lower]
  nfe_nwe.nsnps <- nfe_nwe_snps[nfe_nwe_snps$gene %in% names(nfe_nwe.k),]
  nfe_nwe.snps.k <- nfe_nwe.nsnps[(nfe_nwe.nsnps$nfe_nwe > lower.bound & nfe_nwe.nsnps$nfe_nwe < upper.bound),]
  
  #nfe_nwe.o <- nfe_nwe_snps[order(nfe_nwe_snps$nfe_nwe, decrnfe_nweing = TRUE),]
  nfe_nwe.keep <- na.omit(nfe_nwe.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% nfe_nwe.keep]
  nfe_nwe.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,nfe_nwe.keep)
 
  # ########################################################################
  # ### NFE_ONF ###
  # #######################################################################
  nfe_onf_snps <- no_sex[no_sex$nfe_onf > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  nfe_onf.t <- table(nfe_onf_snps$gene)
  nfe_onf.t <- nfe_onf.t[order(nfe_onf.t, decreasing = TRUE)]
  nfe_onf.k <- nfe_onf.t[nfe_onf.t < mutationburden.upper & nfe_onf.t > mutationburden.lower]
  nfe_onf.nsnps <- nfe_onf_snps[nfe_onf_snps$gene %in% names(nfe_onf.k),]
  nfe_onf.snps.k <- nfe_onf.nsnps[(nfe_onf.nsnps$nfe_onf > lower.bound & nfe_onf.nsnps$nfe_onf < upper.bound),]
  
  #nfe_onf.o <- nfe_onf_snps[order(nfe_onf_snps$nfe_onf, decrnfe_onfing = TRUE),]
  nfe_onf.keep <- na.omit(nfe_onf.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% nfe_onf.keep]
  nfe_onf.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,nfe_onf.keep)
  
  
  # ########################################################################
  # ### NFE_SEU ###
  # #######################################################################
  nfe_seu_snps <- no_sex[no_sex$nfe_seu > 0 & no_sex$nfe_onf < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  nfe_seu.t <- table(nfe_seu_snps$gene)
  nfe_seu.t <- nfe_seu.t[order(nfe_seu.t, decreasing = TRUE)]
  nfe_seu.k <- nfe_seu.t[nfe_seu.t < mutationburden.upper & nfe_seu.t > mutationburden.lower]
  nfe_seu.nsnps <- nfe_seu_snps[nfe_seu_snps$gene %in% names(nfe_seu.k),]
  nfe_seu.snps.k <- nfe_seu.nsnps[(nfe_seu.nsnps$nfe_seu > lower.bound & nfe_seu.nsnps$nfe_seu < upper.bound),]
  
  #nfe_seu.o <- nfe_seu_snps[order(nfe_seu_snps$nfe_seu, decrnfe_seuing = TRUE),]
  nfe_seu.keep <- na.omit(nfe_seu.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% nfe_seu.keep]
  nfe_seu.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,nfe_seu.keep)
  
  # ########################################################################
  # ### NFE_SWE ###
  # #######################################################################
  nfe_swe_snps <- no_sex[no_sex$nfe_swe > 0 & no_sex$nfe_onf < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  nfe_swe.t <- table(nfe_swe_snps$gene)
  nfe_swe.t <- nfe_swe.t[order(nfe_swe.t, decreasing = TRUE)]
  nfe_swe.k <- nfe_swe.t[nfe_swe.t < mutationburden.upper & nfe_swe.t > mutationburden.lower]
  nfe_swe.nsnps <- nfe_swe_snps[nfe_swe_snps$gene %in% names(nfe_swe.k),]
  nfe_swe.snps.k <- nfe_swe.nsnps[(nfe_swe.nsnps$nfe_swe > lower.bound & nfe_swe.nsnps$nfe_swe < upper.bound),]
  
  #nfe_swe.o <- nfe_swe_snps[order(nfe_swe_snps$nfe_swe, decrnfe_sweing = TRUE),]
  nfe_swe.keep <- na.omit(nfe_swe.snps.k$gene)#eep top most enhanced genes 
  temp.expr <- expression[,colnames(expression) %in% nfe_swe.keep]
  nfe_swe.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,nfe_swe.keep)
  
  # ########################################################################
  # ### SAS ###
  # #######################################################################
  sas_snps <- no_sex[no_sex$sas > 0 & no_sex$nfe_onf < 0 & no_sex$nfe_bgr < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  sas.t <- table(sas_snps$gene)
  sas.t <- sas.t[order(sas.t, decreasing = TRUE)]
  sas.k <- sas.t[sas.t < mutationburden.upper & sas.t > mutationburden.lower]
  sas.nsnps <- sas_snps[sas_snps$gene %in% names(sas.k),]
  sas.snps.k <- sas.nsnps[(sas.nsnps$sas > lower.bound & sas.nsnps$sas < upper.bound),]
  
  #sas.o <- sas_snps[order(sas_snps$sas, decrsasing = TRUE),]
  sas.keep <- na.omit(sas.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% sas.keep]
  sas.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,sas.keep)
  
  # ########################################################################
  # ### nfe_bgr ###
  # #######################################################################
  nfe_bgr_snps <- no_sex[no_sex$nfe_bgr > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                       no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  nfe_bgr.t <- table(nfe_bgr_snps$gene)
  nfe_bgr.t <- nfe_bgr.t[order(nfe_bgr.t, decreasing = TRUE)]
  nfe_bgr.k <- nfe_bgr.t[nfe_bgr.t < mutationburden.upper & nfe_bgr.t > mutationburden.lower]
  nfe_bgr.nsnps <- nfe_bgr_snps[nfe_bgr_snps$gene %in% names(nfe_bgr.k),]
  nfe_bgr.snps.k <- nfe_bgr.nsnps[(nfe_bgr.nsnps$nfe_bgr > lower.bound & nfe_bgr.nsnps$nfe_bgr < upper.bound),]
  
  #nfe_bgr.o <- nfe_bgr_snps[order(nfe_bgr_snps$nfe_bgr, decrnfe_bgring = TRUE),]
  nfe_bgr.keep <- na.omit(nfe_bgr.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% nfe_bgr.keep]
  nfe_bgr.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,nfe_bgr.keep)
  
  # ########################################################################
  # ### eas_jpn ###
  # #######################################################################
  eas_jpn_snps <- no_sex[no_sex$eas_jpn > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$nfe_bgr < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  eas_jpn.t <- table(eas_jpn_snps$gene)
  eas_jpn.t <- eas_jpn.t[order(eas_jpn.t, decreasing = TRUE)]
  eas_jpn.k <- eas_jpn.t[eas_jpn.t < mutationburden.upper & eas_jpn.t > mutationburden.lower]
  eas_jpn.nsnps <- eas_jpn_snps[eas_jpn_snps$gene %in% names(eas_jpn.k),]
  eas_jpn.snps.k <- eas_jpn.nsnps[(eas_jpn.nsnps$eas_jpn > lower.bound & eas_jpn.nsnps$eas_jpn < upper.bound),]
  
  #eas_jpn.o <- eas_jpn_snps[order(eas_jpn_snps$eas_jpn, decreas_jpning = TRUE),]
  eas_jpn.keep <- na.omit(eas_jpn.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% eas_jpn.keep]
  eas_jpn.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,eas_jpn.keep)



  
  # ########################################################################
  # ### eas_kor ###
  # #######################################################################
  eas_kor_snps <- no_sex[no_sex$eas_kor > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$nfe_bgr < 0 & no_sex$eas_jpn < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  
  eas_kor.t <- table(eas_kor_snps$gene)
  eas_kor.t <- eas_kor.t[order(eas_kor.t, decreasing = TRUE)]
  eas_kor.k <- eas_kor.t[eas_kor.t < mutationburden.upper & eas_kor.t > mutationburden.lower]
  eas_kor.nsnps <- eas_kor_snps[eas_kor_snps$gene %in% names(eas_kor.k),]
  eas_kor.snps.k <- eas_kor.nsnps[(eas_kor.nsnps$eas_kor > lower.bound & eas_kor.nsnps$eas_kor < upper.bound),]
  
  #eas_kor.o <- eas_kor_snps[order(eas_kor_snps$eas_kor, decreas_koring = TRUE),]
  eas_kor.keep <- na.omit(eas_kor.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% eas_kor.keep]
  eas_kor.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,eas_kor.keep)
  
  
  # ########################################################################
  # ### nfe_est ###
  # #######################################################################
  nfe_est_snps <- no_sex[no_sex$nfe_est > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$nfe_bgr < 0 & no_sex$eas_jpn < 0 & no_sex$fin < 0 & no_sex$eas_kor < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  nfe_est.t <- table(nfe_est_snps$gene)
  nfe_est.t <- nfe_est.t[order(nfe_est.t, decreasing = TRUE)]
  nfe_est.k <- nfe_est.t[nfe_est.t < mutationburden.upper & nfe_est.t > mutationburden.lower]
  nfe_est.nsnps <- nfe_est_snps[nfe_est_snps$gene %in% names(nfe_est.k),]
  nfe_est.snps.k <- nfe_est.nsnps[(nfe_est.nsnps$nfe_est > lower.bound & nfe_est.nsnps$nfe_est < upper.bound),]
  nfe_est.keep <- na.omit(nfe_est.snps.k$gene)
  
  temp.expr <- expression[,colnames(expression) %in% nfe_est.keep]
  nfe_est.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,nfe_est.keep)
  

  # ########################################################################
  # ### asj ###
  # #######################################################################
  asj_snps <- no_sex[no_sex$asj > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$nfe_bgr < 0 & no_sex$eas_jpn < 0 & no_sex$fin < 0 & no_sex$eas_kor < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$nfe_est < 0 & no_sex$oth < 0, ]
  asj.t <- table(asj_snps$gene)
  asj.t <- asj.t[order(asj.t, decreasing = TRUE)]
  asj.k <- asj.t[asj.t < mutationburden.upper & asj.t > mutationburden.lower]
  asj.nsnps <- asj_snps[asj_snps$gene %in% names(asj.k),]
  asj.snps.k <- asj.nsnps[(asj.nsnps$asj > lower.bound & asj.nsnps$asj < upper.bound),]

  asj.keep <- na.omit(asj.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% asj.keep]
  asj.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,asj.keep)
 
  
  # ########################################################################
  # ### oth ###
  # #######################################################################
  oth_snps <- no_sex[no_sex$oth > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                       no_sex$nfe_bgr < 0 & no_sex$eas_jpn < 0 & no_sex$fin < 0 & no_sex$eas_kor < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$nfe_est < 0 & no_sex$asj < 0, ]

  oth.t <- table(oth_snps$gene)
  oth.t <- oth.t[order(oth.t, decreasing = TRUE)]
  oth.k <- oth.t[oth.t < mutationburden.upper & oth.t > mutationburden.lower]
  oth.nsnps <- oth_snps[oth_snps$gene %in% names(oth.k),]
  oth.snps.k <- oth.nsnps[(oth.nsnps$oth > lower.bound & oth.nsnps$oth < upper.bound),]
  
  #oth.o <- oth_snps[order(oth_snps$oth, decrothing = TRUE),]
  oth.keep <- na.omit(oth.snps.k$gene)
  temp.expr <- expression[,colnames(expression) %in% oth.keep]
  oth.keep <- pf_top_varying_genes(temp.expr, keep.num)
  
  dat <- c(dat,oth.keep)
  dat <- dat
  dat <- dat[dat %in% "subtype" == FALSE] # get rid of subtype 
  valid.keep <- unique(dat)
  # valid.genes <- table(dat)
  # valid.keep <- valid.genes[valid.genes > 8] # chose 8 becauses there are 8 european ancestries.
  # valid.keep <-  names(valid.keep)
  
  ## table 
  
  return(valid.keep)
  }
