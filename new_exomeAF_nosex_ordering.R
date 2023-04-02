########################################################################
### EXOME ALLELE FREQUENCY ORDERING ###
#######################################################################

require(tibble)
require(readr)  # for read_csv()
require(dplyr) #data frame handling 
require(tidyr)
require(stringi)

#######################################################################
### FUNCITON TO GET GENES IN NETWORK AND SORT EACH ANCESTRY IN DECREASING ORDER BY AF###
#######################################################################

order_frequencies <- function(include_genes, exome_file){
  #### keep only include genes that are relevant in network ####
  no_sex <- exome_file[exome_file$gene %in% include_genes,] # genes from the human base network 
  
  #### dataframe to hold the genes selected for each ancestry ####
  dat <- list()
  
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

  #afr_snps <- no_sex[no_sex$afr > 0.045,]
  
  #temp <- mtcars[rowSums(mtcars[, names(mtcars) != "afr"]) < 0,]
  # afr.t <- table(afr.dat$afr)
  # afr.t.s <- afr.t[order(afr.t, decreasing = TRUE)]
  # one <- as.numeric(names(afr.t.s[1]))
  # two <- as.numeric(names(afr.t.s[2]))
  # #afr_snps <- afr.dat[(afr.dat$afr == one),]
  # if(one > two){
  #   afr_snps <- afr.dat[(afr.dat$afr >= two & afr.dat$afr <= one),]
  # }else{
  #   afr_snps <- afr.dat[(afr.dat$afr >= one & afr.dat$afr <= two),]
  # }

  dat <- c(dat,unique(afr_snps$gene))
  
  ########################################################################
  ### EAS ###
  #######################################################################
  eas_snps <- no_sex[no_sex$eas > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_onf < 0 & no_sex$amr < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 &
                       no_sex$nfe_est < 0 & no_sex$nfe < 0 & no_sex$fin < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #eas_snps <- no_sex[no_sex$eas > 0.045,]
   # eas.dat <- no_sex[no_sex$eas > 0,]
  # eas.t <- table(eas.dat$eas)
  # eas.t.s <- eas.t[order(eas.t, decreasing = TRUE)]
  # one <- as.numeric(names(eas.t.s[1]))
  # two <- as.numeric(names(eas.t.s[2]))
  # #eas_snps <- eas.dat[(eas.dat$eas == one),]
  # if(one > two){
  #   eas_snps <- eas.dat[(eas.dat$eas >= two & eas.dat$eas <= one),]
  #   #print("here")
  # }else{
  #   eas_snps <- eas.dat[(eas.dat$eas >= one & eas.dat$eas <= two),]
  # }
  dat <- c(dat,unique(eas_snps$gene))
  #
  #
  # ########################################################################
  # ### NFE ###
  # #######################################################################
  nfe_snps <- no_sex[no_sex$nfe > 0 & no_sex$sas < 0 & no_sex$amr < 0 & no_sex$afr < 0 & no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$eas_oea < 0 & no_sex$eas < 0 & 
                       no_sex$asj < 0 & no_sex$oth < 0 & no_sex$fin < 0, ]

  #nfe_snps <- no_sex[no_sex$nfe > 0.045,]
  # nfe.dat <- no_sex[no_sex$nfe > 0,]
  # nfe.t <- table(nfe.dat$nfe)
  # nfe.t.s <- nfe.t[order(nfe.t, decreasing = TRUE)]
  # one <- as.numeric(names(nfe.t.s[1]))
  # two <- as.numeric(names(nfe.t.s[2]))
  # #nfe_snps <- nfe.dat[(nfe.dat$nfe == one),]
  # if(one > two){
  #   nfe_snps <- nfe.dat[(nfe.dat$nfe >= two & nfe.dat$nfe <= one),]
  #   print("here")
  # }else{
  #   nfe_snps <- nfe.dat[(nfe.dat$nfe >= one & nfe.dat$nfe <= two),]
  # }
  dat <- c(dat,unique(nfe_snps$gene))

  # ########################################################################
  # ### AMR ###
  # #######################################################################
  amr_snps <- no_sex[no_sex$amr > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_onf < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 &
                      no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$eas_oea < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$fin < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #amr_snps <- no_sex[no_sex$amr > 0.045,]
  # amr.dat <- no_sex[no_sex$amr > 0,]
  # amr.t <- table(amr.dat$amr)
  # amr.t.s <- amr.t[order(amr.t, decreasing = TRUE)]
  # one <- as.numeric(names(amr.t.s[1]))
  # two <- as.numeric(names(amr.t.s[2]))
  # #amr_snps <- amr.dat[(amr.dat$amr == one),]
  # if(one > two){
  #   amr_snps <- amr.dat[(amr.dat$amr >= two & amr.dat$amr <= one),]
  #   print("here")
  # }else{
  #   amr_snps <- amr.dat[(amr.dat$amr >= one & amr.dat$amr <= two),]
  # }
  dat <- c(dat,unique(amr_snps$gene))
  
  # ########################################################################
  # ### FIN ###
  # #######################################################################
  fin_snps <- no_sex[no_sex$fin > 0 & no_sex$sas < 0 & no_sex$afr < 0 & no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$nfe_seu < 0 & no_sex$nfe < 0 & no_sex$nfe_bgr < 0 & no_sex$nfe_onf < 0 & no_sex$nfe_swe < 0 &
                       no_sex$nfe_nwe < 0 & no_sex$nfe_est < 0 & no_sex$eas_oea < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]


  #fin_snps <- no_sex[no_sex$fin > 0.045,]
  # fin.dat <- no_sex[no_sex$fin > 0,]
  # fin.t <- table(fin.dat$fin)
  # fin.t.s <- fin.t[order(fin.t, decreasing = TRUE)]
  # one <- as.numeric(names(fin.t.s[1]))
  # two <- as.numeric(names(fin.t.s[2]))
  # #fin_snps <- fin.dat[(fin.dat$fin == one),]
  # if(one > two){
  #   fin_snps <- fin.dat[(fin.dat$fin >= two & fin.dat$fin <= one),]
  #   print("here")
  # }else{
  #   fin_snps <- fin.dat[(fin.dat$fin >= one & fin.dat$fin <= two),]
  #}
  dat <- c(dat,unique(fin_snps$gene))
  
  # ########################################################################
  # ### EAS_OEA ###
  # #######################################################################
  eas_oea_snps <- no_sex[no_sex$eas_oea > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_onf < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 &
                       no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0  & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #eas_oea_snps <- no_sex[no_sex$eas_oea > 0.045,]
  # eas_oea.dat <- no_sex[no_sex$eas_oea > 0,]
  # eas_oea.t <- table(eas_oea.dat$eas_oea)
  # eas_oea.t.s <- eas_oea.t[order(eas_oea.t, decreasing = TRUE)]
  # one <- as.numeric(names(eas_oea.t.s[1]))
  # two <- as.numeric(names(eas_oea.t.s[2]))
  # #eas_oea_snps <- eas_oea.dat[(eas_oea.dat$eas_oea == one),]
  # if(one > two){
  #   eas_oea_snps <- eas_oea.dat[(eas_oea.dat$eas_oea >= two & eas_oea.dat$eas_oea <= one),]
  #   print("here")
  # }else{
  #   eas_oea_snps <- eas_oea.dat[(eas_oea.dat$eas_oea >= one & eas_oea.dat$eas_oea <= two),]
  # }
  dat <- c(dat,unique(eas_oea_snps$gene))
  
  # ########################################################################
  # ### NFE_NWE ###
  # #######################################################################
  nfe_nwe_snps <- no_sex[no_sex$nfe_nwe > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_onf < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #nfe_nwe_snps <- no_sex[no_sex$nfe_nwe > 0.045,]
  # nfe_nwe.dat <- no_sex[no_sex$nfe_nwe > 0,]
  # nfe_nwe.t <- table(nfe_nwe.dat$nfe_nwe)
  # nfe_nwe.t.s <- nfe_nwe.t[order(nfe_nwe.t, decreasing = TRUE)]
  # one <- as.numeric(names(nfe_nwe.t.s[1]))
  # two <- as.numeric(names(nfe_nwe.t.s[2]))
  # #nfe_nwe_snps <- nfe_nwe.dat[(nfe_nwe.dat$nfe_nwe == one),]
  # if(one > two){
  #   nfe_nwe_snps <- nfe_nwe.dat[(nfe_nwe.dat$nfe_nwe >= two & nfe_nwe.dat$nfe_nwe <= one),]
  #   print("here")
  # }else{
  #   nfe_nwe_snps <- nfe_nwe.dat[(nfe_nwe.dat$nfe_nwe >= one & nfe_nwe.dat$nfe_nwe <= two),]
  # }
  dat <- c(dat,unique(nfe_nwe_snps$gene))
  
  # ########################################################################
  # ### NFE_ONF ###
  # #######################################################################
  nfe_onf_snps <- no_sex[no_sex$nfe_onf > 0 & no_sex$nfe_seu < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #nfe_onf_snps <- no_sex[no_sex$nfe_onf > 0.045,]
  # nfe_onf.dat <- no_sex[no_sex$nfe_onf > 0,]
  # nfe_onf.t <- table(nfe_onf.dat$nfe_onf)
  # nfe_onf.t.s <- nfe_onf.t[order(nfe_onf.t, decreasing = TRUE)]
  # one <- as.numeric(names(nfe_onf.t.s[1]))
  # two <- as.numeric(names(nfe_onf.t.s[2]))
  # #nfe_onf_snps <- nfe_onf.dat[(nfe_onf.dat$nfe_onf == one),]
  # if(one > two){
  #   nfe_onf_snps <- nfe_onf.dat[(nfe_onf.dat$nfe_onf >= two & nfe_onf.dat$nfe_onf <= one),]
  #   print("here")
  # }else{
  #   nfe_onf_snps <- nfe_onf.dat[(nfe_onf.dat$nfe_onf >= one & nfe_onf.dat$nfe_onf <= two),]
  # }
  dat <- c(dat,unique(nfe_onf_snps$gene))
  
  # ########################################################################
  # ### NFE_SEU ###
  # #######################################################################
  nfe_seu_snps <- no_sex[no_sex$nfe_seu > 0 & no_sex$nfe_onf < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_swe < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #nfe_seu_snps <- no_sex[no_sex$nfe_seu > 0.045,]
  # nfe_seu.dat <- no_sex[no_sex$nfe_seu > 0,]
  # nfe_seu.t <- table(nfe_seu.dat$nfe_seu)
  # nfe_seu.t.s <- nfe_seu.t[order(nfe_seu.t, decreasing = TRUE)]
  # one <- as.numeric(names(nfe_seu.t.s[1]))
  # two <- as.numeric(names(nfe_seu.t.s[2]))
  # #nfe_seu_snps <- nfe_seu.dat[(nfe_seu.dat$nfe_seu == one),]
  # if(one > two){
  #   nfe_seu_snps <- nfe_seu.dat[(nfe_seu.dat$nfe_seu >= two & nfe_seu.dat$nfe_seu <= one),]
  #   print("here")
  # }else{
  #   nfe_seu_snps <- nfe_seu.dat[(nfe_seu.dat$nfe_seu >= one & nfe_seu.dat$nfe_seu <= two),]
  # }
  dat <- c(dat,unique(nfe_seu_snps$gene))
  # ########################################################################
  # ### NFE_SWE ###
  # #######################################################################
  nfe_swe_snps <- no_sex[no_sex$nfe_swe > 0 & no_sex$nfe_onf < 0 & no_sex$nfe_bgr < 0 & no_sex$sas < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #nfe_swe_snps <- no_sex[no_sex$nfe_swe > 0.045,]
  # nfe_swe.dat <- no_sex[no_sex$nfe_swe > 0,]
  # nfe_swe.t <- table(nfe_swe.dat$nfe_swe)
  # nfe_swe.t.s <- nfe_swe.t[order(nfe_swe.t, decreasing = TRUE)]
  # one <- as.numeric(names(nfe_swe.t.s[1]))
  # two <- as.numeric(names(nfe_swe.t.s[2]))
  # #nfe_swe_snps <- nfe_swe.dat[(nfe_swe.dat$nfe_swe == one),]
  # if(one > two){
  #   nfe_swe_snps <- nfe_swe.dat[(nfe_swe.dat$nfe_swe >= two & nfe_swe.dat$nfe_swe <= one),]
  #   print("here")
  # }else{
  #   nfe_swe_snps <- nfe_swe.dat[(nfe_swe.dat$nfe_swe >= one & nfe_swe.dat$nfe_swe <= two),]
  # }
  dat <- c(dat,unique(nfe_swe_snps$gene))
  
  # ########################################################################
  # ### SAS ###
  # #######################################################################
  sas_snps <- no_sex[no_sex$sas > 0 & no_sex$nfe_onf < 0 & no_sex$nfe_bgr < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #sas_snps <- no_sex[no_sex$sas > 0.045,]
  # sas.dat <- no_sex[no_sex$sas > 0,]
  # sas.t <- table(sas.dat$sas)
  # sas.t.s <- sas.t[order(sas.t, decreasing = TRUE)]
  # one <- as.numeric(names(sas.t.s[1]))
  # two <- as.numeric(names(sas.t.s[2]))
  # #sas_snps <- sas.dat[(sas.dat$sas == one),]
  # if(one > two){
  #   sas_snps <- sas.dat[(sas.dat$sas >= two & sas.dat$sas <= one),]
  #   print("here")
  # }else{
  #   sas_snps <- sas.dat[(sas.dat$sas >= one & sas.dat$sas <= two),]
  # }
  dat <- c(dat,unique(sas_snps$gene))
 
  
  # ########################################################################
  # ### nfe_bgr ###
  # #######################################################################
  nfe_bgr_snps <- no_sex[no_sex$nfe_bgr > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                       no_sex$eas_jpn < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #nfe_bgr_snps <- no_sex[no_sex$nfe_bgr > 0.045,]
  
  # nfe_bgr.dat <- no_sex[no_sex$nfe_bgr > 0,]
  # nfe_bgr.t <- table(nfe_bgr.dat$nfe_bgr)
  # nfe_bgr.t.s <- nfe_bgr.t[order(nfe_bgr.t, decreasing = TRUE)]
  # one <- as.numeric(names(nfe_bgr.t.s[1]))
  # two <- as.numeric(names(nfe_bgr.t.s[2]))
  # #nfe_bgr_snps <- nfe_bgr.dat[(nfe_bgr.dat$nfe_bgr == one),]
  # if(one > two){
  #   nfe_bgr_snps <- nfe_bgr.dat[(nfe_bgr.dat$nfe_bgr >= two & nfe_bgr.dat$nfe_bgr <= one),]
  #   print("here")
  # }else{
  #   nfe_bgr_snps <- nfe_bgr.dat[(nfe_bgr.dat$nfe_bgr >= one & nfe_bgr.dat$nfe_bgr <= two),]
  # }
  dat <- c(dat,unique(nfe_bgr_snps$gene))

  
  # ########################################################################
  # ### eas_jpn ###
  # #######################################################################
  eas_jpn_snps <- no_sex[no_sex$eas_jpn > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$nfe_bgr < 0 & no_sex$eas_kor < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #eas_jpn_snps <- no_sex[no_sex$eas_jpn > 0.045,]
  # eas_jpn.dat <- no_sex[no_sex$eas_jpn > 0,]
  # eas_jpn.t <- table(eas_jpn.dat$eas_jpn)
  # eas_jpn.t.s <- eas_jpn.t[order(eas_jpn.t, decreasing = TRUE)]
  # one <- as.numeric(names(eas_jpn.t.s[1]))
  # two <- as.numeric(names(eas_jpn.t.s[2]))
  # #eas_jpn_snps <- eas_jpn.dat[(eas_jpn.dat$eas_jpn == one),]
  # if(one > two){
  #   eas_jpn_snps <- eas_jpn.dat[(eas_jpn.dat$eas_jpn >= two & eas_jpn.dat$eas_jpn <= one),]
  #   print("here")
  # }else{
  #   eas_jpn_snps <- eas_jpn.dat[(eas_jpn.dat$eas_jpn >= one & eas_jpn.dat$eas_jpn <= two),]
  # }
  dat <- c(dat,unique(eas_jpn_snps$gene))

  
  
  # ########################################################################
  # ### eas_kor ###
  # #######################################################################
  eas_kor_snps <- no_sex[no_sex$eas_kor > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$nfe_bgr < 0 & no_sex$eas_jpn < 0 & no_sex$fin < 0 & no_sex$nfe_est < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]
  #
  #eas_kor_snps <- no_sex[no_sex$eas_kor > 0.045,]
  # eas_kor.dat <- no_sex[no_sex$eas_kor > 0,]
  # eas_kor.t <- table(eas_kor.dat$eas_kor)
  # eas_kor.t.s <- eas_kor.t[order(eas_kor.t, decreasing = TRUE)]
  # one <- as.numeric(names(eas_kor.t.s[1]))
  # two <- as.numeric(names(eas_kor.t.s[2]))
  # #eas_kor_snps <- eas_kor.dat[(eas_kor.dat$eas_kor == one),]
  # if(one > two){
  #   eas_kor_snps <- eas_kor.dat[(eas_kor.dat$eas_kor >= two & eas_kor.dat$eas_kor <= one),]
  #   print("here")
  # }else{
  #   eas_kor_snps <- eas_kor.dat[(eas_kor.dat$eas_kor >= one & eas_kor.dat$eas_kor <= two),]
  # }
  dat <- c(dat,unique(eas_kor_snps$gene))
  
  
  # ########################################################################
  # ### nfe_est ###
  # #######################################################################
  nfe_est_snps <- no_sex[no_sex$nfe_est > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$nfe_bgr < 0 & no_sex$eas_jpn < 0 & no_sex$fin < 0 & no_sex$eas_kor < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0 & no_sex$oth < 0, ]

  #nfe_est_snps <- no_sex[no_sex$nfe_est > 0.045,]
  # nfe_est.dat <- no_sex[no_sex$nfe_est > 0,]
  # nfe_est.t <- table(nfe_est.dat$nfe_est)
  # nfe_est.t.s <- nfe_est.t[order(nfe_est.t, decreasing = TRUE)]
  # one <- as.numeric(names(nfe_est.t.s[1]))
  # two <- as.numeric(names(nfe_est.t.s[2]))
  # #nfe_est_snps <- nfe_est.dat[(nfe_est.dat$nfe_est == one),]
  # if(one > two){
  #   nfe_est_snps <- nfe_est.dat[(nfe_est.dat$nfe_est >= two & nfe_est.dat$nfe_est <= one),]
  #   print("here")
  # }else{
  #   nfe_est_snps <- nfe_est.dat[(nfe_est.dat$nfe_est >= one & nfe_est.dat$nfe_est <= two),]
  # }
  dat <- c(dat,unique(nfe_est_snps$gene))
  

  # ########################################################################
  # ### asj ###
  # #######################################################################
  asj_snps <- no_sex[no_sex$asj > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                           no_sex$nfe_bgr < 0 & no_sex$eas_jpn < 0 & no_sex$fin < 0 & no_sex$eas_kor < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$nfe_est < 0 & no_sex$oth < 0, ]

  #asj_snps <- no_sex[no_sex$asj > 0.045,]
  # asj.dat <- no_sex[no_sex$asj > 0,]
  # asj.t <- table(asj.dat$asj)
  # asj.t.s <- asj.t[order(asj.t, decreasing = TRUE)]
  # one <- as.numeric(names(asj.t.s[1]))
  # two <- as.numeric(names(asj.t.s[2]))
  # #asj_snps <- asj.dat[(asj.dat$asj == one),]
  # if(one > two){
  #   asj_snps <- asj.dat[(asj.dat$asj >= two & asj.dat$asj <= one),]
  #   print("here")
  # }else{
  #   asj_snps <- asj.dat[(asj.dat$asj >= one & asj.dat$asj <= two),]
  # }
  dat <- c(dat,unique(asj_snps$gene))
  
  
  # ########################################################################
  # ### oth ###
  # #######################################################################
  oth_snps <- no_sex[no_sex$oth > 0 & no_sex$nfe_onf < 0 & no_sex$sas < 0 & no_sex$nfe_swe < 0 & no_sex$nfe_nwe < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$nfe_seu < 0 & no_sex$eas_oea < 0 &
                       no_sex$nfe_bgr < 0 & no_sex$eas_jpn < 0 & no_sex$fin < 0 & no_sex$eas_kor < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$nfe_est < 0 & no_sex$asj < 0, ]

  #oth_snps <- no_sex[no_sex$oth > 0.045,]
  # oth.dat <- no_sex[no_sex$oth > 0,]
  # oth.t <- table(oth.dat$oth)
  # oth.t.s <- oth.t[order(oth.t, decreasing = TRUE)]
  # one <- as.numeric(names(oth.t.s[1]))
  # two <- as.numeric(names(oth.t.s[2]))
  # #oth_snps <- oth.dat[(oth.dat$oth == one),]
  # if(one > two){
  #   oth_snps <- oth.dat[(oth.dat$oth >= two & oth.dat$oth <= one),]
  #   print("here")
  # }else{
  #   oth_snps <- oth.dat[(oth.dat$oth >= one & oth.dat$oth <= two),]
  # }
  dat <- c(dat,unique(oth_snps$gene))
  dat <- unique(dat)
  
  
  return(dat)
  }

# ########################################################################
# ### NFE_SEU ###
# #######################################################################
# #--- get top n unique genes and associated snps ---# 
# nfe_seu_desc <- nfe_seu_desc[!(is.na(nfe_seu_desc$nfe_seu)),]
# #unique_genes <- unique(nfe_seu_desc$gene)
# #top <- unique_genes[1:gene_num]
# #nfe_seu_genes <- tibble(top)
# #colnames(nfe_seu_genes) <- "genes"
# dat <- rbind(dat, nfe_seu_genes)
# 
# #--- get allele frequency of first snp from last gene to get minimum allele frequency for ancestry ---#
# top_snps <- nfe_seu_desc[nfe_seu_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),] #match return a vector of the positions of (first) matches of its first argument in its second. 
# last_row <- nrow(t.first) #get first snp allele frequency occurence of last gene snp
# min <- t.first[last_row,]$nfe_seu
# min_af.dat[1,1]<-"nfe_seu"
# min_af.dat[1,2] <- min
# print(paste("Minimum allele frequency for nfe_seu is ", min))
# nfe_seu_snps <- nfe_seu_desc[nfe_seu_desc$nfe_seu >= min,]
# write_delim(nfe_seu_snps, file = paste0(base_path, nfe_seu_out), delim = '\t', col_names = TRUE)
# 
# 
# ########################################################################
# ### NFE_BGR ###
# #######################################################################
# #--- get top n unique genes and associated snps ---# 
# nfe_bgr_desc <- nfe_bgr_desc[!(is.na(nfe_bgr_desc$nfe_bgr)),]
# unique_genes <- unique(nfe_bgr_desc$gene)
# top <- unique_genes[1:gene_num]
# nfe_bgr_genes <- tibble(top)
# colnames(nfe_bgr_genes) <- "genes"
# dat <- rbind(dat, nfe_bgr_genes)
# 
# #--- get allele frequency of first snp from last gene to get minimum allele frequency for ancestry ---#
# top_snps <- nfe_bgr_desc[nfe_bgr_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$nfe_bgr
# min_af.dat[2,1]<-"nfe_bgr"
# min_af.dat[2,2] <- min
# print(paste("Minimum allele frequency for nfe_bgr is ", min))
# nfe_bgr_snps <- nfe_bgr_desc[nfe_bgr_desc$nfe_bgr >= min,]
# write_delim(nfe_bgr_snps, file = paste0(base_path, nfe_bgr_out), 
#             delim = '\t', col_names = TRUE)
# 
# 
# 
# ########################################################################
# ### SAS ###
# #######################################################################
# #--- get top n unique genes and associated snps ---# 
# sas_desc <- sas_desc[!(is.na(sas_desc$sas)),]
# unique_genes <- unique(sas_desc$gene)
# top <- unique_genes[1:gene_num]
# sas_genes <- tibble(top)
# colnames(sas_genes) <- "genes"
# dat <- rbind(dat, sas_genes)
# 
# #--- get allele frequency of first snp from last gene to get minimum allele frequency for ancestry ---#
# top_snps <- sas_desc[sas_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$sas
# min_af.dat[4,1]<-"sas"
# min_af.dat[4,2] <- min
# print(paste("Minimum allele frequency for sas is ", min))
# sas_snps <- sas_desc[sas_desc$sas >= min,]
# write_delim(sas_snps, file = paste0(base_path, sas_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### NFE_ONF ###
# #######################################################################
# #take out NA because they get placed in the middle when sorted 
# nfe_onf_desc <- nfe_onf_desc[!(is.na(nfe_onf_desc$nfe_onf)),]
# unique_genes <- unique(nfe_onf_desc$gene)
# top <- unique_genes[1:gene_num]
# nfe_onf_genes <- tibble(top)
# colnames(nfe_onf_genes) <- "genes"
# dat <- rbind(dat, nfe_onf_genes)
# 
# #--- get allele frequency of first snp from last gene to get minimum allele frequency for ancestry ---#
# top_snps <- nfe_onf_desc[nfe_onf_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$nfe_onf
# min_af.dat[5,1]<-"nfe_onf"
# min_af.dat[5,2] <- min
# print(paste("Minimum allele frequency for nfe_onf is ", min))
# nfe_onf_snps <- nfe_onf_desc[nfe_onf_desc$nfe_onf >= min,]
# write_delim(nfe_onf_snps, file = paste0(base_path, nfe_onf_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### AMR ###
# #######################################################################
# amr_desc <- amr_desc[!(is.na(amr_desc$amr)),]
# unique_genes <- unique(amr_desc$gene)
# top <- unique_genes[1:gene_num]
# amr_genes <- tibble(top)
# colnames(amr_genes) <- "genes"
# dat <- rbind(dat, amr_genes)
# 
# #--- get allele frequency of first snp from last gene to get minimum allele frequency for ancestry ---#
# top_snps <- amr_desc[amr_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$amr
# min_af.dat[6,1]<-"amr"
# min_af.dat[6,2] <- min
# print(paste("Minimum allele frequency for amr is ", min))
# amr_snps <- amr_desc[amr_desc$amr >= min,]
# #write snps to file 
# write_delim(amr_snps, file = paste0(base_path,amr_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### EAS ###
# #######################################################################
# eas_desc <- eas_desc[!(is.na(eas_desc$eas)),]
# 
# eas_snps <- eas_desc[eas_desc$eas > 0.00005,]
# unique_genes <- eas_snps$gene
# unique_genes <- unique(eas_snps$gene)
# dat <- rbind(dat, unique_genes)
# 
# 
# 
# unique_genes <- unique(eas_desc$gene)
# top <- unique_genes[1:gene_num]
# eas_genes <- tibble(top)
# colnames(eas_genes) <- "genes"
# dat <- rbind(dat, eas_genes)
# 
# #--- get allele frequency of first snp from last gene to get minimum allele frequency for ancestry ---#
# top_snps <- eas_desc[eas_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$eas
# min_af.dat[7,1]<-"eas"
# min_af.dat[7,2] <- min
# print(paste("Minimum allele frequency for eas is ", min))
# eas_snps <- eas_desc[eas_desc$eas >= min,]
# write_delim(eas_snps, file = paste0(base_path,eas_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### NFE_SWE ###
# #######################################################################
# nfe_swe_desc <- nfe_swe_desc[!(is.na(nfe_swe_desc$nfe_swe)),]
# unique_genes <- unique(nfe_swe_desc$gene)
# top <- unique_genes[1:gene_num]
# nfe_swe_genes <- tibble(top)
# colnames(nfe_swe_genes) <- "genes"
# dat <- rbind(dat, nfe_swe_genes)
# 
# #--- get allele frequency of first snp from last gene to get minimum allele frequency for ancestry ---#
# top_snps <- nfe_swe_desc[nfe_swe_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$nfe_swe
# min_af.dat[8,1]<-"nfe_swe"
# min_af.dat[8,2] <- min
# print(paste("Minimum allele frequency for nfe_swe is ", min))
# nfe_swe_snps <- nfe_swe_desc[nfe_swe_desc$nfe_swe >= min,]
# #write snps to file 
# write_delim(nfe_swe_snps, file = paste0(base_path, nfe_swe_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### NFE_NWE ###
# #######################################################################
# nfe_nwe_desc <- nfe_nwe_desc[!(is.na(nfe_nwe_desc$nfe_nwe)),]
# unique_genes <- unique(nfe_nwe_desc$gene)
# top <- unique_genes[1:gene_num]
# nfe_nwe_genes <- tibble(top)
# colnames(nfe_nwe_genes) <- "genes"
# dat <- rbind(dat, nfe_nwe_genes)
# 
# top_snps <- nfe_nwe_desc[nfe_nwe_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$nfe_nwe
# min_af.dat[9,1]<-"nfe_nwe"
# min_af.dat[9,2] <- min
# print(paste("Minimum allele frequency for nfe_nwe is ", min))
# nfe_nwe_snps <- nfe_nwe_desc[nfe_nwe_desc$nfe_nwe >= min,]
# write_delim(nfe_nwe_snps, file = paste0(base_path, nfe_nwe_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### EAS_JPN ###
# #######################################################################
# eas_jpn_desc <- eas_jpn_desc[!(is.na(eas_jpn_desc$eas_jpn)),]
# unique_genes <- unique(eas_jpn_desc$gene)
# top <- unique_genes[1:gene_num]
# eas_jpn_genes <- tibble(top)
# colnames(eas_jpn_genes) <- "genes"
# dat <- rbind(dat, eas_jpn_genes)
# 
# top_snps <- eas_jpn_desc[eas_jpn_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$eas_jpn
# min_af.dat[10,1]<-"eas_jpn"
# min_af.dat[10,2] <- min
# print(paste("Minimum allele frequency for eas_jpn is ", min))
# eas_jpn_snps <- eas_jpn_desc[eas_jpn_desc$eas_jpn >= min,]
# write_delim(eas_jpn_snps, file = paste0(base_path,eas_jpn_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### EAS_KOR ###
# #######################################################################
# eas_kor_desc <- eas_kor_desc[!(is.na(eas_kor_desc$eas_kor)),]
# unique_genes <- unique(eas_kor_desc$gene)
# top <- unique_genes[1:gene_num]
# eas_kor_genes <- tibble(top)
# colnames(eas_kor_genes) <- "genes"
# dat <- rbind(dat, eas_kor_genes)
# 
# top_snps <- eas_kor_desc[eas_kor_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$eas_kor
# min_af.dat[11,1]<-"eas_kor"
# min_af.dat[11,2] <- min
# print(paste("Minimum allele frequency for eas_kor is ", min))
# eas_kor_snps <- eas_kor_desc[eas_kor_desc$eas_kor >= min,]
# write_delim(eas_kor_snps, file = paste0(base_path,eas_kor_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### EAS_OEA ###
# #######################################################################
# eas_oea_desc <- eas_oea_desc[!(is.na(eas_oea_desc$eas_oea)),]
# unique_genes <- unique(eas_oea_desc$gene)
# top <- unique_genes[1:gene_num]
# eas_oea_genes <- tibble(top)
# colnames(eas_oea_genes) <- "genes"
# dat <- rbind(dat, eas_oea_genes)
# 
# top_snps <- eas_oea_desc[eas_oea_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$eas_oea
# min_af.dat[12,1]<-"eas_oea"
# min_af.dat[12,2] <- min
# print(paste("Minimum allele frequency for eas_oea is ", min))
# eas_oea_snps <- eas_oea_desc[eas_oea_desc$eas_oea >= min,]
# write_delim(eas_oea_snps, file = paste0(base_path,eas_oea_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### NFE_EST ###
# #######################################################################
# nfe_est_desc <- nfe_est_desc[!(is.na(nfe_est_desc$nfe_est)),]
# unique_genes <- unique(nfe_est_desc$gene)
# top <- unique_genes[1:gene_num]
# nfe_est_genes <- tibble(top)
# colnames(nfe_est_genes) <- "genes"
# dat <- rbind(dat, nfe_est_genes)
# 
# top_snps <- nfe_est_desc[nfe_est_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$nfe_est
# min_af.dat[13,1]<-"nfe_est"
# min_af.dat[13,2] <- min
# print(paste("Minimum allele frequency for nfe_est is ", min))
# nfe_est_snps <- nfe_est_desc[nfe_est_desc$nfe_est >= min,]
# #write snps to file 
# write_delim(nfe_est_snps, file = paste0(base_path, nfe_est_out), 
#             delim = '\t', col_names = TRUE)
# 
# 
# ########################################################################
# ### NFE ###
# #######################################################################
# nfe_desc  <- nfe_desc[!(is.na(nfe_desc$nfe)),]
# 
# nfe_snps <- nfe_desc[nfe_desc$nfe > 0.00005,]
# unique_genes <- nfe_snps$gene
# unique_genes <- unique(nfe_snps$gene)
# dat <- rbind(dat, unique_genes)
# 
# 
# unique_genes <- unique(nfe_desc$gene)
# top <- unique_genes[1:gene_num]
# nfe_genes <- tibble(top)
# colnames(nfe_genes) <- "genes"
# dat <- rbind(dat, nfe_genes)
# 
# top_snps <- nfe_desc[nfe_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$nfe
# min_af.dat[14,1]<-"nfe"
# min_af.dat[14,2] <- min
# print(paste("Minimum allele frequency for nfe is ", min))
# nfe_snps <- nfe_desc[nfe_desc$nfe >= min,]
# write_delim(nfe_snps, file = paste0(base_path,nfe_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### FIN ###
# #######################################################################
# fin_desc <- fin_desc[!(is.na(fin_desc$fin)),]
# unique_genes <- unique(fin_desc$gene)
# top <- unique_genes[1:gene_num]
# fin_genes <- tibble(top)
# colnames(fin_genes) <- "genes"
# dat <- rbind(dat, fin_genes)
# 
# top_snps <- fin_desc[fin_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$fin
# min_af.dat[15,1]<-"fin"
# min_af.dat[15,2] <- min
# print(paste("Minimum allele frequency for fin is ", min))
# fin_snps <- fin_desc[fin_desc$fin >= min,]
# #write snps to file 
# write_delim(fin_snps, file = paste0(base_path, fin_out), 
#             delim = '\t', col_names = TRUE)
# 
# ########################################################################
# ### ASJ ###
# #######################################################################
# asj_desc <- asj_desc[!(is.na(asj_desc$asj)),]
# unique_genes <- unique(asj_desc$gene)
# top <- unique_genes[1:gene_num]
# asj_genes <- tibble(top)
# colnames(asj_genes) <- "genes"
# dat <- rbind(dat, asj_genes)
# 
# top_snps <- asj_desc[asj_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$asj
# min_af.dat[16,1]<-"asj"
# min_af.dat[16,2] <- min
# print(paste("Minimum allele frequency for asj is ", min))
# asj_snps <- asj_desc[asj_desc$asj >= min,]
# #write snps to file 
# write_delim(asj_snps, file = paste0(base_path,asj_out), 
#             delim = '\t', col_names = TRUE)
# 
# 
# ########################################################################
# ### OTH ###
# #######################################################################
# oth_desc <- oth_desc[!(is.na(oth_desc$oth)),]
# unique_genes <- unique(oth_desc$gene)
# top <- unique_genes[1:gene_num]
# oth_genes <- tibble(top)
# colnames(oth_genes) <- "genes"
# dat <- rbind(dat, oth_genes)
# 
# top_snps <- oth_desc[oth_desc$gene %in% top,]
# t.first <- top_snps[match(unique(top_snps$gene), top_snps$gene),]
# last_row <- nrow(t.first)
# min <- t.first[last_row,]$oth
# min_af.dat[17,1]<-"oth"
# min_af.dat[17,2] <- min
# print(paste("Minimum allele frequency for oth is ", min))
# oth_snps <- oth_desc[oth_desc$oth >= min,]
# #write snps to file 
# write_delim(oth_snps, file = paste0(base_path,oth_out), 
#             delim = '\t', col_names = TRUE)
# 
# 
# ########################################################################
# ### WRITE OUT MINIMUM ALLELE FREQUENCIES AND GENE LIST ###
# #######################################################################
# write_delim(dat, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/all_ancestry_top_af_genes.tsv", 
#             delim = '\t', col_names = TRUE)
# write_delim(min_af.dat, file = "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/minimum_allele_frequencies.tsv", 
#             delim = '\t', col_names = TRUE)

# network.sig <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/disease_networks/mammary_signature_mapped_network_list_threshold_2.tsv", col_names = FALSE)
# afr_desc <- exomeAF[exomeAF$afr > 0,]
# afr_tab <- table(afr_desc$afr)
# afr.tab.s <- afr_tab[order(afr_tab, decreasing = TRUE)]
# afr <- afr_desc[afr_desc$afr >= 0.0000615157 & afr_desc$afr <= 0.0000615233,]
# 
# 
# net.cut <- network[network$Connection > 0.8,]
# tt <- c(net.cut$Gene1, net.cut$Gene2)
# afr_sig$Variable %in% tt
# net.list <- c(net.cut$Gene1, net.cut$Gene2)
# afr_sig$Variable %in% afr$gene
# 
# afr <- afr[afr$gene %in% net.list,]
# length(unique(afr$gene))
# ge <- c(unique(afr$gene), consensus_genes$`Gene Symbol`)
