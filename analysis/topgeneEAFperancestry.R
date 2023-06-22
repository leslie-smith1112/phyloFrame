############ LOOK AT TOP N EAF GENES PER ANCESTRY ############

library(tidyr)
library(dplyr)

exomeAF <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/mean_enhancedAF_exome.tsv", col_names = TRUE)

dat <- data.frame(matrix(ncol=6, nrow=0))

#ancestry, gene, mutation,EAF
#### AFR #### 
temp.afr <- exomeAF[order(exomeAF$afr, decreasing = TRUE),]
temp.afr <- temp.afr[!(is.na(temp.afr$gene)),]
top50 <- temp.afr[match(unique(temp.afr$gene), temp.afr$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, afr, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "afr"

dat <- rbind(dat, to.add)
colnames(dat) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")

### NFE_SEU ###
temp.nfe_seu <- exomeAF[order(exomeAF$nfe_seu, decreasing = TRUE),]
temp.nfe_seu <- temp.nfe_seu[!(is.na(temp.nfe_seu$gene)),]
top50 <- temp.nfe_seu[match(unique(temp.nfe_seu$gene), temp.nfe_seu$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, nfe_seu, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "nfe_seu"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

#### NFE_BGR###
temp.nfe_bgr <- exomeAF[order(exomeAF$nfe_bgr, decreasing = TRUE),]
temp.nfe_bgr <- temp.nfe_bgr[!(is.na(temp.nfe_bgr$gene)),]
top50 <- temp.nfe_bgr[match(unique(temp.nfe_bgr$gene), temp.nfe_bgr$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, nfe_bgr, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "nfe_bgr"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)
dim(dat)

### SAS ###
temp.sas <- exomeAF[order(exomeAF$sas, decreasing = TRUE),]
temp.sas <- temp.sas[!(is.na(temp.sas$gene)),]
top50 <- temp.sas[match(unique(temp.sas$gene), temp.sas$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, sas, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "sas"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### NFE_ONF ###
temp.nfe_onf <- exomeAF[order(exomeAF$nfe_onf, decreasing = TRUE),]
temp.nfe_onf <- temp.nfe_onf[!(is.na(temp.nfe_onf$gene)),]
top50 <- temp.nfe_onf[match(unique(temp.nfe_onf$gene), temp.nfe_onf$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, nfe_onf, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "nfe_onf"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### AMR ### 
temp.amr <- exomeAF[order(exomeAF$amr, decreasing = TRUE),]
temp.amr <- temp.amr[!(is.na(temp.amr$gene)),]
top50 <- temp.amr[match(unique(temp.amr$gene), temp.amr$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, amr, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "amr"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)
 #### EAS #### 
temp.eas <- exomeAF[order(exomeAF$eas, decreasing = TRUE),]
temp.eas <- temp.eas[!(is.na(temp.eas$gene)),]
top50 <- temp.eas[match(unique(temp.eas$gene), temp.eas$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, eas, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "eas"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### NFE_SWE ### 
temp.nfe_swe <- exomeAF[order(exomeAF$nfe_swe, decreasing = TRUE),]
temp.nfe_swe <- temp.nfe_swe[!(is.na(temp.nfe_swe$gene)),]
top50 <- temp.nfe_swe[match(unique(temp.nfe_swe$gene), temp.nfe_swe$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, nfe_swe, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "nfe_swe"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

#### NFE_NWE #### 
temp.nfe_nwe <- exomeAF[order(exomeAF$nfe_nwe, decreasing = TRUE),]
temp.nfe_nwe <- temp.nfe_nwe[!(is.na(temp.nfe_nwe$gene)),]
top50 <- temp.nfe_nwe[match(unique(temp.nfe_nwe$gene), temp.nfe_nwe$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, nfe_nwe, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "nfe_nwe"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### EAS JPN ### 
temp.eas_jpn <- exomeAF[order(exomeAF$eas_jpn, decreasing = TRUE),]
temp.eas_jpn <- temp.eas_jpn[!(is.na(temp.eas_jpn$gene)),]
top50 <- temp.eas_jpn[match(unique(temp.eas_jpn$gene), temp.eas_jpn$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, eas_jpn, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "eas_jpn"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### EAS KOR #### 
temp.eas_kor <- exomeAF[order(exomeAF$eas_kor, decreasing = TRUE),]
temp.eas_kor <- temp.eas_kor[!(is.na(temp.eas_kor$gene)),]
top50 <- temp.eas_kor[match(unique(temp.eas_kor$gene), temp.eas_kor$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, eas_kor, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "eas_kor"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### EAS OEA ### --
temp.eas_oea <- exomeAF[order(exomeAF$eas_oea, decreasing = TRUE),]
temp.eas_oea <- temp.eas_oea[!(is.na(temp.eas_oea$gene)),]
top50 <- temp.eas_oea[match(unique(temp.eas_oea$gene), temp.eas_oea$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, eas_oea, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "eas_oea"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### NFE EAST ### 
temp.nfe_est <- exomeAF[order(exomeAF$nfe_est, decreasing = TRUE),]
temp.nfe_est <- temp.nfe_est[!(is.na(temp.nfe_est$gene)),]
top50 <- temp.nfe_est[match(unique(temp.nfe_est$gene), temp.nfe_est$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, nfe_est, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "nfe_est"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### NFE ### 
temp.nfe <- exomeAF[order(exomeAF$nfe, decreasing= TRUE),]
temp.nfe <- temp.nfe[!(is.na(temp.nfe$gene)),]
top50 <- temp.nfe[match(unique(temp.nfe$gene), temp.nfe$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, nfe, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "nfe"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### FIN ### 
temp.fin <- exomeAF[order(exomeAF$fin, decreasing = TRUE),]
temp.fin <- temp.fin[!(is.na(temp.fin$gene)),]
top50 <- temp.fin[match(unique(temp.fin$gene), temp.fin$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, fin, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "fin"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

### ASJ ### 
temp.asj <- exomeAF[order(exomeAF$asj, decreasing = TRUE),]
temp.asj <- temp.asj[!(is.na(temp.asj$gene)),]
top50 <- temp.asj[match(unique(temp.asj$gene), temp.asj$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, asj, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "asj"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)
### OTH ### 
temp.oth <- exomeAF[order(exomeAF$oth, decreasing = TRUE),]
temp.oth <- temp.oth[!(is.na(temp.oth$gene)),]
top50 <- temp.oth[match(unique(temp.oth$gene), temp.oth$gene),]
top50 <- top50[1:50,] 
top50$gene
to.add <- top50 %>% dplyr::select(gene, oth, rs_id, ref_allele, alt_allele)
to.add$ancestry <- "oth"
colnames(to.add) <- c("gene","EAF", "rs_id","ref_allele","alt_allele","ancestry")
dat <- rbind(dat, to.add)

write.table(dat, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/top50EAFgeneperancestry.tsv",sep = '\t', col.names = TRUE, row.names = FALSE)


#### get all foxa1 muttions #### 
foxa <- exomeAF[exomeAF$gene == "FOXA1",]
foxs <- foxa[!(is.na(foxa$gene)),]
wide <- gather(foxs, ancestry, EAF, nfe_seu:oth, factor_key=TRUE)
head(wide)
pos <- wide[wide$EAF > 0,]
dim(pos)
was <- pos[order(pos$EAF, decreasing = TRUE),]
write.table(was, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/foxa1EAFperancestry.tsv",sep = '\t', col.names = TRUE, row.names = FALSE)









