#### cosmic analysis of genes = what does phyloFrame consider to be EAF for each ancestry? #### 
library(dplyr) # %>% 
exomeAF <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/mean_enhancedAF_exome.tsv", col_names = TRUE)
ccgenes <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/Cancer_census_genes.tsv")
cc.exome <- exomeAF[exomeAF$gene %in% ccgenes$`Gene Symbol`,]
dim(cc.exome)
####  correct #### 
dec <- TRUE
afr.exome <- cc.exome[order(cc.exome$afr,decreasing = dec),]
afr.exome <- afr.exome[afr.exome$afr > 0.001 ,]
afr.mean <- mean(afr.exome$afr)
cosmic.afr <-  afr.exome[match(unique(afr.exome$gene), afr.exome$gene),]
dim(cosmic.afr)
cosmic.afr$ancestry <- "afr"
cosmic.afr <- cosmic.afr %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, afr, ancestry)
colnames(cosmic.afr) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")

all.cosmic <- cosmic.afr

nfe_seu.exome <- cc.exome[order(cc.exome$nfe_seu,decreasing = dec),]
nfe_seu.exome <- nfe_seu.exome[nfe_seu.exome$nfe_seu > 0.001 ,]
nfe_seu.mean <- mean(nfe_seu.exome$nfe_seu)
cosmic.nfe_seu <-  nfe_seu.exome[match(unique(nfe_seu.exome$gene), nfe_seu.exome$gene),]
dim(cosmic.nfe_seu)
cosmic.nfe_seu$ancestry <- "nfe_seu"
cosmic.nfe_seu <- cosmic.nfe_seu %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, nfe_seu, ancestry)
colnames(cosmic.nfe_seu) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.nfe_seu)

nfe_bgr.exome <- cc.exome[order(cc.exome$nfe_bgr,decreasing = dec),]
nfe_bgr.exome <- nfe_bgr.exome[nfe_bgr.exome$nfe_bgr > 0.001 ,]
nfe_bgr.mean <- mean(nfe_bgr.exome$nfe_bgr)
cosmic.nfe_bgr <-  nfe_bgr.exome[match(unique(nfe_bgr.exome$gene), nfe_bgr.exome$gene),]
dim(cosmic.nfe_bgr)
cosmic.nfe_bgr$ancestry <- "nfe_bgr"
cosmic.nfe_bgr <- cosmic.nfe_bgr %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, nfe_bgr, ancestry)
colnames(cosmic.nfe_bgr) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.nfe_bgr)

sas.exome <- cc.exome[order(cc.exome$sas,decreasing = dec),]
sas.exome <- sas.exome[sas.exome$sas > 0.001 ,]
sas.mean <- mean(sas.exome$sas)
cosmic.sas <-  sas.exome[match(unique(sas.exome$gene), sas.exome$gene),]
dim(cosmic.sas)
cosmic.sas$ancestry <- "sas"
cosmic.sas <- cosmic.sas %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, sas, ancestry)
colnames(cosmic.sas) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.sas)

amr.exome <- cc.exome[order(cc.exome$amr,decreasing = dec),]
amr.exome <- amr.exome[amr.exome$amr > 0.001 ,]
amr.mean <- mean(amr.exome$amr)
cosmic.amr <-  amr.exome[match(unique(amr.exome$gene), amr.exome$gene),]
dim(cosmic.amr)
cosmic.amr$ancestry <- "amr"
cosmic.amr <- cosmic.amr %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, amr, ancestry)
colnames(cosmic.amr) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.amr)

nfe_swe.exome <- cc.exome[order(cc.exome$nfe_swe,decreasing = dec),]
nfe_swe.exome <- nfe_swe.exome[nfe_swe.exome$nfe_swe > 0.001 ,]
nfe_swe.mean <- mean(nfe_swe.exome$nfe_swe)
cosmic.nfe_swe <-  nfe_swe.exome[match(unique(nfe_swe.exome$gene), nfe_swe.exome$gene),]
dim(cosmic.nfe_swe)
cosmic.nfe_swe$ancestry <- "nfe_swe"
cosmic.nfe_swe <- cosmic.nfe_swe %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, nfe_swe, ancestry)
colnames(cosmic.nfe_swe) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.nfe_swe)

nfe_nwe.exome <- cc.exome[order(cc.exome$nfe_nwe,decreasing = dec),]
nfe_nwe.exome <- nfe_nwe.exome[nfe_nwe.exome$nfe_nwe > 0.001 ,]
nfe_nwe.mean <- mean(nfe_nwe.exome$nfe_nwe)
cosmic.nfe_nwe <-  nfe_nwe.exome[match(unique(nfe_nwe.exome$gene), nfe_nwe.exome$gene),]
dim(cosmic.nfe_nwe)
cosmic.nfe_nwe$ancestry <- "nfe_nwe"
cosmic.nfe_nwe <- cosmic.nfe_nwe %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, nfe_nwe, ancestry)
colnames(cosmic.nfe_nwe) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.nfe_nwe)

eas_jpn.exome <- cc.exome[order(cc.exome$eas_jpn,decreasing = dec),]
eas_jpn.exome <- eas_jpn.exome[eas_jpn.exome$eas_jpn > 0.001 ,]
eas_jpn.mean <- mean(eas_jpn.exome$eas_jpn)
cosmic.eas_jpn <-  eas_jpn.exome[match(unique(eas_jpn.exome$gene), eas_jpn.exome$gene),]
dim(cosmic.eas_jpn)
cosmic.eas_jpn$ancestry <- "eas_jpn"
cosmic.eas_jpn <- cosmic.eas_jpn %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, eas_jpn, ancestry)
colnames(cosmic.eas_jpn) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.eas_jpn)

eas_kor.exome <- cc.exome[order(cc.exome$eas_kor,decreasing = dec),]
eas_kor.exome <- eas_kor.exome[eas_kor.exome$eas_kor > 0.001 ,]
eas_kor.mean <- mean(eas_kor.exome$eas_kor)
cosmic.eas_kor <-  eas_kor.exome[match(unique(eas_kor.exome$gene), eas_kor.exome$gene),]
dim(cosmic.eas_kor)
cosmic.eas_kor$ancestry <- "eas_kor"
cosmic.eas_kor <- cosmic.eas_kor %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, eas_kor, ancestry)
colnames(cosmic.eas_kor) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.eas_kor)

eas_oea.exome <- cc.exome[order(cc.exome$eas_oea,decreasing = dec),]
eas_oea.exome <- eas_oea.exome[eas_oea.exome$eas_oea > 0.001 ,]
eas_oea.mean <- mean(eas_oea.exome$eas_oea)
cosmic.eas_oea <-  eas_oea.exome[match(unique(eas_oea.exome$gene), eas_oea.exome$gene),]
dim(cosmic.eas_oea)
cosmic.eas_oea$ancestry <- "eas_oea"
cosmic.eas_oea <- cosmic.eas_oea %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, eas_oea, ancestry)
colnames(cosmic.eas_oea) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.eas_oea)

nfe_est.exome <- cc.exome[order(cc.exome$nfe_est,decreasing = dec),]
nfe_est.exome <- nfe_est.exome[nfe_est.exome$nfe_est > 0.001 ,]
nfe_est.mean <- mean(nfe_est.exome$nfe_est)
cosmic.nfe_est <-  nfe_est.exome[match(unique(nfe_est.exome$gene), nfe_est.exome$gene),]
dim(cosmic.nfe_est)
cosmic.nfe_est$ancestry <- "nfe_est"
cosmic.nfe_est <- cosmic.nfe_est %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, nfe_est, ancestry)
colnames(cosmic.nfe_est) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.nfe_est)

nfe.exome <- cc.exome[order(cc.exome$nfe,decreasing = dec),]
nfe.exome <- nfe.exome[nfe.exome$nfe > 0.001 ,]
nfe.mean <- mean(nfe.exome$nfe)
cosmic.nfe <-  nfe.exome[match(unique(nfe.exome$gene), nfe.exome$gene),]
dim(cosmic.nfe)
cosmic.nfe$ancestry <- "nfe"
cosmic.nfe <- cosmic.nfe %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, nfe, ancestry)
colnames(cosmic.nfe) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.nfe)

fin.exome <- cc.exome[order(cc.exome$fin,decreasing = dec),]
fin.exome <- fin.exome[fin.exome$fin > 0.001 ,]
fin.mean <- mean(fin.exome$fin)
cosmic.fin <-  fin.exome[match(unique(fin.exome$gene), fin.exome$gene),]
dim(cosmic.fin)
cosmic.fin$ancestry <- "fin"
cosmic.fin <- cosmic.fin %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, fin, ancestry)
colnames(cosmic.fin) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.fin)

asj.exome <- cc.exome[order(cc.exome$asj,decreasing = dec),]
asj.exome <- asj.exome[asj.exome$asj > 0.001 ,]
asj.mean <- mean(asj.exome$asj)
cosmic.asj <-  asj.exome[match(unique(asj.exome$gene), asj.exome$gene),]
dim(cosmic.asj)
cosmic.asj$ancestry <- "asj"
cosmic.asj <- cosmic.asj %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, asj, ancestry)
colnames(cosmic.asj) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.asj)

oth.exome <- cc.exome[order(cc.exome$oth,decreasing = dec),]
oth.exome <- oth.exome[oth.exome$oth > 0.001 ,]
oth.mean <- mean(oth.exome$oth)
cosmic.oth <-  oth.exome[match(unique(oth.exome$gene), oth.exome$gene),]
dim(cosmic.oth)
cosmic.oth$ancestry <- "oth"
cosmic.oth <- cosmic.oth %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, oth, ancestry)
colnames(cosmic.oth) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.oth)

nfe_onf.exome <- cc.exome[order(cc.exome$nfe_onf,decreasing = dec),]
nfe_onf.exome <- nfe_onf.exome[nfe_onf.exome$nfe_onf > 0.001 ,]
nfe_onf.mean <- mean(nfe_onf.exome$nfe_onf)
cosmic.nfe_onf <-  nfe_onf.exome[match(unique(nfe_onf.exome$gene), nfe_onf.exome$gene),]
dim(cosmic.nfe_onf)
cosmic.nfe_onf$ancestry <- "nfe_onf"
cosmic.nfe_onf <- cosmic.nfe_onf %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, nfe_onf, ancestry)
colnames(cosmic.nfe_onf) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.nfe_onf)

eas.exome <- cc.exome[order(cc.exome$eas,decreasing = dec),]
eas.exome <- eas.exome[eas.exome$eas > 0.001 ,]
eas.mean <- mean(eas.exome$eas)
cosmic.eas <-  eas.exome[match(unique(eas.exome$gene), eas.exome$gene),]
dim(cosmic.eas)
cosmic.eas$ancestry <- "eas"
cosmic.eas <- cosmic.eas %>% dplyr::select(gene, rs_id, ref_allele, alt_allele, eas, ancestry)
colnames(cosmic.eas) <- c("gene","rs_id","ref_allele","alt_allele","EAF","ancestry")
all.cosmic <- rbind(all.cosmic, cosmic.eas)

library(tidyr)
t <- spread(all.cosmic, key = ancestry, value = EAF)
dim(t)
write.table(all.cosmic, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/all_COSMIC_all_ancestries_top_altallele_EAF.tsv",sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(t, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/all_COSMIC_all_ancestries_top_altallele_EAF_spread.tsv",sep = '\t', col.names = TRUE, row.names = FALSE)

mean.all <- c(afr.mean, amr.mean, asj.mean, eas.mean, eas_jpn.mean, eas_kor.mean, eas_oea.mean, fin.mean, nfe.mean, nfe_bgr.mean, nfe_est.mean, nfe_nwe.mean, 
              nfe_onf.mean, nfe_seu.mean, nfe_swe.mean, oth.mean, sas.mean )

new.all <- all.cosmic %>% dplyr::select(gene, EAF, ancestry)
t.new <- spread(new.all, key = ancestry, value = EAF)
head(t.new)
to.add <- c("mean",mean.all)

testing <- rbind(t.new,to.add)

final <- rbind(t.new, to.add)

testing$afr <- paste0(final$afr, " - ", testing$afr)
testing$amr <- paste0(final$amr, " - ", testing$amr)
testing$asj <- paste0(final$asj, " - ", testing$asj)
testing$eas <- paste0(final$eas, " - ", testing$eas)
testing$eas_jpn <- paste0(final$eas_jpn, " - ", testing$eas_jpn)
testing$eas_kor <- paste0(final$eas_kor, " - ", testing$eas_kor)
testing$eas_oea <- paste0(final$eas_oea, " - ", testing$eas_oea)
testing$fin <- paste0(final$fin, " - ", testing$fin)
testing$nfe <- paste0(final$nfe, " - ", testing$nfe)
testing$nfe_bgr <- paste0(final$nfe_bgr, " - ", testing$nfe_bgr)
testing$nfe_est <- paste0(final$nfe_est, " - ", testing$nfe_est)
testing$nfe_nwe <- paste0(final$nfe_nwe, " - ", testing$nfe_nwe)
testing$nfe_onf <- paste0(final$nfe_onf, " - ", testing$nfe_onf)
testing$nfe_seu <- paste0(final$nfe_seu, " - ", testing$nfe_seu)
testing$nfe_swe <- paste0(final$nfe_swe, " - ", testing$nfe_swe)
testing$oth <- paste0(final$oth, " - ", testing$oth)
testing$sas <- paste0(final$sas, " - ", testing$sas)

tail(testing)
dim(testing)
testing <- testing[-710,]
tail(testing)

dim(testing)

testing <- rbind(testing, to.add)
tail(testing)

write.table(testing, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/all_COSMIC_all_ancestries_RANGE_MEAN.tsv",sep = '\t', col.names = TRUE, row.names = FALSE)


tempo <- all.cosmic %>% dplyr::group_by(gene, ancestry) %>% dplyr::summarise(across(c(afr, amr), sum))


l <- reshape(all.cosmic, idvar = c("gene","ancestry"), timevar = "EAF", direction = "wide")
cast(all.cosmic, gene ~ ancestry)

library(reshape)
all.cosmic <- unique(all.cosmic)
write.table(all.cosmic, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/all_COSMIC_all_ancestries_top_EAF.tsv",sep = '\t', col.names = TRUE, row.names = FALSE)


temp.exome <-cc.exome[cc.exome]




