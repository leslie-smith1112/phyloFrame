# 
# COSMIC gene enhanced allele frequency (EAF) comparisons - e.g. t-test AFR vs EUR EAF - correlated? Correlated in genes in PF vs BM?
#   % questions to answer: this is the lollupop plot we discussed 
# % are the COSMIC genes EAF enriched in any populations? Which ones are most commonly enriched? Are there drugs (PharmGKB lookup) that target those genes, and if so are they specific to a population or do they work in many?
#   % of the genes that are COSMIC but also EAF high in some populations - look at the altered alleles in those populations & predict (ExPecTo model at HB website) the functional impact of that alteration on the gene. 
# Are those genes likely functionally altered in some populations but not others?
# 
#change per disease
dir <- "/breastph08TESTAF001_1_network02_051_neigh2_top30VarFinal"
disease <- "breast"
ccgenes <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/Cancer_census_genes.tsv")
exomeAF <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/mean_enhancedAF_exome.tsv", col_names = TRUE)
cc.exome <- exomeAF[exomeAF$gene %in% ccgenes$`Gene Symbol`,]

#eafs for ceratin genes 
somatic <- ccgenes[grepl("breast", ccgenes$`Tumour Types(Somatic)`, ignore.case = TRUE,fixed= TRUE),]
germlin <- ccgenes[grepl("breast", ccgenes$`Tumour Types(Germline)`, ignore.case = TRUE,fixed= TRUE),]
cc.breast <- rbind(somatic, germlin)

#afr has higher enhancment in cosmic genes than nfe and eas 
t.test(cc.exome$nfe, cc.exome$afr)
t.test(cc.exome$nfe, cc.exome$eas)
t.test(cc.exome$eas, cc.exome$afr)

afr.eeas <- c(cc.exome$eas,cc.exome$afr) 
t.test(afr.eeas, cc.exome$nfe)

# analysis of benchmark and phyloframe disease signatures 
pf.eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/eur/")
pf.afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/afr/")
pf.admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/admixed/")
pf.eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/eas/")
pf.mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/phyloFrame/mixed/")

eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/eur/")
afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/afr/")
admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/admixed/")
eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/eas/")
mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,dir,"/model_runs/mixed/")
#read in all signatures and bind them togather by ancestry 
eur.num <- 1:17
pf.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
pf.eur.myfiles <- lapply(pf.eur.names, readr::read_tsv)
afr.num <- 1:2
pf.afr.names <- paste0(pf.afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
pf.afr.myfiles <- lapply(pf.afr.names, readr::read_tsv)
eas.num <- 1
pf.eas.names <- paste0(pf.eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt") 
pf.eas.myfiles <- lapply(pf.eas.names, readr::read_tsv)
admixed.num <- 1
pf.admixed.names <- paste0(pf.admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt")
pf.admixed.myfiles <- lapply(pf.admixed.names, readr::read_tsv)
mixed.num <- 1:6
pf.mixed.names <- paste0(pf.mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt")
pf.mixed.myfiles <- lapply(pf.mixed.names, readr::read_tsv)

df.eur.names <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
eur.myfiles <- lapply(df.eur.names, readr::read_tsv)
df.afr.names <- paste0(afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
afr.myfiles <- lapply(df.afr.names, readr::read_tsv)
df.eas.names <- paste0(eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt") 
eas.myfiles <- lapply(df.eas.names, readr::read_tsv)
df.admixed.names <- paste0(admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt")
admixed.myfiles <- lapply(df.admixed.names, readr::read_tsv)
df.mixed.names <- paste0(mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt")
mixed.myfiles <- lapply(df.mixed.names, readr::read_tsv)

eur.sig <- do.call(rbind,pf.eur.myfiles)
afr.sig <- do.call(rbind,pf.afr.myfiles)
eas.sig <- do.call(rbind,pf.eas.myfiles)
admixed.sig <- do.call(rbind, pf.admixed.myfiles)
mixed.sig <- do.call(rbind, pf.mixed.myfiles)
all.pf.sig <- rbind(eur.sig, afr.sig, eas.sig, admixed.sig, mixed.sig)
#keep only genes in the cosmic 
pf.cc.genes <- cc.exome[ cc.exome$gene %in% all.pf.sig$Variable,]

bm.eur.sig <- do.call(rbind,eur.myfiles)
bm.afr.sig <- do.call(rbind,afr.myfiles)
bm.eas.sig <- do.call(rbind,eas.myfiles)
bm.admixed.sig <- do.call(rbind, admixed.myfiles)
bm.mixed.sig <- do.call(rbind, mixed.myfiles)
all.bm.sig <- rbind(bm.eur.sig, bm.afr.sig, bm.eas.sig, bm.admixed.sig, bm.mixed.sig)
#keep only genes in the cosmic 
bm.cc.genes <- cc.exome[ cc.exome$gene %in% all.bm.sig$Variable,]

#t.test pf vs benchmark cosmic signature genes 
t.test(pf.cc.genes$sas, bm.cc.genes$sas)
t.test(pf.cc.genes$amr, bm.cc.genes$amr)
t.test(pf.cc.genes$fin, bm.cc.genes$fin)
t.test(pf.cc.genes$nfe_seu, bm.cc.genes$nfe_seu)
t.test(pf.cc.genes$nfe_bgr, bm.cc.genes$nfe_bgr)
t.test(pf.cc.genes$afr, bm.cc.genes$afr)
t.test(pf.cc.genes$nfe_onf, bm.cc.genes$nfe_onf)
t.test(pf.cc.genes$eas, bm.cc.genes$eas)
t.test(pf.cc.genes$nfe_swe, bm.cc.genes$nfe_swe)
t.test(pf.cc.genes$nfe_nwe, bm.cc.genes$nfe_nwe)
t.test(pf.cc.genes$eas_jpn, bm.cc.genes$eas_jpn)
t.test(pf.cc.genes$eas_oea, bm.cc.genes$eas_oea)
t.test(pf.cc.genes$nfe_est, bm.cc.genes$nfe_est)
t.test(pf.cc.genes$nfe, bm.cc.genes$nfe)
t.test(pf.cc.genes$asj, bm.cc.genes$asj)
t.test(pf.cc.genes$oth, bm.cc.genes$oth)


## across all cosmic genes 
mean(cc.exome$nfe_seu)
mean(cc.exome$nfe_bgr)
mean(cc.exome$afr)
mean(cc.exome$sas)
mean(cc.exome$nfe_onf)
mean(cc.exome$amr)
mean(cc.exome$eas)
mean(cc.exome$nfe_swe)
mean(cc.exome$nfe_nwe)
mean(cc.exome$eas_jpn)
mean(cc.exome$eas_kor)
mean(cc.exome$eas_oea)
mean(cc.exome$nfe_est)
mean(cc.exome$nfe)
mean(cc.exome$fin)
mean(cc.exome$asj)


# get mean of cosmic genes in pf and benchmark signatures for cosmic genes
dat <- data.frame(matrix(ncol = 2, nrow = 0))
##EUR
for(i in 1:length(pf.eur.myfiles)){
  current <- pf.eur.myfiles[i][[1]]$Variable
  common <- length(intersect(current, ccgenes$`Gene Symbol`))
  current.bm <- eur.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, ccgenes$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

##AFR
for(i in 1:length(pf.afr.myfiles)){
  current <- pf.afr.myfiles[i][[1]]$Variable
  common <- length(intersect(current, ccgenes$`Gene Symbol`))
  current.bm <- afr.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, ccgenes$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

##EAS
for(i in 1:length(pf.eas.myfiles)){
  current <- pf.eas.myfiles[i][[1]]$Variable
  common <- length(intersect(current, ccgenes$`Gene Symbol`))
  current.bm <- eas.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, ccgenes$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}


##ADMIXED
for(i in 1:length(pf.admixed.myfiles)){
  current <- pf.admixed.myfiles[i][[1]]$Variable
  common <- length(intersect(current, ccgenes$`Gene Symbol`))
  current.bm <- admixed.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, ccgenes$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

##MIXED
for(i in 1:length(pf.mixed.myfiles)){
  current <- pf.mixed.myfiles[i][[1]]$Variable
  common <- length(intersect(current, ccgenes$`Gene Symbol`))
  current.bm <- mixed.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, ccgenes$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

t.test(dat$phyloFrame, dat$benchmark)
mean(dat$phyloFrame)
mean(dat$benchmark)


## breast specific cosmic genes 
dat <- data.frame(matrix(ncol = 2, nrow = 0))
##EUR
for(i in 1:length(pf.eur.myfiles)){
  current <- pf.eur.myfiles[i][[1]]$Variable
  common <- length(intersect(current, cc.breast$`Gene Symbol`))
  current.bm <- eur.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, cc.breast$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

##AFR
for(i in 1:length(pf.afr.myfiles)){
  current <- pf.afr.myfiles[i][[1]]$Variable
  common <- length(intersect(current, cc.breast$`Gene Symbol`))
  current.bm <- afr.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, cc.breast$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

##EAS
for(i in 1:length(pf.eas.myfiles)){
  current <- pf.eas.myfiles[i][[1]]$Variable
  common <- length(intersect(current, cc.breast$`Gene Symbol`))
  current.bm <- eas.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, cc.breast$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}


##ADMIXED
for(i in 1:length(pf.admixed.myfiles)){
  current <- pf.admixed.myfiles[i][[1]]$Variable
  common <- length(intersect(current, cc.breast$`Gene Symbol`))
  current.bm <- admixed.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, cc.breast$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

##MIXED
for(i in 1:length(pf.mixed.myfiles)){
  current <- pf.mixed.myfiles[i][[1]]$Variable
  common <- length(intersect(current, cc.breast$`Gene Symbol`))
  current.bm <- mixed.myfiles[i][[1]]$Variable
  common.bm <- length(intersect(current.bm, cc.breast$`Gene Symbol`))
  to.add <- c(common, common.bm)
  dat <- rbind(dat,to.add)
  colnames(dat)<- c("phyloFrame","benchmark")
}

t.test(dat$phyloFrame, dat$benchmark)
mean(dat$phyloFrame)
mean(dat$benchmark)
dim(cc.breast) 

## overall mutation burden of enhanced mutations
#afr: 2495346 mutations, 31287 genes 
#eur: 6387300 mutations, 31775 genes 
#eas: 2201382 mutations, 31084 genes 
afr <- exomeAF[exomeAF$afr > 0,]
nfe <- exomeAF[exomeAF$nfe >0,]
eas <- exomeAF[exomeAF$eas > 0,]
t.test(nfe$nfe,afr$afr) # mean x: 0.0005892796 mean y: 0.0060728146; 95% confidence: 0.005443883 0.005523187;t = 271.05, df = 2632226, p-value < 2.2e-16
####NOTE NOT SURE WHY THE CONFIDENCE INTERVAL IS NEGATIVE BELOW### 
t.test(nfe$nfe, eas$eas) # mean x: 0.0005892796 mean y: 0.0037775465; 95% confidence: -0.003222482 -0.003154051;t = -182.63, df = 2365671, p-value < 2.2e-16
t.test(afr$afr, eas$eas) # mean x: 0.006072815 mean y:0.003777546; 95% confidence: 0.002243696 0.002346840; t = 87.231, df = 4660175, p-value < 2.2e-16

afr.h <- exomeAF[exomeAF$afr > 0.5,]#803 mutations 587 genes 
nfe.h <- exomeAF[exomeAF$nfe >0.5,] # 254 mutations in 179 genes 
eas.h <- exomeAF[exomeAF$eas > 0.5,] # 277 mutations in 198 genes 


t.test(nfe$nfe.,afr.h$afr.h ) # mean x: 0.0005892796 mean y: 0.0060728146; 95% confidence: 0.005443883 0.005523187;t = 271.05, df = 2632226, p-value < 2.2e-16
####NOTE NOT SURE WHY THE CONFIDENCE INTERVAL IS NEGATIVE BELOW### 
t.test(nfe$nfe, eas$eas) # mean x: 0.0005892796 mean y: 0.0037775465; 95% confidence: -0.003222482 -0.003154051;t = -182.63, df = 2365671, p-value < 2.2e-16
t.test(afr$afr, eas$eas) 
t.test(exomeAF$nfe, exomeAF$afr)

cc.pos.nfe<- cc.exome[cc.exome$nfe > 0,]# 326851     710 genes; range: 9.125000e-10 7.708333e-01 <- mutation burdens of positibely enhanced mutations 
cc.pos.afr<- cc.exome[cc.exome$afr > 0,]# 125466     707 genes; range: 1.688750e-08 8.146301e-01
cc.pos.eas<- cc.exome[cc.exome$eas > 0,]# 110490     709 genes; range: 2.812500e-09 6.627535e-01

cc.pos.nfe<- cc.exome[cc.exome$nfe > 0.5,]# 6 mutaions ; 5 genes  
cc.pos.afr<- cc.exome[cc.exome$afr > 0.5,]# 43 mutations, 30 genes 
cc.pos.eas<- cc.exome[cc.exome$eas > 0.5,] #4 mutations, 3 genes 
t.test(cc.pos.nfe$nfe, cc.pos.afr$afr)# mean of x: 0.6420603, mean of y: 0.5835454, 95% confidence -0.03721279  0.15424267, t = 1.4803, df = 6.2684, p-value = 0.1872
################################################# COSMIC gene enhanced allele frequency (EAF) comparisons - e.g. t-test AFR vs EUR EAF - correlated? ##################
## with all mutatiins in cosmic genes ## 
t.test(cc.exome$nfe, cc.exome$afr) ### mean of x: -2.870513e-05 y: 3.390320e-04; t = -22.561, df = 1087484, p-value < 2.2e-16
t.test(cc.exome$nfe, cc.exome$eas) ### mean of x: -2.870513e-05  y: 3.811005e-05; t = -5.4071, df = 1275172, p-value = 6.407e-08; 95 % confidence:-9.103438e-05 -4.259598e-05
t.test(cc.exome$eas, cc.exome$afr) ### mean of x: 3.811005e-05 y: 3.390320e-04; t = -15.931, df = 1571751, p-value < 2.2e-16; 95 % confidence:-0.000337944 -0.000263900
## with only positive EAF for each ancestry ## 
t.test(cc.pos.nfe$nfe, cc.pos.afr$afr) ### mean of x: 0.0004514993 y: 0.0053096807; t = -57.846, df = 129927, p-value < 2.2e-16
t.test(cc.pos.nfe$nfe, cc.pos.eas$eas) ### mean of x: 0.0004514993 y: 0.0031014318 ; t = -39.807, df = 116834, p-value < 2.2e-16; 95 % confidence:-0.002780409 -0.002519456
t.test(cc.pos.afr$afr, cc.pos.eas$eas) ### mean of x: 0.005309681 y: 0.003101432 ; t = 20.829, df = 229316, p-value < 2.2e-16; 95 % confidence:0.002000457 0.002416041

#TAKE AWAY: in cosmic genes east asian and afr have higher avg EAF than nfe, even while nfe have more than double the enhanced mutations.
################################################# COSMIC gene enhanced allele frequency (EAF) comparisons - Correlated in genes in PF vs BM? ##################
## read in all PF and BM signatures ## 
## PF ## 
pf.genes <- exomeAF[exomeAF$gene %in% all.pf.sig$Variable,]
bm.genes <- exomeAF[exomeAF$gene %in% all.bm.sig$Variable,]
bm.genes <- bm.genes[bm.genes$gene %in% ccgenes$`Gene Symbol`,]
pf.genes <- pf.genes[pf.genes$gene %in% ccgenes$`Gene Symbol`,]
t.test(pf.genes$sas, bm.genes$sas)

##PF ## 
t.test(pf.cc.genes$nfe, pf.cc.genes$afr) ###mean of x: 6.770325e-06 y: 4.437332e-04; t = -7.7439, df = 93158, p-value = 9.741e-15; 95 percent confidence interval: -0.0005475592 -0.0003263665
t.test(pf.cc.genes$nfe, pf.cc.genes$eas) ###mean of x: 6.770325e-06 y: 2.198268e-05; t = -1.3989, df = 109220, p-value = 0.1619; 95 percent confidence interval: -1.315891e-04  2.198268e-05
t.test(pf.cc.genes$afr, pf.cc.genes$eas) ###mean of x: 4.437332e-04 y: 6.157353e-05; t = 5.8765, df = 138765, p-value = 4.2e-09; 95 percent confidence interval: 0.0002546983 0.0005096211
## BM ## 
t.test(bm.cc.genes$nfe, bm.cc.genes$afr) ###mean of x: -6.887963e-06 y: 2.791282e-04; t = -8.3532, df = 250112, p-value < 2.2e-16; 95 percent confidence interval: -0.0003531263 -0.0002189059
t.test(bm.cc.genes$nfe, bm.cc.genes$eas) ###mean of x: -6.887963e-06  y: -1.018985e-05; t = 0.12333, df = 291827, p-value = 0.9018; 95 percent confidence interval: -4.917200e-05  5.577578e-05
t.test(bm.cc.genes$afr, bm.cc.genes$eas) ###mean of x: 2.791282e-04  y: -1.018985e-05; t = 7.2863, df = 355008, p-value = 3.193e-13; 95 percent confidence interval: 0.000211493 0.000367143
t.test(pf.cc.genes$sas, bm.cc.genes$sas)

## small check for cosmic gene specifically in a signature ## 
eur.pf.cc.genes <- cc.exome[ cc.exome$gene %in% eur.sig$Variable,]
eur.bm.cc.genes <- cc.exome[ cc.exome$gene %in% bm.eur.sig$Variable,]

t.test(eur.pf.cc.genes$nfe, eur.pf.cc.genes$afr) ###mean of x: 1.534737e-05 y: 1.629969e-04; t = -2.7618, df = 77718, p-value = 0.005749; 95 percent confidence interval: -2.524323e-04 -4.286668e-05
t.test(eur.bm.cc.genes$nfe, eur.bm.cc.genes$afr) ###mean of x: 8.860610e-07 y: 2.769743e-04 ; t = -6.2115, df = 151907, p-value = 5.26e-10; 95 percent confidence interval: -0.0003632046 -0.0001889719

### nott ignificant ## 
t.test(pf.cc.genes$eas, bm.cc.genes$eas)
t.test(pf.cc.genes$nfe, bm.cc.genes$nfe)
t.test(pf.cc.genes$afr, bm.cc.genes$afr)
t.test(pf.cc.genes$sas, bm.cc.genes$sas)#
t.test(pf.cc.genes$amr, bm.cc.genes$amr)#
t.test(pf.cc.genes$nfe_seu, bm.cc.genes$nfe_seu)
t.test(pf.cc.genes$nfe_bgr, bm.cc.genes$nfe_bgr)
t.test(pf.cc.genes$nfe_onf, bm.cc.genes$nfe_onf)
t.test(pf.cc.genes$eas_jpn, bm.cc.genes$eas_jpn)
t.test(pf.cc.genes$fin, bm.cc.genes$fin)#
t.test(pf.cc.genes$asj, bm.cc.genes$asj)
t.test(pf.cc.genes$nfe_est, bm.cc.genes$nfe_est)
t.test(pf.cc.genes$nfe_swe, bm.cc.genes$nfe_swe)#
t.test(pf.cc.genes$nfe_nwe, bm.cc.genes$nfe_nwe)#
t.test(pf.cc.genes$eas_kor, bm.cc.genes$eas_kor)#
t.test(pf.cc.genes$eas_oea, bm.cc.genes$eas_oea)
t.test(pf.cc.genes$oth, bm.cc.genes$oth)

all.bm.sig
t.test(pf.cc.genes$sas, bm.cc.genes$sas)
t.test(pf.cc.genes$eas_kor, bm.cc.genes$eas_kor)
t.test(pf.cc.genes$amr, bm.cc.genes$amr)

t.test(pf.cc.genes$asj, bm.cc.genes$asj)

################################################# are the COSMIC genes EAF enriched in any populations? ############################################## 
## how to define enriched - mutatoin murden or most enhanced mutations? ## 
eas.exome <- cc.exome[order(cc.exome$eas, decreasing = TRUE),]
eas.exome$gene[1:100]
eas.exome <- eas.exome[eas.exome$eas > 0 ,]
mean(eas.exome$eas) ## avg EAF: 0.003101432
range(eas.exome$eas) ## range: 2.812500e-09 6.627535e-01
eas.genes <- unique(eas.exome$gene)
length(eas.genes) # 709 gnees with at least 1 enhanced mutation 
eas.mut.burden <- table(eas.exome$gene)
eas.mut.burden <- eas.mut.burden[order(eas.mut.burden,decreasing = TRUE)]
eas.mut <- eas.mut.burden[1:100]

c(eas.genes, eas.mut.burden)
### top mutated genes in east asian ###
# MUC16    MUC4  RNF213   LRP1B   KMT2D   BIRC6   KMT2C    FAT1 PDE4DIP   NCOR2   TRRAP   CSMD3    FAT3    FAT4   FANCA 
# 2112    2070     821     764     739     716     716     672     671     647     625     597     576     576     574 
# POLE   AKAP9    TSC2   PTPRD     ATM 
# 574     529     505     484     477 

### top 50 most enhanced genes in east asian ### 
# [1] "RHOA"     "PTK6"     "NUTM2D"   "MGMT"     "PTPN11"   "FANCA"    "CRNKL1"   "NOTCH1"   "RGPD3"    "ETNK1"   
# [11] "SH2B3"    "MACC1"    "SFPQ"     "WT1"      "ARHGAP26" "TEC"      "BCOR"     "SIRPA"    "ALK"      "RUNX1"   
# [21] "RABEP1"   "EPAS1"    "STK11"    "NCKIPSD"  "SEPT9"    "PLCG1"    "CDC73"    "FLT4"     "GAS7"     "FHIT"    
# [31] "TCF3"     "LATS2"    "GNAS"     "PIK3R1"   "NF2"      "HLF"      "CUX1"     "BCR"      "TSC1"     "NRG1"    
# [41] "ASXL1"    "MYO5A"    "HIST1H3B" "RAC1"     "PIK3CB"   "STRN"     "LATS1"    "FES"      "LPP"      "ETV5"

afr.exome <- cc.exome[order(cc.exome$afr, decreasing = TRUE),]
afr.exome$gene[1:50]
afr.exome <- afr.exome[afr.exome$afr > 0 ,]
mean(afr.exome$afr) ## avg EAF: 0.005309681
range(afr.exome$afr) ## range: 1.688750e-08 8.146301e-01
afr.genes <- unique(afr.exome$gene)
length(afr.genes) # 707 gnees with at least 1 enhanced mutation 
afr.mut.burden <- table(afr.exome$gene)
afr.mut.burden <- afr.mut.burden[order(afr.mut.burden,decreasing = TRUE)]
afr.mut.burden[1:20]
afr.mut <- afr.mut.burden[1:100] 
afr.genes
### top mutated genes in afr ###
# MUC4   MUC16  RNF213   KMT2C    FAT1   LRP1B   KMT2D PDE4DIP   CSMD3   NCOR2   BIRC6    POLE    FAT3   FANCA   TRRAP 
# 3493    2641     946     864     855     853     826     793     703     702     688     688     679     676     659 
# FAT4   AKAP9    PCM1  NOTCH1    TSC2 
# 628     604     587     570     569 

### top 50 most enhanced genes in afr ### 
# [1] "AMER1"    "RUNX1"    "MSN"      "PTK6"     "MSN"      "FAM47C"   "FAM47C"   "PTK6"     "DDX5"     "FLNA"    
# [11] "DDX5"     "DDX5"     "DDX5"     "MLLT3"    "LPP"      "FANCD2"   "TP53"     "LMNA"     "RPL22"    "KDM6A"   
# [21] "KDM6A"    "BTK"      "ERC1"     "ABL1"     "KIAA1549" "TCF7L2"   "CDKN1B"   "TRRAP"    "RANBP2"   "BCOR"    
# [31] "TCEA1"    "AR"       "PRPF40B"  "RANBP2"   "FAT3"     "BRAF"     "RANBP2"   "BRAF"     "ATRX"     "LMNA"    
# [41] "TSC2"     "ATRX"     "TRRAP"    "ABL1"     "MGMT"     "CREB3L1"  "SIX1"     "SEPT9"    "TSC2"     "MUC4"
# 

eur.exome <- cc.exome[order(cc.exome$nfe, decreasing = TRUE),]
eur.exome$gene[1:50]
eur.exome <- eur.exome[eur.exome$nfe > 0 ,]
mean(eur.exome$nfe) ## avg EAF: 0.0004514993
range(eur.exome$nfe) ## range: 9.125000e-10 7.708333e-01
eur.genes <- unique(eur.exome$gene)
length(eur.genes) # 710 gnees with at least 1 enhanced mutation 
eur.mut.burden <- table(eur.exome$gene)
eur.mut.burden <- eur.mut.burden[order(eur.mut.burden,decreasing = TRUE)]
eur.mut.burden[1:20]

### top mutated genes in eur ###
# MUC16    MUC4   LRP1B   KMT2D  RNF213   KMT2C    FAT1   BIRC6   CSMD3    FAT4 PDE4DIP    FAT3    POLE   AKAP9   NCOR2 
# 6935    4220    2368    2330    2317    2086    1983    1976    1952    1878    1850    1791    1788    1758    1703 
# FANCA   TRRAP     ATM    PCM1   BRCA2 
# 1694    1664    1548    1393    1379 

### top 50 most enhanced genes in eur ### 
# [1] "LASP1"    "MYCL"     "ABL1"     "STRN"     "MAP3K1"   "ABL1"     "MGMT"     "SEPT9"    "FHIT"     "DCC"     
# [11] "MGMT"     "MUC1"     "MUC1"     "GAS7"     "GNA11"    "PER1"     "RHOA"     "CXCR4"    "BCL9"     "KIF5B"   
# [21] "DNAJB1"   "NTRK1"    "TAL1"     "ZNF429"   "LRP1B"    "ASPSCR1"  "LEPROTL1" "FLT3"     "PER1"     "PER1"    
# [31] "CBLB"     "USP6"     "WRN"      "PCM1"     "USP6"     "FLNA"     "EZR"      "GAS7"     "PCM1"     "DROSHA"  
# [41] "TFRC"     "MGMT"     "AKT1"     "MAP3K1"   "PML"      "NCOA2"    "FHIT"     "MYO5A"    "MAP3K1"   "FANCG" 

## mutated 
mutated <- c() #hold genes most mutation burde 
highest <- c()#holds genes with highest EAF 


afr <- cc.exome[order(cc.exome$afr, decreasing = TRUE),]
afr<-afr[1:100,]
nfe_seu <- cc.exome[order(cc.exome$nfe_seu, decreasing = TRUE),]
nfe_seu <- nfe_seu[1:100,]
nfe_bgr <- cc.exome[order(cc.exome$nfe_bgr, decreasing = TRUE),]
nfe_bgr <- nfe_bgr[1:100,]
sas <- cc.exome[order(cc.exome$sas, decreasing = TRUE),]
sas <- sas[1:100,]
nfe_onf <- cc.exome[order(cc.exome$nfe_onf, decreasing = TRUE),]
nfe_onf<-nfe_onf[1:100,]
amr <- cc.exome[order(cc.exome$amr, decreasing = TRUE),]
amr <- amr[1:100,]
eas <- cc.exome[order(cc.exome$eas, decreasing = TRUE),]
eas <- eas[1:100,]
nfe_swe <- cc.exome[order(cc.exome$nfe_swe, decreasing = TRUE),]
nfe_swe <- nfe_swe[1:100,]
nfe_nwe <- cc.exome[order(cc.exome$nfe_nwe, decreasing = TRUE),]
nfe_nwe <- nfe_nwe[1:100,]
eas_jpn <- cc.exome[order(cc.exome$eas_jpn, decreasing = TRUE),]
eas_jpn <- eas_jpn[1:100,]
eas_kor <- cc.exome[order(cc.exome$eas_kor, decreasing = TRUE),]
eas_kor <- eas_kor[1:100,]
eas_oea <- cc.exome[order(cc.exome$eas_oea, decreasing = TRUE),]
eas_oea <- eas_oea[1:100,]
nfe_est <- cc.exome[order(cc.exome$nfe_est, decreasing = TRUE),]
nfe_est <- nfe_est[1:100,]
nfe <- cc.exome[order(cc.exome$nfe, decreasing = TRUE),]
nfe <- nfe[1:100,]
fin <- cc.exome[order(cc.exome$fin, decreasing = TRUE),]
fin <- fin[1:100,]
asj<- cc.exome[order(cc.exome$asj, decreasing = TRUE),]
asj <- asj[1:100,]

head(afr$gene)

afr1 <- c(nfe_seu$gene, nfe_bgr$gene, sas$gene, nfe_onf$gene, amr$gene, eas$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene, eas_kor$gene,
          eas_oea$gene, nfe_est$gene, nfe$gene, fin$gene, asj$gene)

nfe_seu1 <- c(afr$gene, sas$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene, eas$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene, eas_kor$gene,
              eas_oea$gene, nfe_est$gene,  fin$gene, asj$gene)

sas1 <- c(afr$gene, nfe_bgr$gene, nfe_seu$gene, nfe_onf$gene, amr$gene, eas$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene, eas_kor$gene,
          eas_oea$gene, nfe_est$gene, nfe$gene, fin$gene, asj$gene)

nfe_bgr1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_onf$gene, amr$gene, eas$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene, eas_kor$gene,
              eas_oea$gene, nfe_est$gene,  fin$gene, asj$gene)

nfe_onf1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, amr$gene, eas$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene, eas_kor$gene,
              eas_oea$gene, nfe_est$gene, fin$gene, asj$gene)

amr1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, eas$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene, eas_kor$gene,
          eas_oea$gene, nfe_est$gene, nfe$gene, fin$gene, asj$gene)

eas1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene, nfe_swe$gene, nfe_nwe$gene, 
          nfe_est$gene, nfe$gene, fin$gene, asj$gene)

nfe_swe1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene, eas$gene, nfe_nwe$gene, eas_jpn$gene, eas_kor$gene,
              eas_oea$gene, nfe_est$gene,  fin$gene, asj$gene)

nfe_nwe1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene, eas$gene, nfe_swe$gene, eas_jpn$gene, eas_kor$gene,
              eas_oea$gene, nfe_est$gene, fin$gene, asj$gene)

eas_jpn1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene,  nfe_swe$gene, nfe_nwe$gene, eas_kor$gene,
              eas_oea$gene, nfe_est$gene, nfe$gene, fin$gene, asj$gene)

eas_kor1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene,
              eas_oea$gene, nfe_est$gene, nfe$gene, fin$gene, asj$gene)

eas_oea1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene,
              eas_kor$gene, nfe_est$gene, nfe$gene, fin$gene, asj$gene)

nfe_est1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene, eas$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene,
              eas_kor$gene, eas_oea$gene,  fin$gene, asj$gene)

nfe1 <- c(afr$gene, sas$gene, amr$gene, eas$gene, eas_jpn$gene,
          eas_kor$gene, eas_oea$gene,  fin$gene, asj$gene)

fin1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene, eas$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene,
          eas_kor$gene, eas_oea$gene, nfe_est$gene, nfe$gene, asj$gene)

asj1 <- c(afr$gene, sas$gene, nfe_seu$gene, nfe_bgr$gene, nfe_onf$gene, amr$gene, eas$gene, nfe_swe$gene, nfe_nwe$gene, eas_jpn$gene,
          eas_kor$gene, eas_oea$gene, nfe_est$gene, fin$gene, nfe$gene)


afr.genes <- afr[!(afr$gene %in% afr1),]
afr.genes <- unique(afr.genes$gene)
afr.genes
[1] "AMER1"    "MSN"      "FAM47C"   "DDX5"     "MLLT3"    "TP53"     "LMNA"     "RPL22"    "ERC1"     "KIAA1549" "TCF7L2"   "CDKN1B"   "RANBP2"  
"TCEA1"    "PRPF40B"  "BRAF"     "ATRX"     "TSC2"     "CREB3L1" 
[20] "SIX1"     "ACKR3"    "SDHA"     "PDCD1LG2" "GPHN"     "ATP1A1"   "FOXO4"    "SRGAP3"   "HOXD13"   "STAT3"    "POLD1"    "MLLT1"    "BCL6"

gene <- cc.exome[cc.exome$gene == "BRAF",]
t.test(gene$afr, gene$nfe)
t.test(gene$afr, gene$eas)
t.test(gene$afr, gene$sas)
t.test(gene$afr, gene$fin)
t.test(gene$afr, gene$asj)

nfe_seu.genes <- nfe_seu[!(nfe_seu$gene %in% nfe_seu1),]
nfe_seu.genes <- unique(nfe_seu.genes$gene)
nfe_seu.genes
"ZNRF3"
gene <- nfe_seu[nfe_seu$gene == "ZNRF3",]
t.test(gene$afr, gene$nfe)
t.test(gene$afr, gene$eas)
t.test(gene$afr, gene$sas)
t.test(gene$afr, gene$fin)
t.test(gene$afr, gene$asj)


nfe_bgr.genes <- nfe_bgr[!(nfe_bgr$gene %in% nfe_bgr1),]
nfe_bgr.genes <- unique(nfe_bgr.genes$gene)
nfe_bgr.genes
nfe_bgr.genes 
"DAXX"   "FBLN2"  "POLG"   "ABL2"   "RARA"   "PTPN6"  "POT1"   "TSHR"   "CHST11" "PABPC1"

sas.genes <- sas[!(sas$gene %in% sas1),]
sas.genes <- unique(sas.genes$gene)
sas.genes
"CBL"      "RAD21"    "FUS"      "NCOA4"    "TERT"     "SFRP4"    "FIP1L1"   "CNTRL"    "AXIN1"    "NFIB"     "PRDM16"   "ABI1"     "ARNT"     "TNFRSF14" "ID3"      "CYLD"     "PTPRD"    "CSF3R"    "STAG2"

nfe_onf.genes <- nfe_onf[!(nfe_onf$gene %in% nfe_onf1),]
nfe_onf.genes <- unique(nfe_onf.genes$gene)
nfe_onf.genes
"SNX29"

amr.genes <- amr[!(amr$gene %in% amr1),]
amr.genes <- unique(amr.genes$gene)
amr.genes
[1] "MAF"      "EP300"    "ACVR1B"   "NONO"     "CTNNA2"   "PDGFRA"   "EXT2"     "ARHGEF12" "HNF1A"    "FGFR2"    "PRDM2"    "FOXL2"    "DCAF12L2" "TCF12"    "HMGA1"    "BCL9L"    "MDS2"     "FKBP9"    "DDX6"    
[20] "NUMA1"    "KMT2D"    "PAX3"     "KDR"      "DEK"      "BAZ1A"    "USP8"     "YWHAE"    "CPEB3"    "EXT1"     "MUTYH" 

eas.genes <- eas[!(eas$gene %in% eas1),]
eas.genes <- unique(eas.genes$gene)
eas.genes
"PTPN11"   "FANCA"    "ETNK1"    "SFPQ"     "TEC"      "SIRPA"    "RABEP1"   "STK11"    "PLCG1"    "CDC73"    "FLT4"     "TCF3"     "PIK3R1"   "NF2"      "HLF"      "ASXL1"    "HIST1H3B" "RAC1"     "PIK3CB"  
[20] "LATS1"

nfe_swe.genes <- nfe_swe[!(nfe_swe$gene %in% nfe_swe1),]
nfe_swe.genes <- unique(nfe_swe.genes$gene)
nfe_swe.genes
"CRLF2" "PTCH1" "MLF1" 

nfe_nwe.genes <- nfe_nwe[!(nfe_nwe$gene %in% nfe_nwe1),]
nfe_nwe.genes <- unique(nfe_nwe.genes$gene)
nfe_nwe.genes
"PAFAH1B2" "PDGFRB"

eas_jpn.genes <- eas_jpn[!(eas_jpn$gene %in% eas_jpn1),]
eas_jpn.genes <- unique(eas_jpn.genes$gene)
eas_jpn.genes
"CBFA2T3" "KLF6"    "TRIM27"  "SMAD3"   "SGK1"    "CDK4"    "CLIP1"   "NIN"     "PATZ1"   "LYN"     "ATF1"    "N4BP2"   "TPM3"    "ETV5" 

eas_kor.genes <- eas_kor[!(eas_kor$gene %in% eas_kor1),]
eas_kor.genes <- unique(eas_kor.genes$gene)
eas_kor.genes
"NCOA1"  "WIF1"   "RB1"    "RECQL4"

eas_oea.genes <- eas_oea[!(eas_oea$gene %in% eas_oea1),]
eas_oea.genes <- unique(eas_oea.genes$gene)
eas_oea.genes
"NF2"    "ASXL1"  "LATS1"  "PIK3CB" "SDC4" 

nfe_est.genes <- nfe_est[!(nfe_est$gene %in% nfe_est1),]
nfe_est.genes <- unique(nfe_est.genes$gene)
nfe_est.genes
"CTNNA1"   "LSM14A"   "JAK3"     "U2AF1"    "PTPRK"    "ARHGEF10" "BAX"      "ROS1" 

nfe.genes <- nfe[!(nfe$gene %in% nfe1),]
nfe.genes <- unique(nfe.genes$gene)
nfe.genes
"ZNF429"   "CBLB"     "PML"      "NCOA2"    "FAT1"     "ACSL6"    "MYH9"     "GNAQ"     "POLR2A"   "PAFAH1B2" "ESR1"     "NBEA"     "GLI1"     "ZNRF3"    "SNX29" 

fin.genes <- fin[!(fin$gene %in% fin1),]
fin.genes <- unique(fin.genes$gene)
fin.genes
"BIRC3"  "NACA"   "EPHA3"  "SETBP1" "BARD1"  "IL7R"   "CCND1"  "SKI"    "PRKACA" "FCRL4"  "IKZF3"  "MAP2K2"

asj.genes <- asj[!(asj$gene %in% asj1),]
asj.genes <- unique(asj.genes$gene)
asj.genes
"ZMYM2"  "USP9X"  "AKAP9"  "REL"    "TRIM33" "SND1" 


for(i in 1:17){
  current <- i+5
  exome.t <- cc.exome[order(cc.exome[,current], decreasing = TRUE),]
  exome <- exome.t[exome.t[,current] > 0,]
  help <- table(exome$gene)
  exome.n <- exome[1:50,]
  highest <- c(highest,exome.n$gene)
  mut.burden <- table(exome$gene)
  mut.burden <- mut.burden[order(mut.burden,decreasing = TRUE)]
  mut.burden <- mut.burden[1:50]
  mutated <- c(mutated, names(mut.burden) )
  
}


head(mutated)
head(highest)
top 100 mutation burden only from enhanced > 0 
tem <- table(mutated)
tem <- tem[order(tem, decreasing = TRUE)]

hihghest enahnceent 
what <- table(highest)
what <- what[order(what, decreasing = TRUE)]
head(what)


temp <- "afr"
afr.exome <- cc.exome[order(cc.exome$afr, decreasing = TRUE),]
afr.exome$gene[1:50]
afr.exome <- afr.exome[afr.exome$afr > 0 ,]
mean(afr.exome$afr) ## avg EAF: 0.005309681
range(afr.exome$afr) ## range: 1.688750e-08 8.146301e-01
afr.genes <- unique(afr.exome$gene)
length(afr.genes) # 707 gnees with at least 1 enhanced mutation 
afr.mut.burden <- table(afr.exome$gene)
afr.mut.burden <- afr.mut.burden[order(afr.mut.burden,decreasing = TRUE)]
afr.mut.burden[1:20]
afr.mut <- afr.mut.burden[1:100] 
afr.genes



## enriched in any populations) 
mean(cc.exome$nfe_seu)
mean(cc.exome$nfe_bgr)
mean(cc.exome$afr)
mean(cc.exome$sas)
mean(cc.exome$nfe_onf)
mean(cc.exome$amr)
mean(cc.exome$eas)
mean(cc.exome$nfe_swe)
mean(cc.exome$nfe_nwe)
mean(cc.exome$eas_jpn)
mean(cc.exome$eas_kor)
mean(cc.exome$eas_oea)
mean(cc.exome$nfe_est)
mean(cc.exome$nfe)
mean(cc.exome$fin)
mean(cc.exome$asj)

dim(cc.exome[cc.exome$nfe_seu > 0,])
dim(cc.exome[cc.exome$afr > 0,])
dim(cc.exome[cc.exome$nfe > 0,])



old.cc.exome <- cc.exome

cc.exome <- cc.exome[cc.exome]
################################################# Which ones are most commonly enriched? ############################################## 
common.enhanced <- eur.exome[eur.exome$gene %in% afr.exome$gene,]
common.enhanced <- common.enhanced[common.enhanced$gene %in% eas.exome$gene,]

### take top enriched genes from each each ancestry and see which ones are in common between them ###
eur.hold <- eur.exome$gene[1:200]
afr.hold <- afr.exome$gene[1:200]
eas.hold <- eas.exome$gene[1:200]
all.1 <- eur.hold[eur.hold %in% afr.hold]
all.2 <- all.1[all.1 %in% eas.hold]
unique(all.2)
#"MGMT"  "SEPT9" "GAS7"  "MUC4"  "EPAS1" "CUX1"  "NSD1" 

################################################# Are there drugs (PharmGKB lookup) that target those genes, and if so are they specific to a population or do they work in many? ############################################## 
eur.mutation <- eur.exome[1:200,]
eur.MGMT <-eur.mutation[eur.mutation$gene %in% "MGMT",]
afr.mutation <- afr.exome[1:200,]
eas.mutation <- eas.exome[1:200,]

################################################# % of the genes that are COSMIC but also EAF high in some populations - look at the altered alleles in those populations & predict (ExPecTo model at HB website) the functional impact of that alteration on the gene. ############################################## 
eas.top.20.mut <- eas.exome[1:20,]

afr.top.20.mut <- afr.exome[1:20,]

eur.top.20.mut <- eur.exome[1:20,]
################################################# Are those genes likely functionally altered in some populations but not others? ############################################## 








model.stats <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(model.stats) <- c("Model","Test_Samples","Tstat","Df","pValue","Metric", "Disease")

disease <- "uterine"
model <- "mixed"
ancestry.sample <- "mixed_ancestry"
ancestry.model <- "mixed"

temp.list <- list.files(paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/PAPER_RUNS/samples/",ancestry.sample,"/"))
num <- length(temp.list)
num

### AFR 
df.path <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/PAPER_RUNS/model_runs/",ancestry.model)
pf.path <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/PAPER_RUNS/model_runs/phyloFrame/",ancestry.model)
mod.num <- list(1:num)
mod.num <- unlist(mod.num)

# df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"_metrics.tsv") #CHANGE HERE
# pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"_metrics.tsv")#CHANGE HERE
df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"afr_metrics.tsv") #CHANGE HERE
pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"afr_metrics.tsv")#CHANGE HERE
myfiles <- lapply(df.names, readr::read_tsv)
pffiles <- lapply(pf.names, readr::read_tsv)

my.row <- paste0(ancestry.model,"_model_",mod.num)

# - for benchmark model - #
afr.mod <- do.call("rbind", myfiles)
afr.roc <- afr.mod[afr.mod$.metric == "roc_auc",]
afr.roc$model <- my.row
afr.acc <- afr.mod[afr.mod$.metric == "accuracy",]
afr.acc$model <- my.row
afr.sens<- afr.mod[afr.mod$.metric == "sens",]
afr.sens$model <- my.row
afr.spec <- afr.mod[afr.mod$.metric == "spec",]
afr.spec$model <- my.row
afr.prec <- afr.mod[afr.mod$.metric == "precision",]
afr.prec$model <- my.row
afr.recall <- afr.mod[afr.mod$.metric == "recall",]
afr.recall$model <- my.row
afr.f.meas <- afr.mod[afr.mod$.metric == "f_meas",]
afr.f.meas$model <- my.row
afr.kap <- afr.mod[afr.mod$.metric == "kap",]
afr.kap$model <- my.row
afr.mcc <- afr.mod[afr.mod$.metric == "mcc",]
afr.mcc$model <- my.row

# - phyloframe model - # 
pf.mod <- do.call("rbind", pffiles)
pf.roc <- pf.mod[pf.mod$.metric == "roc_auc",]
pf.roc$model <- my.row
pf.acc <- pf.mod[pf.mod$.metric == "accuracy",]
pf.acc$model <- my.row
pf.sens<- pf.mod[pf.mod$.metric == "sens",]
pf.sens$model <- my.row
pf.spec <- pf.mod[pf.mod$.metric == "spec",]
pf.spec$model <- my.row
pf.prec <- pf.mod[pf.mod$.metric == "precision",]
pf.prec$model <- my.row
pf.recall <- pf.mod[pf.mod$.metric == "recall",]
pf.recall$model <- my.row
pf.f.meas <- pf.mod[pf.mod$.metric == "f_meas",]
pf.f.meas$model <- my.row
pf.kap <- pf.mod[pf.mod$.metric == "kap",]
pf.kap$model <- my.row
pf.mcc <- pf.mod[pf.mod$.metric == "mcc",]
pf.mcc$model <- my.row

afr.roc.t <- t.test(afr.roc$.estimate, pf.roc$.estimate, var.equal = FALSE)
temp.list <- c(model,"afr", afr.roc.t$statistic, afr.roc.t$parameter, afr.roc.t$p.value, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
afr.acc.t <- t.test(afr.acc$.estimate, pf.acc$.estimate, var.equal = FALSE)
temp.list <- c(model,"afr", afr.acc.t$statistic, afr.acc.t$parameter, afr.acc.t$p.value, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
afr.sens.t <- t.test(afr.sens$.estimate, pf.sens$.estimate, var.equal = FALSE)
temp.list <- c(model,"afr", afr.sens.t$statistic, afr.sens.t$parameter, afr.sens.t$p.value, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
afr.spec.t <- t.test(afr.spec$.estimate, pf.spec$.estimate, var.equal = FALSE)
temp.list <- c(model,"afr", afr.spec.t$statistic, afr.spec.t$parameter, afr.spec.t$p.value, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
afr.prec.t <- t.test(afr.prec$.estimate, pf.prec$.estimate, var.equal = FALSE)
temp.list <- c(model,"afr", afr.prec.t$statistic, afr.prec.t$parameter, afr.prec.t$p.value, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
afr.recall.t <- t.test(afr.recall$.estimate, pf.recall$.estimate, var.equal = FALSE)
temp.list <- c(model,"afr", afr.recall.t$statistic, afr.recall.t$parameter, afr.recall.t$p.value, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
afr.f.meas.t <- t.test(afr.f.meas$.estimate, pf.f.meas$.estimate, var.equal = FALSE)
temp.list <- c(model,"afr", afr.f.meas.t$statistic, afr.f.meas.t$parameter, afr.f.meas.t$p.value, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
afr.kap.t <- t.test(afr.kap$.estimate, pf.kap$.estimate, var.equal = FALSE)
temp.list <- c(model,"afr", afr.kap.t$statistic, afr.kap.t$parameter, afr.kap.t$p.value, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
afr.mcc.t <- t.test(afr.mcc$.estimate, pf.mcc$.estimate, var.equal = FALSE)
temp.list <- c(model,"afr", afr.mcc.t$statistic, afr.mcc.t$parameter, afr.mcc.t$p.value, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats

# ## if not enough samples 
temp.list <- c(model,"afr", NA, NA, NA, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"afr", NA, NA, NA, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"afr", NA, NA, NA, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"afr", NA, NA, NA, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"afr", NA, NA, NA, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"afr", NA, NA, NA, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"afr", NA, NA, NA, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"afr", NA, NA, NA, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"afr", NA, NA, NA, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats

## ADMIXED 
df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"admixed_metrics.tsv") #CHANGE HERE
pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"admixed_metrics.tsv")#CHANGE HERE
myfiles <- lapply(df.names, readr::read_tsv)
pffiles <- lapply(pf.names, readr::read_tsv)

#- for benchmark model - #
admixed.mod <- do.call("rbind", myfiles)
admixed.roc <- admixed.mod[admixed.mod$.metric == "roc_auc",]
admixed.roc$model <- my.row
admixed.acc <- admixed.mod[admixed.mod$.metric == "accuracy",]
admixed.acc$model <- my.row
admixed.sens<- admixed.mod[admixed.mod$.metric == "sens",]
admixed.sens$model <- my.row
admixed.spec <- admixed.mod[admixed.mod$.metric == "spec",]
admixed.spec$model <- my.row
admixed.prec <- admixed.mod[admixed.mod$.metric == "precision",]
admixed.prec$model <- my.row
admixed.recall <- admixed.mod[admixed.mod$.metric == "recall",]
admixed.recall$model <- my.row
admixed.f.meas <- admixed.mod[admixed.mod$.metric == "f_meas",]
admixed.f.meas$model <- my.row
admixed.kap <- admixed.mod[admixed.mod$.metric == "kap",]
admixed.kap$model <- my.row
admixed.mcc <- admixed.mod[admixed.mod$.metric == "mcc",]
admixed.mcc$model <- my.row

#- phyloframe model -# 
pf.mod <- do.call("rbind", pffiles)
pf.roc <- pf.mod[pf.mod$.metric == "roc_auc",]
pf.roc$model <- my.row
pf.acc <- pf.mod[pf.mod$.metric == "accuracy",]
pf.acc$model <- my.row
pf.sens<- pf.mod[pf.mod$.metric == "sens",]
pf.sens$model <- my.row
pf.spec <- pf.mod[pf.mod$.metric == "spec",]
pf.spec$model <- my.row
pf.prec <- pf.mod[pf.mod$.metric == "precision",]
pf.prec$model <- my.row
pf.recall <- pf.mod[pf.mod$.metric == "recall",]
pf.recall$model <- my.row
pf.f.meas <- pf.mod[pf.mod$.metric == "f_meas",]
pf.f.meas$model <- my.row
pf.kap <- pf.mod[pf.mod$.metric == "kap",]
pf.kap$model <- my.row
pf.mcc <- pf.mod[pf.mod$.metric == "mcc",]
pf.mcc$model <- my.row

admixed.roc.t <- t.test(admixed.roc$.estimate, pf.roc$.estimate, var.equal = FALSE)
temp.list <- c(model,"admixed", admixed.roc.t$statistic, admixed.roc.t$parameter, admixed.roc.t$p.value, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
admixed.acc.t <- t.test(admixed.acc$.estimate, pf.acc$.estimate, var.equal = FALSE)
temp.list <- c(model,"admixed", admixed.acc.t$statistic, admixed.acc.t$parameter, admixed.acc.t$p.value, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
admixed.sens.t <- t.test(admixed.sens$.estimate, pf.sens$.estimate, var.equal = FALSE)
temp.list <- c(model,"admixed", admixed.sens.t$statistic, admixed.sens.t$parameter, admixed.sens.t$p.value, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
admixed.spec.t <- t.test(admixed.spec$.estimate, pf.spec$.estimate, var.equal = FALSE)
temp.list <- c(model,"admixed", admixed.spec.t$statistic, admixed.spec.t$parameter, admixed.spec.t$p.value, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
admixed.prec.t <- t.test(admixed.prec$.estimate, pf.prec$.estimate, var.equal = FALSE)
temp.list <- c(model,"admixed", admixed.prec.t$statistic, admixed.prec.t$parameter, admixed.prec.t$p.value, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
admixed.recall.t <- t.test(admixed.recall$.estimate, pf.recall$.estimate, var.equal = FALSE)
temp.list <- c(model,"admixed", admixed.recall.t$statistic, admixed.recall.t$parameter, admixed.recall.t$p.value, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
admixed.f.meas.t <- t.test(admixed.f.meas$.estimate, pf.f.meas$.estimate, var.equal = FALSE)
temp.list <- c(model,"admixed", admixed.f.meas.t$statistic, admixed.f.meas.t$parameter, admixed.f.meas.t$p.value, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
admixed.kap.t <- t.test(admixed.kap$.estimate, pf.kap$.estimate, var.equal = FALSE)
temp.list <- c(model,"admixed", admixed.kap.t$statistic, admixed.kap.t$parameter, admixed.kap.t$p.value, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
admixed.mcc.t <- t.test(admixed.mcc$.estimate, pf.mcc$.estimate, var.equal = FALSE)
temp.list <- c(model,"admixed", admixed.mcc.t$statistic, admixed.mcc.t$parameter, admixed.mcc.t$p.value, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats

temp.list <- c(model,"Admixed", NA, NA, NA, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Admixed", NA, NA, NA, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Admixed", NA, NA, NA, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Admixed", NA, NA, NA, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Admixed", NA, NA, NA, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Admixed", NA, NA, NA, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Admixed", NA, NA, NA, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Admixed", NA, NA, NA, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Admixed", NA, NA, NA, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats

## EAS
# df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"_metrics.tsv") #CHANGE HERE
# pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"_metrics.tsv")#CHANGE HERE
df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eas_metrics.tsv") #CHANGE HERE
pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"eas_metrics.tsv")#CHANGE HERE
myfiles <- lapply(df.names, readr::read_tsv)
pffiles <- lapply(pf.names, readr::read_tsv)

#- for benchmark model - #
eas.mod <- do.call("rbind", myfiles)
eas.roc <- eas.mod[eas.mod$.metric == "roc_auc",]
eas.roc$model <- my.row
eas.acc <- eas.mod[eas.mod$.metric == "accuracy",]
eas.acc$model <- my.row
eas.sens<- eas.mod[eas.mod$.metric == "sens",]
eas.sens$model <- my.row
eas.spec <- eas.mod[eas.mod$.metric == "spec",]
eas.spec$model <- my.row
eas.prec <- eas.mod[eas.mod$.metric == "precision",]
eas.prec$model <- my.row
eas.recall <- eas.mod[eas.mod$.metric == "recall",]
eas.recall$model <- my.row
eas.f.meas <- eas.mod[eas.mod$.metric == "f_meas",]
eas.f.meas$model <- my.row
eas.kap <- eas.mod[eas.mod$.metric == "kap",]
eas.kap$model <- my.row
eas.mcc <- eas.mod[eas.mod$.metric == "mcc",]
eas.mcc$model <- my.row

#- phyloframe model -# 
pf.mod <- do.call("rbind", pffiles)
pf.roc <- pf.mod[pf.mod$.metric == "roc_auc",]
pf.roc$model <- my.row
pf.acc <- pf.mod[pf.mod$.metric == "accuracy",]
pf.acc$model <- my.row
pf.sens<- pf.mod[pf.mod$.metric == "sens",]
pf.sens$model <- my.row
pf.spec <- pf.mod[pf.mod$.metric == "spec",]
pf.spec$model <- my.row
pf.prec <- pf.mod[pf.mod$.metric == "precision",]
pf.prec$model <- my.row
pf.recall <- pf.mod[pf.mod$.metric == "recall",]
pf.recall$model <- my.row
pf.f.meas <- pf.mod[pf.mod$.metric == "f_meas",]
pf.f.meas$model <- my.row
pf.kap <- pf.mod[pf.mod$.metric == "kap",]
pf.kap$model <- my.row
pf.mcc <- pf.mod[pf.mod$.metric == "mcc",]
pf.mcc$model <- my.row

eas.roc.t <- t.test(eas.roc$.estimate, pf.roc$.estimate, var.equal = FALSE)
temp.list <- c(model,"eas", eas.roc.t$statistic, eas.roc.t$parameter, eas.roc.t$p.value, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eas.acc.t <- t.test(eas.acc$.estimate, pf.acc$.estimate, var.equal = FALSE)
temp.list <- c(model,"eas", eas.acc.t$statistic, eas.acc.t$parameter, eas.acc.t$p.value, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eas.sens.t <- t.test(eas.sens$.estimate, pf.sens$.estimate, var.equal = FALSE)
temp.list <- c(model,"eas", eas.sens.t$statistic, eas.sens.t$parameter, eas.sens.t$p.value, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eas.spec.t <- t.test(eas.spec$.estimate, pf.spec$.estimate, var.equal = FALSE)
temp.list <- c(model,"eas", eas.spec.t$statistic, eas.spec.t$parameter, eas.spec.t$p.value, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eas.prec.t <- t.test(eas.prec$.estimate, pf.prec$.estimate, var.equal = FALSE)
temp.list <- c(model,"eas", eas.prec.t$statistic, eas.prec.t$parameter, eas.prec.t$p.value, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eas.recall.t <- t.test(eas.recall$.estimate, pf.recall$.estimate, var.equal = FALSE)
temp.list <- c(model,"eas", eas.recall.t$statistic, eas.recall.t$parameter, eas.recall.t$p.value, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eas.f.meas.t <- t.test(eas.f.meas$.estimate, pf.f.meas$.estimate, var.equal = FALSE)
temp.list <- c(model,"eas", eas.f.meas.t$statistic, eas.f.meas.t$parameter, eas.f.meas.t$p.value, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eas.kap.t <- t.test(eas.kap$.estimate, pf.kap$.estimate, var.equal = FALSE)
temp.list <- c(model,"eas", eas.kap.t$statistic, eas.kap.t$parameter, eas.kap.t$p.value, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eas.mcc.t <- t.test(eas.mcc$.estimate, pf.mcc$.estimate, var.equal = FALSE)
temp.list <- c(model,"eas", eas.mcc.t$statistic, eas.mcc.t$parameter, eas.mcc.t$p.value, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats

temp.list <- c(model,"Eas", NA, NA, NA, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eas", NA, NA, NA, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eas", NA, NA, NA, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eas", NA, NA, NA, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eas", NA, NA, NA, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eas", NA, NA, NA, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eas", NA, NA, NA, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eas", NA, NA, NA, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eas", NA, NA, NA, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats

## EUR
df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eur_metrics.tsv") #CHANGE HERE
pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"eur_metrics.tsv")#CHANGE HERE
myfiles <- lapply(df.names, readr::read_tsv)
pffiles <- lapply(pf.names, readr::read_tsv)

#- for benchmark model - #
eur.mod <- do.call("rbind", myfiles)
eur.roc <- eur.mod[eur.mod$.metric == "roc_auc",]
eur.roc$model <- my.row
eur.acc <- eur.mod[eur.mod$.metric == "accuracy",]
eur.acc$model <- my.row
eur.sens<- eur.mod[eur.mod$.metric == "sens",]
eur.sens$model <- my.row
eur.spec <- eur.mod[eur.mod$.metric == "spec",]
eur.spec$model <- my.row
eur.prec <- eur.mod[eur.mod$.metric == "precision",]
eur.prec$model <- my.row
eur.recall <- eur.mod[eur.mod$.metric == "recall",]
eur.recall$model <- my.row
eur.f.meas <- eur.mod[eur.mod$.metric == "f_meas",]
eur.f.meas$model <- my.row
eur.kap <- eur.mod[eur.mod$.metric == "kap",]
eur.kap$model <- my.row
eur.mcc <- eur.mod[eur.mod$.metric == "mcc",]
eur.mcc$model <- my.row

#- phyloframe model -# 
pf.mod <- do.call("rbind", pffiles)
pf.roc <- pf.mod[pf.mod$.metric == "roc_auc",]
pf.roc$model <- my.row
pf.acc <- pf.mod[pf.mod$.metric == "accuracy",]
pf.acc$model <- my.row
pf.sens<- pf.mod[pf.mod$.metric == "sens",]
pf.sens$model <- my.row
pf.spec <- pf.mod[pf.mod$.metric == "spec",]
pf.spec$model <- my.row
pf.prec <- pf.mod[pf.mod$.metric == "precision",]
pf.prec$model <- my.row
pf.recall <- pf.mod[pf.mod$.metric == "recall",]
pf.recall$model <- my.row
pf.f.meas <- pf.mod[pf.mod$.metric == "f_meas",]
pf.f.meas$model <- my.row
pf.kap <- pf.mod[pf.mod$.metric == "kap",]
pf.kap$model <- my.row
pf.mcc <- pf.mod[pf.mod$.metric == "mcc",]
pf.mcc$model <- my.row

eur.roc.t <- t.test(eur.roc$.estimate, pf.roc$.estimate, var.equal = FALSE)
temp.list <- c(model,"eur", eur.roc.t$statistic, eur.roc.t$parameter, eur.roc.t$p.value, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eur.acc.t <- t.test(eur.acc$.estimate, pf.acc$.estimate, var.equal = FALSE)
temp.list <- c(model,"eur", eur.acc.t$statistic, eur.acc.t$parameter, eur.acc.t$p.value, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eur.sens.t <- t.test(eur.sens$.estimate, pf.sens$.estimate, var.equal = FALSE)
temp.list <- c(model,"eur", eur.sens.t$statistic, eur.sens.t$parameter, eur.sens.t$p.value, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eur.spec.t <- t.test(eur.spec$.estimate, pf.spec$.estimate, var.equal = FALSE)
temp.list <- c(model,"eur", eur.spec.t$statistic, eur.spec.t$parameter, eur.spec.t$p.value, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eur.prec.t <- t.test(eur.prec$.estimate, pf.prec$.estimate, var.equal = FALSE)
temp.list <- c(model,"eur", eur.prec.t$statistic, eur.prec.t$parameter, eur.prec.t$p.value, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eur.recall.t <- t.test(eur.recall$.estimate, pf.recall$.estimate, var.equal = FALSE)
temp.list <- c(model,"eur", eur.recall.t$statistic, eur.recall.t$parameter, eur.recall.t$p.value, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eur.f.meas.t <- t.test(eur.f.meas$.estimate, pf.f.meas$.estimate, var.equal = FALSE)
temp.list <- c(model,"eur", eur.f.meas.t$statistic, eur.f.meas.t$parameter, eur.f.meas.t$p.value, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eur.kap.t <- t.test(eur.kap$.estimate, pf.kap$.estimate, var.equal = FALSE)
temp.list <- c(model,"eur", eur.kap.t$statistic, eur.kap.t$parameter, eur.kap.t$p.value, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
eur.mcc.t <- t.test(eur.mcc$.estimate, pf.mcc$.estimate, var.equal = FALSE)
temp.list <- c(model,"eur", eur.mcc.t$statistic, eur.mcc.t$parameter, eur.mcc.t$p.value, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats


temp.list <- c(model,"Eur", NA, NA, NA, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eur", NA, NA, NA, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eur", NA, NA, NA, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eur", NA, NA, NA, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eur", NA, NA, NA, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eur", NA, NA, NA, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eur", NA, NA, NA, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eur", NA, NA, NA, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Eur", NA, NA, NA, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats


## MIXED
df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"mixed_metrics.tsv") #CHANGE HERE
pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"mixed_metrics.tsv")#CHANGE HERE
myfiles <- lapply(df.names, readr::read_tsv)
pffiles <- lapply(pf.names, readr::read_tsv)

#- for benchmark model - #
mixed.mod <- do.call("rbind", myfiles)
mixed.roc <- mixed.mod[mixed.mod$.metric == "roc_auc",]
mixed.roc$model <- my.row
mixed.acc <- mixed.mod[mixed.mod$.metric == "accuracy",]
mixed.acc$model <- my.row
mixed.sens<- mixed.mod[mixed.mod$.metric == "sens",]
mixed.sens$model <- my.row
mixed.spec <- mixed.mod[mixed.mod$.metric == "spec",]
mixed.spec$model <- my.row
mixed.prec <- mixed.mod[mixed.mod$.metric == "precision",]
mixed.prec$model <- my.row
mixed.recall <- mixed.mod[mixed.mod$.metric == "recall",]
mixed.recall$model <- my.row
mixed.f.meas <- mixed.mod[mixed.mod$.metric == "f_meas",]
mixed.f.meas$model <- my.row
mixed.kap <- mixed.mod[mixed.mod$.metric == "kap",]
mixed.kap$model <- my.row
mixed.mcc <- mixed.mod[mixed.mod$.metric == "mcc",]
mixed.mcc$model <- my.row

#- phyloframe model -# 
pf.mod <- do.call("rbind", pffiles)
pf.roc <- pf.mod[pf.mod$.metric == "roc_auc",]
pf.roc$model <- my.row
pf.acc <- pf.mod[pf.mod$.metric == "accuracy",]
pf.acc$model <- my.row
pf.sens<- pf.mod[pf.mod$.metric == "sens",]
pf.sens$model <- my.row
pf.spec <- pf.mod[pf.mod$.metric == "spec",]
pf.spec$model <- my.row
pf.prec <- pf.mod[pf.mod$.metric == "precision",]
pf.prec$model <- my.row
pf.recall <- pf.mod[pf.mod$.metric == "recall",]
pf.recall$model <- my.row
pf.f.meas <- pf.mod[pf.mod$.metric == "f_meas",]
pf.f.meas$model <- my.row
pf.kap <- pf.mod[pf.mod$.metric == "kap",]
pf.kap$model <- my.row
pf.mcc <- pf.mod[pf.mod$.metric == "mcc",]
pf.mcc$model <- my.row

mixed.roc.t <- t.test(mixed.roc$.estimate, pf.roc$.estimate, var.equal = FALSE)
temp.list <- c(model,"mixed", mixed.roc.t$statistic, mixed.roc.t$parameter, mixed.roc.t$p.value, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
mixed.acc.t <- t.test(mixed.acc$.estimate, pf.acc$.estimate, var.equal = FALSE)
temp.list <- c(model,"mixed", mixed.acc.t$statistic, mixed.acc.t$parameter, mixed.acc.t$p.value, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
mixed.sens.t <- t.test(mixed.sens$.estimate, pf.sens$.estimate, var.equal = FALSE)
temp.list <- c(model,"mixed", mixed.sens.t$statistic, mixed.sens.t$parameter, mixed.sens.t$p.value, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
mixed.spec.t <- t.test(mixed.spec$.estimate, pf.spec$.estimate, var.equal = FALSE)
temp.list <- c(model,"mixed", mixed.spec.t$statistic, mixed.spec.t$parameter, mixed.spec.t$p.value, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
mixed.prec.t <- t.test(mixed.prec$.estimate, pf.prec$.estimate, var.equal = FALSE)
temp.list <- c(model,"mixed", mixed.prec.t$statistic, mixed.prec.t$parameter, mixed.prec.t$p.value, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
mixed.recall.t <- t.test(mixed.recall$.estimate, pf.recall$.estimate, var.equal = FALSE)
temp.list <- c(model,"mixed", mixed.recall.t$statistic, mixed.recall.t$parameter, mixed.recall.t$p.value, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
mixed.f.meas.t <- t.test(mixed.f.meas$.estimate, pf.f.meas$.estimate, var.equal = FALSE)
temp.list <- c(model,"mixed", mixed.f.meas.t$statistic, mixed.f.meas.t$parameter, mixed.f.meas.t$p.value, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
mixed.kap.t <- t.test(mixed.kap$.estimate, pf.kap$.estimate, var.equal = FALSE)
temp.list <- c(model,"mixed", mixed.kap.t$statistic, mixed.kap.t$parameter, mixed.kap.t$p.value, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
mixed.mcc.t <- t.test(mixed.mcc$.estimate, pf.mcc$.estimate, var.equal = FALSE)
temp.list <- c(model,"mixed", mixed.mcc.t$statistic, mixed.mcc.t$parameter, mixed.mcc.t$p.value, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats

temp.list <- c(model,"Mixed", NA, NA, NA, "AUC", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Mixed", NA, NA, NA, "Accuracy", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Mixed", NA, NA, NA, "Sensitivity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Mixed", NA, NA, NA, "Specificity", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Mixed", NA, NA, NA, "Precision", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Mixed", NA, NA, NA, "Recall", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Mixed", NA, NA, NA, "Fmeas", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Mixed", NA, NA, NA, "Kap", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats
temp.list <- c(model,"Mixed", NA, NA, NA, "Mcc", disease)
model.stats[nrow(model.stats) + 1,] <- temp.list
model.stats


save.incase <- model.stats

#breast <- model.stats
#thyroid <- model.stats
uterine <- model.stats

write.table(uterine, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/stats/all_anc_ttest_uterine_models.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
d


########## 
#median of ancestry
estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)

## do for each ancestry - median 
afr <- estimated_ancestry[estimated_ancestry$consensus_ancestry == "afr",]
dim(afr)
tt <- table(afr$consensus_ancestry, afr$tumor_type)
median(tt)

amr <- estimated_ancestry[estimated_ancestry$consensus_ancestry == "amr",]
amr.tt <- table(amr$consensus_ancestry, amr$tumor_type)
median(amr.tt)

eas <- estimated_ancestry[estimated_ancestry$consensus_ancestry == "eas",]
eas.tt <- table(eas$consensus_ancestry, eas$tumor_type)
median(eas.tt)

eur <- estimated_ancestry[estimated_ancestry$consensus_ancestry == "eur",]
eur.tt <- table(eur$consensus_ancestry, eur$tumor_type)
median(eur.tt)

sas <- estimated_ancestry[estimated_ancestry$consensus_ancestry == "sas",]
sas.tt <- table(sas$consensus_ancestry, sas$tumor_type)
median(sas.tt)

admixed <- c("admix","afr_admix","eas_admix","eur_admix","sas_admix")
admix <- estimated_ancestry[estimated_ancestry$consensus_ancestry %in% admixed,]
admix.tt <- table(admix$consensus_ancestry, admix$tumor_type)
median(admix.tt)

all <- table(estimated_ancestry$tumor_type,estimated_ancestry$consensus_ancestry)
dat <- as.data.frame(all)
head(dat)
sum(dat$Freq) # total samples = 10678 

sum(dat$Freq[dat$Var2 == "afr"])
tbl.all <- table(estimated_ancestry$consensus_ancestry)
total <- 10678
afr <- 651
admix <- 343 + 7 + 68 + 24 + 12 
amr <- 41
eas <- 669 
eur <- 8836
sas <- 27 

(afr/total)*100
(admix/total)*100
(amr/total)*100
(eas/total)*100
(eur/total)*100
(sas/total)*100

tbl.all <- table(estimated_ancestry$tumor_type)
dat.all <- as.data.frame(tbl.all)
head(dat.all)
eur.t <- estimated_ancestry[estimated_ancestry$consensus_ancestry == "eur",]

eur.all <- table(eur.t$tumor_type)
eur.all <- as.data.frame(eur.all)
head(eur.all)
colnames(eur.all) <- c("Var1", "eur_num")

eur.percent <- merge(dat.all, eur.all, by = "Var1")
head(eur.percent)
eur.percent$percent <- (eur.percent$eur_num/ eur.percent$Freq)*100
head(eur.percent)
eur.percent
range(eur.percent$percent)
tbale(lapply(c("eur","afr","amr","eas","sas"), function(i){temp <-estimated_ancestry[estimated_ancestry$consensus_ancestry == i,] tt <- table(afr$consensus_ancestry, afr$tumor_type)
median(tt)})
