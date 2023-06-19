### ENHANCED ALLELE FREQUENCY CALCULATION FOR GNOAMD EXOME DATA ###
require(readr)  # for read_csv()
require(dplyr) #data frame handling 
library(tidyr)
library(stringi) #string mnutations 

annotated_genes <- "/home/leslie.smith1/blue_kgraim/leslie.smith1/old_phyloFrame/preprocessing/gnomad_exome_parsed_all_ancestries.tsv"
genes <- readr::read_delim(annotated_genes, delim = '\t',col_names = FALSE) 
##for exome 
colnames(genes) <- c("chrom","position","rs_id","ref_allele", "alt_allele",
                     "AC_amr_male","AN_amr_male","AF_amr_male",
                     "AC_amr_female","AN_amr_female","AF_amr_female",
                     "AC_eas_male","AN_eas_male", "AF_eas_male",
                     "AC_eas_female","AN_eas_female","AF_eas_female",
                     "AC_nfe_male","AN_nfe_male", "AF_nfe_male",
                     "AC_nfe_female","AN_nfe_female", "AF_nfe_female",
                     "AC_asj_male","AN_asj_male","AF_asj_male",
                     "AC_asj_female","AN_asj_female", "AF_asj_female",
                     "AC_oth_male","AN_oth_male", "AF_oth_male",
                     "AC_oth_female","AN_oth_female","AF_oth_female",
                     "AC_fin_male","AN_fin_male","AF_fin_male",
                     "AC_fin_female","AN_fin_female", "AF_fin_female",
                     "AC_afr_male","AN_afr_male", "AF_afr_male",
                     "AC_afr_female","AN_afr_female", "AF_afr_female",
                     "AC_sas_male","AN_sas_male", "AF_sas_male",
                     "AC_sas_female","AN_sas_female","AF_sas_female",
                     "AC_nfe_seu", "AN_nfe_seu","AF_nfe_seu",
                     "AC_nfe_bgr",	"AN_nfe_bgr",	"AF_nfe_bgr",
                     "AC_afr", "AN_afr",	"AF_afr",
                     "AC_sas",	"AN_sas",	"AF_sas",
                     "AC_nfe_onf",	"AN_nfe_onf",	"AF_nfe_onf",
                     "AC_amr",	"AN_amr",	"AF_amr",
                     "AC_eas",	"AN_eas",	"AF_eas",
                     "AC_nfe_swe",	"AN_nfe_swe",	"AF_nfe_swe",
                     "AC_nfe_nwe",	"AN_nfe_nwe",	"AF_nfe_nwe",
                     "AC_eas_jpn",	"AN_eas_jpn",	"AF_eas_jpn",
                     "AC_eas_kor",	"AN_eas_kor",	"AF_eas_kor",
                     "AC_eas_oea",	"AN_eas_oea",	"AF_eas_oea",
                     "AC_nfe_est",	"AN_nfe_est",	"AF_nfe_est",
                     "AC_nfe",	"AN_nfe",	"AF_nfe",
                     "AC_fin",	"AN_fin",	"AF_fin",
                     "AC_asj",	"AN_asj",	"AF_asj",
                     "AC_oth",	"AN_oth",	"AF_oth",
                     "gene","consequence","type", "distance")

## - keep only the actual AC, AN, or AF from the file to use in EAF calculation - ## 
genes$AC_amr_male <- stringr::str_remove(genes$AC_amr_male,"AC_amr_male=" )
genes$AC_amr_male <- as.numeric(genes$AC_amr_male)
genes$AN_amr_male <- stringr::str_remove(genes$AN_amr_male,"AN_amr_male=" )
genes$AN_amr_male <- as.numeric(genes$AN_amr_male)
genes$AF_amr_male <- stringr::str_remove(genes$AF_amr_male,"AF_amr_male=" )
genes$AF_amr_male <- as.numeric(genes$AF_amr_male)

genes$AC_amr_female <- stringr::str_remove(genes$AC_amr_female,"AC_amr_female=" )
genes$AC_amr_female <- as.numeric(genes$AC_amr_female)
genes$AN_amr_female <- stringr::str_remove(genes$AN_amr_female,"AN_amr_female=" )
genes$AN_amr_female <- as.numeric(genes$AN_amr_female)
genes$AF_amr_female <- stringr::str_remove(genes$AF_amr_female,"AF_amr_female=" )
genes$AF_amr_female <- as.numeric(genes$AF_amr_female)

genes$AC_eas_male <- stringr::str_remove(genes$AC_eas_male,"AC_eas_male=" )
genes$AC_eas_male <- as.numeric(genes$AC_eas_male)
genes$AN_eas_male <- stringr::str_remove(genes$AN_eas_male,"AN_eas_male=" )
genes$AN_eas_male <- as.numeric(genes$AN_eas_male)
genes$AF_eas_male <- stringr::str_remove(genes$AF_eas_male,"AF_eas_male=" )
genes$AF_eas_male <- as.numeric(genes$AF_eas_male)

genes$AC_eas_female <- stringr::str_remove(genes$AC_eas_female,"AC_eas_female=" )
genes$AC_eas_female <- as.numeric(genes$AC_eas_female)
genes$AN_eas_female <- stringr::str_remove(genes$AN_eas_female,"AN_eas_female=" )
genes$AN_eas_female <- as.numeric(genes$AN_eas_female)
genes$AF_eas_female <- stringr::str_remove(genes$AF_eas_female,"AF_eas_female=" )
genes$AF_eas_female <- as.numeric(genes$AF_eas_female)

genes$AC_nfe_male <- stringr::str_remove(genes$AC_nfe_male,"AC_nfe_male=" )
genes$AC_nfe_male <- as.numeric(genes$AC_nfe_male)
genes$AN_nfe_male <- stringr::str_remove(genes$AN_nfe_male,"AN_nfe_male=" )
genes$AN_nfe_male <- as.numeric(genes$AN_nfe_male)
genes$AF_nfe_male <- stringr::str_remove(genes$AF_nfe_male,"AF_nfe_male=" )
genes$AF_nfe_male <- as.numeric(genes$AF_nfe_male)

genes$AC_nfe_female <- stringr::str_remove(genes$AC_nfe_female,"AC_nfe_female=" )
genes$AC_nfe_female <- as.numeric(genes$AC_nfe_female)
genes$AN_nfe_female <- stringr::str_remove(genes$AN_nfe_female,"AN_nfe_female=" )
genes$AN_nfe_female <- as.numeric(genes$AN_nfe_female)
genes$AF_nfe_female <- stringr::str_remove(genes$AF_nfe_female,"AF_nfe_female=" )
genes$AF_nfe_female <- as.numeric(genes$AF_nfe_female)

genes$AC_asj_male <- stringr::str_remove(genes$AC_asj_male,"AC_asj_male=" )
genes$AC_asj_male <- as.numeric(genes$AC_asj_male)
genes$AN_asj_male <- stringr::str_remove(genes$AN_asj_male,"AN_asj_male=" )
genes$AN_asj_male <- as.numeric(genes$AN_asj_male)
genes$AF_asj_male <- stringr::str_remove(genes$AF_asj_male,"AF_asj_male=" )
genes$AF_asj_male <- as.numeric(genes$AF_asj_male)

genes$AC_asj_female <- stringr::str_remove(genes$AC_asj_female,"AC_asj_female=" )
genes$AC_asj_female <- as.numeric(genes$AC_asj_female)
genes$AN_asj_female <- stringr::str_remove(genes$AN_asj_female,"AN_asj_female=" )
genes$AN_asj_female <- as.numeric(genes$AN_asj_female)
genes$AF_asj_female <- stringr::str_remove(genes$AF_asj_female,"AF_asj_female=" )
genes$AF_asj_female <- as.numeric(genes$AF_asj_female)

genes$AC_oth_male <- stringr::str_remove(genes$AC_oth_male,"AC_oth_male=" )
genes$AC_oth_male <- as.numeric(genes$AC_oth_male)
genes$AN_oth_male <- stringr::str_remove(genes$AN_oth_male,"AN_oth_male=" )
genes$AN_oth_male <- as.numeric(genes$AN_oth_male)
genes$AF_oth_male <- stringr::str_remove(genes$AF_oth_male,"AF_oth_male=" )
genes$AF_oth_male <- as.numeric(genes$AF_oth_male)

genes$AC_oth_female <- stringr::str_remove(genes$AC_oth_female,"AC_oth_female=" )
genes$AC_oth_female <- as.numeric(genes$AC_oth_female)
genes$AN_oth_female <- stringr::str_remove(genes$AN_oth_female,"AN_oth_female=" )
genes$AN_oth_female <- as.numeric(genes$AN_oth_female)
genes$AF_oth_female <- stringr::str_remove(genes$AF_oth_female,"AF_oth_female=" )
genes$AF_oth_female <- as.numeric(genes$AF_oth_female)

genes$AC_fin_male <- stringr::str_remove(genes$AC_fin_male,"AC_fin_male=" )
genes$AC_fin_male <- as.numeric(genes$AC_fin_male)
genes$AN_fin_male <- stringr::str_remove(genes$AN_fin_male,"AN_fin_male=" )
genes$AN_fin_male <- as.numeric(genes$AN_fin_male)
genes$AF_fin_male <- stringr::str_remove(genes$AF_fin_male,"AF_fin_male=" )
genes$AF_fin_male <- as.numeric(genes$AF_fin_male)

genes$AC_fin_female <- stringr::str_remove(genes$AC_fin_female,"AC_fin_female=" )
genes$AC_fin_female <- as.numeric(genes$AC_fin_female)
genes$AN_fin_female <- stringr::str_remove(genes$AN_fin_female,"AN_fin_female=" )
genes$AN_fin_female <- as.numeric(genes$AN_fin_female)
genes$AF_fin_female <- stringr::str_remove(genes$AF_fin_female,"AF_fin_female=" )
genes$AF_fin_female <- as.numeric(genes$AF_fin_female)

genes$AC_afr_male <- stringr::str_remove(genes$AC_afr_male,"AC_afr_male=" )
genes$AC_afr_male <- as.numeric(genes$AC_afr_male)
genes$AN_afr_male <- stringr::str_remove(genes$AN_afr_male,"AN_afr_male=" )
genes$AN_afr_male <- as.numeric(genes$AN_afr_male)
genes$AF_afr_male <- stringr::str_remove(genes$AF_afr_male,"AF_afr_male=" )
genes$AF_afr_male <- as.numeric(genes$AF_afr_male)

genes$AC_afr_female <- stringr::str_remove(genes$AC_afr_female,"AC_afr_female=" )
genes$AC_afr_female <- as.numeric(genes$AC_afr_female)
genes$AN_afr_female <- stringr::str_remove(genes$AN_afr_female,"AN_afr_female=" )
genes$AN_afr_female <- as.numeric(genes$AN_afr_female)
genes$AF_afr_female <- stringr::str_remove(genes$AF_afr_female,"AF_afr_female=" )
genes$AF_afr_female <- as.numeric(genes$AF_afr_female)

genes$AC_sas_male <- stringr::str_remove(genes$AC_sas_male,"AC_sas_male=" )
genes$AC_sas_male <- as.numeric(genes$AC_sas_male)
genes$AN_sas_male <- stringr::str_remove(genes$AN_sas_male,"AN_sas_male=" )
genes$AN_sas_male <- as.numeric(genes$AN_sas_male)
genes$AF_sas_male <- stringr::str_remove(genes$AF_sas_male,"AF_sas_male=" )
genes$AF_sas_male <- as.numeric(genes$AF_sas_male)

genes$AC_sas_female <- stringr::str_remove(genes$AC_sas_female,"AC_sas_female=" )
genes$AC_sas_female <- as.numeric(genes$AC_sas_female)
genes$AN_sas_female <- stringr::str_remove(genes$AN_sas_female,"AN_sas_female=" )
genes$AN_sas_female <- as.numeric(genes$AN_sas_female)
genes$AF_sas_female <- stringr::str_remove(genes$AF_sas_female,"AF_sas_female=" )
genes$AF_sas_female <- as.numeric(genes$AF_sas_female)

#####NONSEPCIFIC SEX ANCESTRIES #####
genes$AC_nfe_seu <- stringr::str_remove(genes$AC_nfe_seu,"AC_nfe_seu=" )
genes$AC_nfe_seu <- as.numeric(genes$AC_nfe_seu)
genes$AN_nfe_seu <- stringr::str_remove(genes$AN_nfe_seu,"AN_nfe_seu=" )
genes$AN_nfe_seu <- as.numeric(genes$AN_nfe_seu)
genes$AF_nfe_seu <- stringr::str_remove(genes$AF_nfe_seu,"AF_nfe_seu=" )
genes$AF_nfe_seu <- as.numeric(genes$AF_nfe_seu)

genes$AC_nfe_bgr <- stringr::str_remove(genes$AC_nfe_bgr,"AC_nfe_bgr=" )
genes$AC_nfe_bgr <- as.numeric(genes$AC_nfe_bgr)
genes$AN_nfe_bgr <- stringr::str_remove(genes$AN_nfe_bgr,"AN_nfe_bgr=" )
genes$AN_nfe_bgr <- as.numeric(genes$AN_nfe_bgr)
genes$AF_nfe_bgr <- stringr::str_remove(genes$AF_nfe_bgr,"AF_nfe_bgr=" )
genes$AF_nfe_bgr <- as.numeric(genes$AF_nfe_bgr)

genes$AC_afr <- stringr::str_remove(genes$AC_afr,"AC_afr=" )
genes$AC_afr <- as.numeric(genes$AC_afr)
genes$AN_afr <- stringr::str_remove(genes$AN_afr,"AN_afr=" )
genes$AN_afr <- as.numeric(genes$AN_afr)
genes$AF_afr <- stringr::str_remove(genes$AF_afr,"AF_afr=" )
genes$AF_afr <- as.numeric(genes$AF_afr)

genes$AC_sas <- stringr::str_remove(genes$AC_sas,"AC_sas=" )
genes$AC_sas <- as.numeric(genes$AC_sas)
genes$AN_sas <- stringr::str_remove(genes$AN_sas,"AN_sas=" )
genes$AN_sas <- as.numeric(genes$AN_sas)
genes$AF_sas <- stringr::str_remove(genes$AF_sas,"AF_sas=" )
genes$AF_sas <- as.numeric(genes$AF_sas)

genes$AC_nfe_onf <- stringr::str_remove(genes$AC_nfe_onf,"AC_nfe_onf=" )
genes$AC_nfe_onf <- as.numeric(genes$AC_nfe_onf)
genes$AN_nfe_onf <- stringr::str_remove(genes$AN_nfe_onf,"AN_nfe_onf=" )
genes$AN_nfe_onf <- as.numeric(genes$AN_nfe_onf)
genes$AF_nfe_onf <- stringr::str_remove(genes$AF_nfe_onf,"AF_nfe_onf=" )
genes$AF_nfe_onf <- as.numeric(genes$AF_nfe_onf)

genes$AC_amr <- stringr::str_remove(genes$AC_amr,"AC_amr=" )
genes$AC_amr <- as.numeric(genes$AC_amr)
genes$AN_amr <- stringr::str_remove(genes$AN_amr,"AN_amr=" )
genes$AN_amr <- as.numeric(genes$AN_amr)
genes$AF_amr <- stringr::str_remove(genes$AF_amr,"AF_amr=" )
genes$AF_amr <- as.numeric(genes$AF_amr)

genes$AC_eas <- stringr::str_remove(genes$AC_eas,"AC_eas=" )
genes$AC_eas <- as.numeric(genes$AC_eas)
genes$AN_eas <- stringr::str_remove(genes$AN_eas,"AN_eas=" )
genes$AN_eas <- as.numeric(genes$AN_eas)
genes$AF_eas <- stringr::str_remove(genes$AF_eas,"AF_eas=" )
genes$AF_eas <- as.numeric(genes$AF_eas)

genes$AC_nfe_swe <- stringr::str_remove(genes$AC_nfe_swe,"AC_nfe_swe=" )
genes$AC_nfe_swe <- as.numeric(genes$AC_nfe_swe)
genes$AN_nfe_swe <- stringr::str_remove(genes$AN_nfe_swe,"AN_nfe_swe=" )
genes$AN_nfe_swe <- as.numeric(genes$AN_nfe_swe)
genes$AF_nfe_swe <- stringr::str_remove(genes$AF_nfe_swe,"AF_nfe_swe=" )
genes$AF_nfe_swe <- as.numeric(genes$AF_nfe_swe)

genes$AC_nfe_nwe <- stringr::str_remove(genes$AC_nfe_nwe,"AC_nfe_nwe=" )
genes$AC_nfe_nwe <- as.numeric(genes$AC_nfe_nwe)
genes$AN_nfe_nwe <- stringr::str_remove(genes$AN_nfe_nwe,"AN_nfe_nwe=" )
genes$AN_nfe_nwe <- as.numeric(genes$AN_nfe_nwe)
genes$AF_nfe_nwe <- stringr::str_remove(genes$AF_nfe_nwe,"AF_nfe_nwe=" )
genes$AF_nfe_nwe <- as.numeric(genes$AF_nfe_nwe)

genes$AC_eas_jpn <- stringr::str_remove(genes$AC_eas_jpn,"AC_eas_jpn=" )
genes$AC_eas_jpn <- as.numeric(genes$AC_eas_jpn)
genes$AN_eas_jpn <- stringr::str_remove(genes$AN_eas_jpn,"AN_eas_jpn=" )
genes$AN_eas_jpn <- as.numeric(genes$AN_eas_jpn)
genes$AF_eas_jpn <- stringr::str_remove(genes$AF_eas_jpn,"AF_eas_jpn=" )
genes$AF_eas_jpn <- as.numeric(genes$AF_eas_jpn)

genes$AC_eas_kor <- stringr::str_remove(genes$AC_eas_kor,"AC_eas_kor=" )
genes$AC_eas_kor <- as.numeric(genes$AC_eas_kor)
genes$AN_eas_kor <- stringr::str_remove(genes$AN_eas_kor,"AN_eas_kor=" )
genes$AN_eas_kor <- as.numeric(genes$AN_eas_kor)
genes$AF_eas_kor <- stringr::str_remove(genes$AF_eas_kor,"AF_eas_kor=" )
genes$AF_eas_kor <- as.numeric(genes$AF_eas_kor)

genes$AC_eas_oea <- stringr::str_remove(genes$AC_eas_oea,"AC_eas_oea=" )
genes$AC_eas_oea <- as.numeric(genes$AC_eas_oea)
genes$AN_eas_oea <- stringr::str_remove(genes$AN_eas_oea,"AN_eas_oea=" )
genes$AN_eas_oea <- as.numeric(genes$AN_eas_oea)
genes$AF_eas_oea <- stringr::str_remove(genes$AF_eas_oea,"AF_eas_oea=" )
genes$AF_eas_oea <- as.numeric(genes$AF_eas_oea)

genes$AC_nfe_est <- stringr::str_remove(genes$AC_nfe_est,"AC_nfe_est=" )
genes$AC_nfe_est <- as.numeric(genes$AC_nfe_est)
genes$AN_nfe_est <- stringr::str_remove(genes$AN_nfe_est,"AN_nfe_est=" )
genes$AN_nfe_est <- as.numeric(genes$AN_nfe_est)
genes$AF_nfe_est <- stringr::str_remove(genes$AF_nfe_est,"AF_nfe_est=" )
genes$AF_nfe_est <- as.numeric(genes$AF_nfe_est)

genes$AC_nfe <- stringr::str_remove(genes$AC_nfe,"AC_nfe=" )
genes$AC_nfe <- as.numeric(genes$AC_nfe)
genes$AN_nfe <- stringr::str_remove(genes$AN_nfe,"AN_nfe=" )
genes$AN_nfe <- as.numeric(genes$AN_nfe)
genes$AF_nfe <- stringr::str_remove(genes$AF_nfe,"AF_nfe=" )
genes$AF_nfe <- as.numeric(genes$AF_nfe)

genes$AC_fin <- stringr::str_remove(genes$AC_fin,"AC_fin=" )
genes$AC_fin <- as.numeric(genes$AC_fin)
genes$AN_fin <- stringr::str_remove(genes$AN_fin,"AN_fin=" )
genes$AN_fin <- as.numeric(genes$AN_fin)
genes$AF_fin <- stringr::str_remove(genes$AF_fin,"AF_fin=" )
genes$AF_fin <- as.numeric(genes$AF_fin)

genes$AC_asj <- stringr::str_remove(genes$AC_asj,"AC_asj=" )
genes$AC_asj <- as.numeric(genes$AC_asj)
genes$AN_asj <- stringr::str_remove(genes$AN_asj,"AN_asj=" )
genes$AN_asj <- as.numeric(genes$AN_asj)
genes$AF_asj <- stringr::str_remove(genes$AF_asj,"AF_asj=" )
genes$AF_asj <- as.numeric(genes$AF_asj)

genes$AC_oth <- stringr::str_remove(genes$AC_oth,"AC_oth=" )
genes$AC_oth <- as.numeric(genes$AC_oth)
genes$AN_oth <- stringr::str_remove(genes$AN_oth,"AN_oth=" )
genes$AN_oth <- as.numeric(genes$AN_oth)
genes$AF_oth <- stringr::str_remove(genes$AF_oth,"AF_oth=" )
genes$AF_oth <- as.numeric(genes$AF_oth)


##CREATE TABLE TO PUT FREQUENCIES IN AFTER CALCULATIONS 
no_sex <- genes %>% dplyr::select(chrom, position, rs_id, ref_allele, alt_allele,AF_nfe_seu, AF_nfe_bgr, AF_afr, AF_sas, AF_nfe_onf, AF_amr, AF_eas,	AF_nfe_swe, AF_nfe_nwe,	
                                  AF_eas_jpn, AF_eas_kor,AF_eas_oea, AF_nfe_est, AF_nfe,	AF_fin,	AF_asj,	AF_oth,
                                  gene,consequence,type, distance)

colnames(no_sex) <- c("chrom","position","rs_id","ref_allele", "alt_allele","nfe_seu","nfe_bgr", "afr", "sas", "nfe_onf", "amr",	
                      "eas",	"nfe_swe", "nfe_nwe", "eas_jpn", "eas_kor",	"eas_oea", "nfe_est", "nfe", "fin", "asj", "oth","gene","consequence","type", "distance")


#####NONE SEX SPECIFIC ENHANCED #################

no_sex$nfe_seu <- (genes$AF_nfe_seu - ((genes$AF_nfe_bgr+ genes$AF_afr+ genes$AF_sas+ genes$AF_nfe_onf+ genes$AF_amr+ genes$AF_eas+ genes$AF_nfe_swe+ genes$AF_nfe_nwe+ 
                                          genes$AF_eas_jpn+ genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16)                 )


##NFE BGR 
no_sex$nfe_bgr <- (genes$AF_nfe_bgr - ((genes$AF_nfe_seu+ genes$AF_afr+ genes$AF_sas+ genes$AF_nfe_onf+ genes$AF_amr+ genes$AF_eas+ genes$AF_nfe_swe+ genes$AF_nfe_nwe+ 
                                          genes$AF_eas_jpn+ genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))


##AFR 
no_sex$afr <- (genes$AF_afr - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_sas+ genes$AF_nfe_onf+ genes$AF_amr+ genes$AF_eas+ genes$AF_nfe_swe+ genes$AF_nfe_nwe+ 
                                  genes$AF_eas_jpn+ genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16)) 


##SAS
no_sex$sas <- (genes$AF_sas - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_nfe_onf+ genes$AF_amr+ genes$AF_eas+ genes$AF_nfe_swe+ genes$AF_nfe_nwe+ 
                                  genes$AF_eas_jpn+ genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))  

##NFE ONF
no_sex$nfe_onf <- (genes$AF_nfe_onf - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_amr+ genes$AF_eas+ genes$AF_nfe_swe+ genes$AF_nfe_nwe+ 
                                          genes$AF_eas_jpn+ genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16)) 

##AMR
##NFE ONF
no_sex$amr <- (genes$AF_amr - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf+ genes$AF_eas+ genes$AF_nfe_swe+ genes$AF_nfe_nwe+ 
                                  genes$AF_eas_jpn+ genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))

##EAS 
no_sex$eas <- (genes$AF_eas - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr+ genes$AF_nfe_swe+ genes$AF_nfe_nwe+ 
                                  genes$AF_eas_jpn+ genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))


##NFE_SWE
no_sex$nfe_swe <- (genes$AF_nfe_swe - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_nwe+ 
                                          genes$AF_eas_jpn+ genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))


##NFE_NWE
no_sex$nfe_nwe <- (genes$AF_nfe_nwe - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_swe + 
                                          genes$AF_eas_jpn+ genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))

##EAS_JPN
no_sex$eas_jpn <- (genes$AF_eas_jpn - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_swe + 
                                          genes$AF_nfe_nwe + genes$AF_eas_kor+ genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))


##EAS-KOR
no_sex$eas_kor <- (genes$AF_eas_kor - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_swe + 
                                          genes$AF_nfe_nwe + genes$AF_eas_jpn + genes$AF_eas_oea+ genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))


##EAS OEA 
no_sex$eas_oea <- (genes$AF_eas_oea - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_swe + 
                                          genes$AF_nfe_nwe + genes$AF_eas_jpn + genes$AF_eas_kor + genes$AF_nfe_est+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))


##NFE EST 
no_sex$nfe_est <- (genes$AF_nfe_est - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_swe + 
                                          genes$AF_nfe_nwe + genes$AF_eas_jpn + genes$AF_eas_kor + genes$AF_eas_oea+ genes$AF_nfe+ genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))


##NFE
no_sex$nfe <- (genes$AF_nfe - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_swe + 
                                  genes$AF_nfe_nwe + genes$AF_eas_jpn + genes$AF_eas_kor + genes$AF_eas_oea + genes$AF_nfe_est + genes$AF_fin+ genes$AF_asj+ genes$AF_oth)/16))

##FIN
no_sex$fin <- (genes$AF_fin - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_swe + 
                                  genes$AF_nfe_nwe + genes$AF_eas_jpn + genes$AF_eas_kor + genes$AF_eas_oea + genes$AF_nfe_est + genes$AF_nfe + genes$AF_asj + genes$AF_oth)/16))

##ASJ
no_sex$asj <- (genes$AF_asj - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_swe + 
                                  genes$AF_nfe_nwe + genes$AF_eas_jpn + genes$AF_eas_kor + genes$AF_eas_oea + genes$AF_nfe_est + genes$AF_nfe + genes$AF_fin + genes$AF_oth)/16))


##OTH
no_sex$oth <- (genes$AF_oth - ((genes$AF_nfe_seu + genes$AF_nfe_bgr + genes$AF_afr + genes$AF_sas + genes$AF_nfe_onf + genes$AF_amr + genes$AF_eas + genes$AF_nfe_swe + 
                                  genes$AF_nfe_nwe + genes$AF_eas_jpn + genes$AF_eas_kor + genes$AF_eas_oea + genes$AF_nfe_est + genes$AF_nfe + genes$AF_fin + genes$AF_asj)/16))

  

write.table(no_sex, "/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/mean_enhancedAF_exome.tsv", sep='\t', col.names=T, row.names=F, quote=F)



