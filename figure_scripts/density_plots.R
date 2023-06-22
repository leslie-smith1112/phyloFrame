
###### FIGURE SUPP 4 AND DENISTY PLOTS FOR AACR2023 ######## 

library(tidyr)
library(ggplot2)
exomeAF <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/preprocessing/mean_enhancedAF_exome.tsv", col_names = TRUE)
census <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/cosmic_genes/Census_allWed Nov 2 01 58 37 2022.tsv", col_names = TRUE)
head(census)
census_genes <- census$`Gene Symbol`

data.long <- gather(exomeAF,ancestry, EAF, nfe_seu:oth, factor_key = TRUE)

head(data.long)
data.long$ancestry <- factor(data.long$ancestry, levels = c("afr", "amr", "asj","eas", "eas_jpn", 
                                                            "eas_kor" ,"eas_oea", "fin", "nfe", "nfe_bgr", "nfe_est","nfe_nwe","nfe_onf","nfe_seu", "nfe_swe","oth", "sas"))
e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.7) + theme_minimal() + scale_x_continuous(limits = c(-0.00005, .00025)) + scale_y_continuous(limits = c(0, 100000))

e+scale_fill_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#238B45", "#41AB5D","#74C476","#CE1256",
                             "#CA4136","#67000D", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#825283"))  + scale_color_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#238B45", "#41AB5D","#74C476","#CE1256",
                                                                                                                                                           "#CA4136","#67000D", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#825283")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")

#plot COSMIC genes EAF for only afr, eas, and nfe becasue those are the primary ancestries in analysis 

#keep only geens that are also in the  COSMIC dataset
cosmicAF <- exomeAF[exomeAF$gene %in% census_genes,]
cosmicAf <- cosmicAF %>% dplyr::select(nfe,afr,eas)

#get data in format for plotting
data.long <- gather(cosmicAf,ancestry, EAF, nfe_seu:oth, factor_key = TRUE)

e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.7) + theme_minimal() + scale_x_continuous(limits = c(-0.00005, .00025)) + scale_y_continuous(limits = c(0, 500))

e+scale_fill_manual(values=c("#1F619E","#496849","#CA4136"))  + scale_color_manual(values=c("#1F619E","#496849","#CA4136"))





### FOR AACR ### 
cosmicAFA <- exomeAF[exomeAF$gene %in% afr.sig,]

afr_sig_snps <- exomeAF[exomeAF$gene %in% afr$Variable,] #22958    26
afr_snp <- afr_sig_snps %>% dplyr::select(gene,afr)
afr_snp$ancestry <- "afr"
ggplot(afr_snp, aes(x = afr, fill = ancestry)) + geom_density(alpha = 0.5) + scale_x_continuous(limits = c(-0.00005, .00008)) + labs(subtitle = "Afr Samples Signature with Mean Enhanced Allele Frequency")


ggplot(cosmicAFA, aes(x = afr )) + geom_density( color = "#80B1D3",fill="#80B1D3") + scale_x_continuous(limits = c(-0.00005, .00008)) +
  theme_minimal() + coord_cartesian(ylim=c(0,280000)) # +scale_y_continuous(trans='log10')

ggplot(cosmicAF, aes(x = eas )) + geom_density( color = "#31A354",fill="#31A354") + 
  scale_x_continuous(limits = c(-0.00005, .00008)) + theme_minimal() + coord_cartesian(ylim=c(0,280000))

ggplot(cosmicAF, aes(x = nfe )) + geom_density( color = "#DE2D26",fill="#DE2D26") + scale_x_continuous(limits = c(-0.00005, .00008)) + 
  theme_minimal() + coord_cartesian(ylim=c(0,280000))

temp <- cosmicAF %>% select(c(afr,eas,nfe))

data.long <- gather(temp,ancestry, EAF, afr:nfe, factor_key = TRUE)

#color = c("#80B1D3","#31A354","#DE2D26"),fill=c("#80B1D3","#31A354","#DE2D26"))
ggplot(data.long, aes(x = EAF, fill = ancestry )) + geom_density() + scale_x_continuous(limits = c(-0.00005, .00008)) +geom_density(alpha = 0.5)+
  theme_minimal() + coord_cartesian(ylim=c(0,280000)) +scale_fill_manual( values = c("#80B1D3","#31A354","#DE2D26"))

afr.sorted <- cosmicAF[order(cosmicAF$afr, decreasing = TRUE),]
afr.max <- afr.sorted[1,]

data_long <- gather(afr.max, ancestry, EAF, nfe_seu:oth, factor_key=TRUE)

eas.sorted <- cosmicAF[order(cosmicAF$eas, decreasing = TRUE),]
eas.max <- eas.sorted[1,]

data_long <- gather(eas.max, ancestry, EAF, nfe_seu:oth, factor_key=TRUE)

nfe.sorted <- cosmicAF[order(cosmicAF$nfe, decreasing = TRUE),]
nfe.max <- nfe.sorted[1,]

data_long <- gather(nfe.max, ancestry, EAF, nfe_seu:oth, factor_key=TRUE)

test <- cosmicAF[1,]
data_long <- gather(test, ancestry, EAF, nfe_seu:oth, factor_key=TRUE)

temp <- exomeAF %>% dplyr::select(-chrom, -position, -rs_id, -ref_allele, -alt_allele, -gene, -consequence, -type, -distance)
data.long <- gather()

## density plot to demonstrate variation in EAF across ancestries for certain snps
data_long$ancestry <- factor(data_long$ancestry, levels = c("afr", "amr", "asj","eas", "eas_jpn", 
"eas_kor" ,"eas_oea", "fin", "nfe", "nfe_bgr", "nfe_est","nfe_nwe","nfe_onf","nfe_seu", "nfe_swe","oth", "sas"))

e <- ggplot(data = data_long, aes(x=ancestry, y = EAF, fill = ancestry)) + geom_bar(stat = "identity") + theme_minimal()
e+scale_fill_manual(values=c("#80B1D3","#F781BF","#B2ABD2","#005A32","#238B45", "#41AB5D","#74C476","#CE1256",
                                      "#67000D", "#A50F15","#CB181D", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#FFD92F"
                                      )) + coord_cartesian(ylim=c(0,1))

e

#get colors
display.brewer.all()
library(RColorBrewer)
colorpat <- "Purples"
display.brewer.pal(n = 8, name = colorpat)
brewer.pal(n = 9, name = colorpat)

