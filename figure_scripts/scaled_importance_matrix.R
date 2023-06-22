#FIGURE NOT IN PAPER - USED IN PIPELINE OUTPUT
#############################################################################
#SCRIPT FOR GENE WEIGHT COMPARISON BETWEEN BENCHMARK AND PHYLOFRAME SIGNATURES. 
#############################################################################
library(pheatmap)
library(grid)

get.matrix <- function(disease, work.dir, plot.dir){
  dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/",disease,"/",work.dir,"/model_runs/")
  
  eur.num <- length(list.files(paste0(dir,"eur")))
  eur.num <- list(1:eur.num)
  eur.num <- unlist(eur.num)
  afr.num <- length(list.files(paste0(dir,"afr")))
  afr.num <- list(1:afr.num)
  afr.num <- unlist(afr.num)
  eas.num <- length(list.files(paste0(dir,"eas")))
  eas.num <- list(1:eas.num)
  eas.num <- unlist(eas.num)
  admixed.num <- length(list.files(paste0(dir,"admixed")))
  admixed.num <- list(1:admixed.num)
  admixed.num <- unlist(admixed.num)
  mixed.num <- length(list.files(paste0(dir,"mixed")))
  mixed.num <- list(1:mixed.num)
  mixed.num <- unlist(mixed.num)
  #############################################################################
  
  eur.count <- length(eur.num)
  afr.count <- length(afr.num)
  eas.count <- length(eas.num)
  admixed.count <- length(admixed.num)
  mixed.count <- length(mixed.num)
  
  eur.model.names <- paste0("eur_model_",eur.num)
  afr.model.names <- paste0("afr_model_",afr.num)
  eas.model.names <- paste0("eas_model_",eas.num)
  admixed.names <- paste0("admixed_",admixed.num)  
  mixed.names <- paste0("mixed_model_", mixed.num)
  
  ## read in signatures from phyloframe and bechmark ancestry models 
  pf.eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/eur/")
  pf.afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/afr/")
  pf.admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/admixed/")
  pf.eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/eas/")
  pf.mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/phyloFrame/mixed/")
  
  eur.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/eur/")
  afr.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/afr/")
  admixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/admixed/")
  eas.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/eas/")
  mixed.dir <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",work.dir,"/model_runs/mixed/")
  
  pf.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt") 
  pf.eur.myfiles <- lapply(pf.eur.names, readr::read_tsv)
  pf.afr.names <- paste0(pf.afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt") 
  pf.afr.myfiles <- lapply(pf.afr.names, readr::read_tsv)
  pf.eas.names <- paste0(pf.eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt") 
  pf.eas.myfiles <- lapply(pf.eas.names, readr::read_tsv)
  pf.admixed.names <- paste0(pf.admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt") 
  pf.admixed.myfiles <- lapply(pf.admixed.names, readr::read_tsv)
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
  
  ############################### MATRIX CREATION ##############################################
  
  ## get list of all genes 
  
  
  add.model.names <- function(anc.list, ancestry, pf){
    dat <- data.frame(matrix(ncol = 5, nrow = 0))
    if(pf == TRUE){
      for (i in 1:length(anc.list)) {
        temp.dat <- anc.list[i]
        temp.dat <- as.data.frame(temp.dat)
        class(temp.dat)
        model.type <- rep(paste0("pf.",ancestry,".model.",i), nrow(temp.dat))
        temp.dat$model.type <- model.type
        temp.dat$scaled <- (temp.dat$Importance - min(temp.dat$Importance))/(max(temp.dat$Importance) - min(temp.dat$Importance))
        dat <- rbind(dat, temp.dat)
      }
      dat <- dat %>% dplyr::select(c(Variable, model.type, scaled))
    }else{
      for (i in 1:length(anc.list)) {
        temp.dat <- anc.list[i]
        temp.dat <- as.data.frame(temp.dat)
        class(temp.dat)
        model.type <- rep(paste0(ancestry,".model.",i), nrow(temp.dat))
        temp.dat$model.type <- model.type
        temp.dat$scaled <- (temp.dat$Importance - min(temp.dat$Importance))/(max(temp.dat$Importance) - min(temp.dat$Importance))
        dat <- rbind(dat, temp.dat)
        
      }
      dat <- dat %>% dplyr::select(c(Variable, model.type, scaled))
    }
    return(dat)
  }
  pf.eur.dat <- add.model.names(pf.eur.myfiles,"eur", TRUE)  ## repeat this for every ancestry
  pf.afr.dat <- add.model.names(pf.afr.myfiles,"afr", TRUE)
  pf.eas.dat <- add.model.names(pf.eas.myfiles,"eas", TRUE)
  pf.admixed.dat <- add.model.names(pf.admixed.myfiles,"admixed", TRUE)
  pf.mixed.dat <- add.model.names(pf.mixed.myfiles,"mixed", TRUE)
  
  eur.dat <- add.model.names(eur.myfiles,"eur", FALSE)  ## repeat this for every ancestry
  afr.dat <- add.model.names(afr.myfiles,"afr", FALSE)
  eas.dat <- add.model.names(eas.myfiles,"eas", FALSE)
  admixed.dat <- add.model.names(admixed.myfiles,"admixed", FALSE)
  mixed.dat <- add.model.names(mixed.myfiles,"mixed", FALSE)
  
  ##### MATRIX FOR ALL PF AND BENCHMARK #### 
  ## get list of all genes
  all.genes <- unique(c(pf.eur.dat$Variable, pf.afr.dat$Variable, pf.eas.dat$Variable, pf.admixed.dat$Variable, pf.mixed.dat$Variable,
                        eur.dat$Variable, afr.dat$Variable, eas.dat$Variable, admixed.dat$Variable, mixed.dat$Variable))
  all.dataframes <- rbind(pf.eur.dat, pf.afr.dat, pf.eas.dat, pf.admixed.dat, pf.mixed.dat, eur.dat, afr.dat, eas.dat, admixed.dat, mixed.dat)
#  all.dataframes <- rbind(eur.dat, afr.dat, eas.dat, admixed.dat, mixed.dat)
  
  tt <- all.dataframes %>% spread(Variable, scaled)
  tt[is.na(tt)] <- -1
  tt <- column_to_rownames(tt, "model.type")
  tt.mat <- as.matrix(tt)
  tt.mat[1:5,1:5]
  
  anc.names  <- c(rep("admixed",length(admixed.num)), rep("afr", length(afr.num)),
                            rep("eas",length(eas.num)), rep("eur",length(eur.num)), rep("mixed", length(mixed.num)))
  annotation <- data.frame(ancestry = c(anc.names, anc.names))
  rownames(annotation) <- rownames(tt.mat)
  
  heatmap <- pheatmap(
    tt.mat,
    legend_breaks = c(0, 0.5, 1),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    # annotation_colors = the.color,
    # annotation_col = for.annotation,
    annotation_row = annotation,
    main = "Gene Model Importance",
    fontsize_row = 5,
    fontsize_col = 5,
    colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(25),
    #grid.text("xlabel example", y=-0.07, gp=gpar(fontsize=16))
  )
  # 
  # 
  #heatmap
  png(paste0(plot.dir,disease, "gene_importance_all.png"),width = 1200, height = 1000)
  print(heatmap)
  dev.off()
  
  #### ONLY PHYLOFRAME ### 
  ##### MATRIX FOR ALL PF AND BENCHMARK #### 
  ## get list of all genes
  all.genes <- unique(c(pf.eur.dat$Variable, pf.afr.dat$Variable, pf.eas.dat$Variable, pf.admixed.dat$Variable, pf.mixed.dat$Variable))
  all.dataframes <- rbind(pf.eur.dat, pf.afr.dat, pf.eas.dat, pf.admixed.dat, pf.mixed.dat)
  #  all.dataframes <- rbind(eur.dat, afr.dat, eas.dat, admixed.dat, mixed.dat)
  
  tt <- all.dataframes %>% spread(Variable, scaled)
  tt[is.na(tt)] <- -1
  tt <- column_to_rownames(tt, "model.type")
  tt.mat <- as.matrix(tt)
  tt.mat[1:5,1:5]
  
  anc.names  <- c(rep("admixed",length(admixed.num)), rep("afr", length(afr.num)),
                  rep("eas",length(eas.num)), rep("eur",length(eur.num)), rep("mixed", length(mixed.num)))
  annotation <- data.frame(ancestry = anc.names)
  rownames(annotation) <- rownames(tt.mat)
  
  heatmap <- pheatmap(
    tt.mat,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    # annotation_colors = the.color,
    # annotation_col = for.annotation,
    annotation_row = annotation,
    main = "Gene Model Importance",
    fontsize_row = 5,
    fontsize_col = 5,
    colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(25),
    #grid.text("xlabel example", y=-0.07, gp=gpar(fontsize=16))
  )
  # 
  # 
  #heatmap
  png(paste0(plot.dir,disease, "gene_importance_phyloFrame.png"),width = 1200, height = 1000)
  print(heatmap)
  dev.off()
  
  #### GET BENCHMARK ### 
  ##### MATRIX FOR ALL PF AND BENCHMARK #### 
  ## get list of all genes
  all.genes <- unique(c(eur.dat$Variable, afr.dat$Variable, eas.dat$Variable, admixed.dat$Variable, mixed.dat$Variable))
  all.dataframes <- rbind(eur.dat, afr.dat, eas.dat, admixed.dat, mixed.dat)
  #  all.dataframes <- rbind(eur.dat, afr.dat, eas.dat, admixed.dat, mixed.dat)
  
  tt <- all.dataframes %>% spread(Variable, scaled)
  tt[is.na(tt)] <- -1
  tt <- column_to_rownames(tt, "model.type")
  tt.mat <- as.matrix(tt)
  tt.mat[1:5,1:5]
  
  anc.names  <- c(rep("admixed",length(admixed.num)), rep("afr", length(afr.num)),
                  rep("eas",length(eas.num)), rep("eur",length(eur.num)), rep("mixed", length(mixed.num)))
  annotation <- data.frame(ancestry = anc.names)
  rownames(annotation) <- rownames(tt.mat)
  
  heatmap <- pheatmap(
    tt.mat,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    # annotation_colors = the.color,
    # annotation_col = for.annotation,
    annotation_row = annotation,
    main = "Gene Model Importance",
    fontsize_row = 5,
    fontsize_col = 5,
    colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(25),
    #grid.text("xlabel example", y=-0.07, gp=gpar(fontsize=16))
  )
  # 
  # 
  #heatmap
  png(paste0(plot.dir,disease, "gene_importance_benchmark.png"),width = 1200, height = 1000)
  print(heatmap)
  dev.off()
}

############ EUR ############
## define lists for number of models in each ancestry for help reading in files ## 

