#### GRAPHING MODEL PEFORMANCE #### 
#### THESE ARE THE SCATTERPLOT PHYLOFRAME VS BENCHMARK FOR ALL BATCHES IN EACH MODEL ## 
require(tibble)
require(readr)  # for read_csv()
require(dplyr) #data frame handling 
require(tidyr)
require(stringi)
library(ggplot2)
set.seed(1234)
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/scaled_importance_matrix.R")
#disease <- "thyroid"
#ancestry.sample.mixed <- "e"
#ancestry.sample <- "eur"
#ancestry.model <- "eur"
#dir <- "pen1_version1_3"
#run(disease, dir, ancestry.sample, ancestry.model, 0.5)
# 
# 
# 
# 
#run(disease, dir, ancestry.sample, ancestry.model, mixture)

run <- function(disease, dir, ancestry.sample, ancestry.model, mixture){
  temp.list <- list.files(paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/samples/",ancestry.sample,"/"))
  num <- length(temp.list)
  num
  
  all <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(all) <- c("model.num", "ancestry", "benchmark", "phyloFrame")
  
  df.path <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",dir,"/model_runs/",ancestry.model)
  pf.path <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",dir,"/model_runs/phyloFrame/",ancestry.model)
  
  mod.num <- list(1:num)
  mod.num <- unlist(mod.num)
  
  main.test <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"_metrics.tsv")
  df.main.test <- paste0(df.path,"/model_",mod.num,"/model_", mod.num,"_metrics.tsv")
  main.file <- readr::read_tsv(main.test)
  df.main.file <- readr::read_tsv(df.main.test)
  pf.main.file <- main.file[main.file$.metric == "roc_auc",]
  df.main.file <- df.main.file[df.main.file$.metric == "roc_auc",]
  model.num <- paste0("model_",mod.num)
  ancestry <- rep(ancestry.model, num) #CHANGE HERE AND BELOW 
  main.tests <- data.frame("model.num" = model.num, "ancestry" = ancestry, "benchmark" = df.main.file$.estimate, "phyloFrame" = pf.main.file$.estimate)
  
  
  if(!(ancestry.sample == "afr")){
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"afr_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"afr_metrics.tsv")#CHANGE HERE
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    # - for benchmark model - #
    afr.mod <- do.call("rbind", myfiles)
    afr.mod <- afr.mod[afr.mod$.metric == "roc_auc",]
    
    # - phyloframe model - # 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("afr", num) #CHANGE HERE AND BELOW 
    new.afr.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                             "benchmark" = afr.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.afr.df)
  }else if(num >1){
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"afr_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"afr_metrics.tsv")#CHANGE HERE
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    # - for benchmark model - #
    afr.mod <- do.call("rbind", myfiles)
    afr.mod <- afr.mod[afr.mod$.metric == "roc_auc",]
    
    # - phyloframe model - # 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("afr", num) #CHANGE HERE AND BELOW 
    new.afr.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                             "benchmark" = afr.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.afr.df)
  }else{
    print("No conditions met")
  }
  
  if(!(ancestry.sample == "admixed")){
    
    ###ADMIXED 
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"admixed_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"admixed_metrics.tsv")#CHANGE HERE
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    #- for benchmark model - #
    admixed.mod <- do.call("rbind", myfiles)
    admixed.mod <- admixed.mod[admixed.mod$.metric == "roc_auc",]
    
    #- phyloframe model -# 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("admixed", num) #CHANGE HERE AND BELOW 
    new.admixed.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                                 "benchmark" = admixed.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.admixed.df)
  } else if(num > 1){
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"admixed_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"admixed_metrics.tsv")#CHANGE HERE
    
    
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    #- for benchmark model - #
    admixed.mod <- do.call("rbind", myfiles)
    admixed.mod <- admixed.mod[admixed.mod$.metric == "roc_auc",]
    
    #- phyloframe model -# 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("admixed", num) #CHANGE HERE AND BELOW 
    new.admixed.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                                 "benchmark" = admixed.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.admixed.df)
  }else{
    print("ERROR")
  }
  
  if(!(ancestry.sample == "eas")){
    ####EAS 
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eas_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"eas_metrics.tsv")#CHANGE HERE
    
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    #- for benchmark model - #
    eas.mod <- do.call("rbind", myfiles)
    eas.mod <- eas.mod[eas.mod$.metric == "roc_auc",]
    #- phyloframe model -# 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("eas", num) #CHANGE HERE AND BELOW 
    new.eas.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                             "benchmark" = eas.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.eas.df)
    
  }else if(num >1){
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eas_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"eas_metrics.tsv")#CHANGE HERE
    
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    #- for benchmark model - #
    eas.mod <- do.call("rbind", myfiles)
    eas.mod <- eas.mod[eas.mod$.metric == "roc_auc",]
    #- phyloframe model -# 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("eas", num) #CHANGE HERE AND BELOW 
    new.eas.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                             "benchmark" = eas.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.eas.df)
  }else{
    print("ERROR")
  }
  
  if(!(ancestry.sample == "eur")){
    #### EUR 
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eur_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"eur_metrics.tsv")#CHANGE HERE
    
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    #- for benchmark model - #
    eur.mod <- do.call("rbind", myfiles)
    eur.mod <- eur.mod[eur.mod$.metric == "roc_auc",]
    #- phyloframe model -# 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("eur", num) #CHANGE HERE AND BELOW 
    new.eur.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                             "benchmark" = eur.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.eur.df)
    
  }else if(num >1){
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"eur_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"eur_metrics.tsv")#CHANGE HERE
    
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    #- for benchmark model - #
    eur.mod <- do.call("rbind", myfiles)
    eur.mod <- eur.mod[eur.mod$.metric == "roc_auc",]
    #- phyloframe model -# 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("eur", num) #CHANGE HERE AND BELOW 
    new.eur.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                             "benchmark" = eur.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.eur.df)
  }else {
    print("ERROR")
  }
  
  if(!(ancestry.sample == "mixed")){
    ####MIXED 
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"mixed_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"mixed_metrics.tsv")#CHANGE HERE
    
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    #- for benchmark model - #
    mixed.mod <- do.call("rbind", myfiles)
    mixed.mod <- mixed.mod[mixed.mod$.metric == "roc_auc",]
    #- phyloframe model -# 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("mixed", num) #CHANGE HERE AND BELOW 
    new.mixed.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                               "benchmark" = mixed.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.mixed.df)
  }else if(num > 1){
    df.names <- paste0(df.path,"/model_",mod.num,"/model_",mod.num,"mixed_metrics.tsv") #CHANGE HERE
    pf.names <- paste0(pf.path,"/model_",mod.num,"/model_",mod.num,"mixed_metrics.tsv")#CHANGE HERE
    
    myfiles <- lapply(df.names, readr::read_tsv)
    pffiles <- lapply(pf.names, readr::read_tsv)
    
    #- for benchmark model - #
    mixed.mod <- do.call("rbind", myfiles)
    mixed.mod <- mixed.mod[mixed.mod$.metric == "roc_auc",]
    #- phyloframe model -# 
    pf.mod <- do.call("rbind", pffiles)
    pf.mod <- pf.mod[pf.mod$.metric == "roc_auc",]
    
    model.num <- paste0("model_",mod.num)
    ancestry <- rep("mixed", num) #CHANGE HERE AND BELOW 
    new.mixed.df <- data.frame("model.num" = model.num, "ancestry" = ancestry, 
                               "benchmark" = mixed.mod$.estimate, "phyloFrame"= pf.mod$.estimate)
    all <- rbind(all, new.mixed.df)
  }else {
    print("ERROR")
  }
  #all.a <- rbind(all, main.tests)
  all.a <- all
  df.path <- paste0("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/", disease,"/",dir,"/plots/")
  dir.create(df.path)
  get.matrix(disease, dir, df.path)
  #sometimes plots do better in vector graphics
  
  # make plot
  
  
  # png(file= paste0(df.path,ancestry.model, ".png"),
  #     width=600, height=350)
  p <- ggplot(all, aes(x=benchmark, y=phyloFrame, color=ancestry)) + geom_point(alpha = 5/10) + ggtitle(paste0(disease," models trained on ",num, " batches of ", ancestry.model," samples, penalty ",mixture)) +
    scale_color_manual(values = c("afr" = "blue",
                                  "admixed"="orange",
                                  "eas"="green",
                                  "eur" = "red",
                                  "mixed" = "purple")) + xlim(0,1) + ylim(0,1) + geom_abline()
  #p
  png(paste0(df.path,ancestry.model, "model_performance.png"),width = 800, height = 600)
  print(p)
  dev.off()
  # dev.off()
  
  
  #### boxplots ####
  library(data.table)
  # all ancestries 
  bench.dat <- data.frame("roc_auc" = all.a$benchmark, "model" = rep("becnhmark",nrow(all.a)))
  pf.dat <- data.frame("roc_auc" = all.a$phyloFrame, "model" = rep("phyloFrame",nrow(all.a)))
  all.dat <- rbind(bench.dat, pf.dat)
  
  
  png(paste0(df.path,ancestry.model, "all_boxplot.png"),width = 800, height = 600)
  tt <- boxplot(roc_auc ~ model,all.dat, main = paste0(disease, " ",ancestry.model," Benchmark vs PhyloFrame All"))
  recordedPlot = recordPlot()
  #print(tt)
  dev.off()
  
  pf.only <- data.frame("roc_auc" = all.a$phyloFrame, "ancestry" = all.a$ancestry)
  png(paste0(df.path,ancestry.model, "phyloFrame_boxplot.png"),width = 800, height = 600)
  hi <- boxplot(roc_auc ~ ancestry,pf.only, main = paste0(disease, " ",ancestry.model," PhyloFrame"))
  recordedPlot = recordPlot()
  #print(hi)
  dev.off()
  
  bench.only <- data.frame("roc_auc" = all.a$benchmark, "ancestry" = all.a$ancestry)
  png(paste0(df.path,ancestry.model, "benchmark_boxplot.png"),width = 800, height = 600)
  bye <- boxplot(roc_auc ~ ancestry,bench.only, main = paste0(disease, " ",ancestry.model," Benchmark"))
  recordedPlot = recordPlot()
  dev.off()

}
