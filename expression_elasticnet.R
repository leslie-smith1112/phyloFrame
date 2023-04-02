#BRCA expression data PhyloFrame 
library(tidymodels)
library(workflows) #tidy models package for bundling model specs 
library(parsnip) #for modeling 
library(workflows) #put model in workflow 
require(readr)  # for read_csv()
require(dplyr) #data frame handling 
library(tidyr)
library(stringi) #string mutations 
library(rpart.plot)
library(vip)
set.seed(4831)

######################################################################################################
# ## PHYLOFRAME PREPROCESSING BEFORE ELASTICNET RUN + ELASTICNET RUN##
######################################################################################################
## -- set up output directory -- ##
set.directory <- function(in.dir){
  if(!(file.exists(in.dir))){
    dir.create(in.dir)
  }
}

## expects matrix with samples as column names 
## -- trim expression matrix to the selected samples and genes for current run --  ##
## result is genes x sample 
trim.expr.matrix <- function(expression_matrix, include_genes, include_samples){
  if(!(is.null(include_genes)))
  {
    expression <- expression_matrix[expression_matrix$Hugo_Symbol %in% include_genes,]
  }else{
    no <- is.na(expression_matrix$Hugo_Symbol)
    expression <- expression_matrix[!no,]
  }
  dim(expression)
  if(!(is.null(include_samples))){
    #### used to narrow samples for runs with even eas, afr, and euro samples ####
    include_samples <- c("Hugo_Symbol",include_samples)
    final.expression <- expression[,colnames(expression) %in% include_samples]
  }else{
    final.expression <- subset(expression, select =  -Entrez_Gene_Id) #get rid of entrez id
  }
  return(final.expression)
}

## expects matrix with samples as row names 
trim.expr.matrix_OG <- function(expression_matrix, include_samples){
    #### used to narrow samples for runs with even eas, afr, and euro samples ####
    final.expression <- expression_matrix[rownames(expression_matrix) %in% include_samples,]
  return(final.expression)
}

elasticnet.run <- function(in.matrix, directory, out_file, en.mix, seed=4831){
  print(en.mix)
  set.seed(seed)
  
  split <- initial_split(in.matrix, strata = subtype)
  train <- training(split)
  test <- testing(split)
  
  ## -- check that there is no intersection between train and test set -- ##
  length(intersect(rownames(train), rownames(test)))
  
  #validation set 
  val_set <- validation_split(train,
                              strata = subtype,
                              prop = 0.80)
  ### logistic regression ###
  ### set model - we will tune the penalty mixture of elasticnet is passed in ###
  ###ELASTIC NET ####
  lr_mod <-
    logistic_reg(penalty = tune(), mixture = en.mix) %>%
    set_engine("glmnet") #provides variable importance scores

  ### create recipe, remove indicator values that only contain 0 and noramlize the predictors ###
  rec <-
    recipe(subtype ~ ., data = train) 
    #%>%
    # step_zv(all_predictors()) %>%
    # step_normalize(all_predictors())
  
  #### create workflow ###
  lr_wf <- workflow() %>%
    add_model(lr_mod) %>%
    add_recipe(rec)
  
  #grid for tuning 
  lr_reg_grid <- tibble(penalty = 10^seq(-4,-1,length.out = 30))
  
  lr_res <-
    lr_wf %>% 
    tune_grid(val_set,
              grid = lr_reg_grid,
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(yardstick::roc_auc)) #first use roc metric to get tune9) parameter 
  
  ### view top models predicted with tuning ###
  lr_res %>%
    show_best("roc_auc", n = 15) %>%
    arrange(penalty)
  
  ### printing out all tested penalty to get correct slice number ###
  print(lr_res %>% 
          collect_metrics(), n = 30)
  
  ### get the slice with the best roc_auc ###
  best <- lr_res %>% 
    collect_metrics()
  
  ### selecting the penalty with best ROC, if there is a tie one is randomly selected ###
  w.score <- max(best$mean)
  chosen <- best[best$mean == w.score,]
  row <- sample(1:nrow(chosen), 1)
  the.one <- chosen[row,]
  penal <- the.one$penalty
  
  lr_best <-
    lr_res %>% 
    collect_metrics() %>%
    arrange(penalty) %>% 
    dplyr::slice(row)
  
  #------commented out for library issues, only code for confusion matrix------------- After tuning to get confusion matrix ---------------#
  ####TEMPORARILY COMMENTED OUT BECAUSE OF SPEC MASKING ISSUES #####
  # lr_res <-
  #   lr_wf %>%
  #   tune_grid(val_set,
  #             grid = lr_reg_grid,
  #             control = control_grid(save_pred = TRUE), #save model predictions
  #             metrics = yardstick::metric_set(spec)) #run this second round to get specificity
  # print("best")
  # #for spec
  # lr_auc <-
  #   lr_res %>%
  #   collect_predictions(parameters = lr_best)
  # 
  # predicted <- lr_auc$.pred_class
  # actual <- lr_auc$subtype
  # confusion_matrix <- table(predicted,actual) #want this printed out
  # out_matrix <- paste0(directory, "/", out_file,"_confusion_matrix.txt")
  # write.table(confusion_matrix, file = out_matrix, col.names = T)

  
  #========================LAST FIT =======================================
  #last fit 
  last_lr_mod <- 
    logistic_reg(penalty = penal, mixture = en.mix) %>% 
    set_engine("glmnet", importance = "impurity") #provides variable importance scores 
  
  #last workflow
  last_wf <- 
    lr_wf %>%
    update_model(last_lr_mod)
  
  last_fitt <- last_fit(last_wf, split) #split is our train and test data from beginning - it trains on both train set and validation set this time and test on test set  
  pred <- collect_metrics(last_fitt)
  write.table(pred, paste0(directory, "/",out_file,"_metrics.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  desc.scores <- last_fitt %>%
    extract_fit_engine() %>%
    vi()
  desc.scores <- desc.scores[desc.scores$Importance > 0,]
  print(desc.scores)
  ## -- added run information to check genes in ancestry signatures -- ## 
  #out_scores <- paste0(directory, "/", out_file,"_general_info.txt")
  print(paste0("Writing to ",out_file))
  write.table(desc.scores, file = paste0(directory, "/", out_file,"_all_sig.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

  #look at model
  my_model <- extract_workflow(last_fitt)
  dat <- tidy(my_model)
  write.table(dat, paste0(directory, "/",out_file,"_model_coefficients.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  write_rds(my_model, paste0(directory,"/",out_file, "_EN_model.rds"))
  return(my_model)
}

## for each model fit, call this on each ancestry ## 
model.metrics <- function(my_model, dat, directory, out_file, subtype1, subtype2)
{
  pred_prob <- predict(my_model, dat, type = "prob")
  pred_class <- predict(my_model, dat, type = "class")
  results <- dat %>% dplyr::select(subtype) %>% bind_cols(pred_class, pred_prob) 
  results <- rownames_to_column(results, "sample_id")
  write.table(results, paste0(directory,"/", out_file,"_results.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  # new.res1 <- bind_cols(dat$subtype, results$.pred_Basal, results$.pred_Luminal, results$.pred_class)
  # colnames(new.res) <- c("subtype", ".pred_Basal", ".pred_Luminal", ".pred_class")
  new.res <- bind_cols(dat$subtype, results[,4], results[,5], results$.pred_class)#TODO will probably get an error here
  colnames(new.res) <- c("subtype", paste0(".pred_",subtype1), paste0(".pred_",subtype2), ".pred_class")#TODO will probably get an error here
  confusion <- conf_mat(results, truth = subtype,estimate = .pred_class)
  confusion
  conf.df <- as.data.frame(confusion$table)

  write_delim(conf.df, paste0(directory,"/", out_file,"_confusion_matrix.tsv"), delim  = "\t")
  #auc <- roc_auc(new.res, truth = subtype, estimate = .pred_Basal)
  auc <- roc_auc(new.res, truth = subtype, estimate = paste0(".pred_",subtype1)) #TODO will probably get an error here
  #- all can be put in 1 matrix - #
  senss <- yardstick::sens(results, truth = subtype, estimate = .pred_class)
  specc <- yardstick::spec(results, truth = subtype, estimate = .pred_class)
  acc <- accuracy(results, truth = subtype, estimate = .pred_class)
  prec <- yardstick::precision(results, truth = subtype,estimate = .pred_class)
  re <- yardstick::recall(results, truth = subtype,estimate = .pred_class)
  f <- yardstick::f_meas(results, truth = subtype,estimate = .pred_class)
  kapp <- kap(results, truth = subtype,estimate = .pred_class)
  mccc <- mcc(results, truth = subtype,estimate = .pred_class)
  metrics <- rbind(auc, acc, senss, specc, prec, re, f, kapp, mccc)
  print(metrics)
  write.table(metrics,paste0(directory,"/", out_file,"_metrics.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  #return(metrics)
}  
  
  
  

