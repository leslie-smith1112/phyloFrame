require(readr)  # for read_csv()
require(dplyr) #data frame handling 
library(tidyr)
library(stringi) #string mnutations 
library(tidymodels)
library(workflows) #tidy models package for bundling model specs 
library(parsnip) #for modeling 
library(workflows) #put model in workflow 
library(rpart.plot)
library(vip)
library(igraph)
#phyloFrame(network, exomeAF, dz.signature, consensus_genes, expression, samples, clinical, out, "breast")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/new_exomeAF_nosex_ordering.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/expression_elasticnet.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/diseases/breast.R")


get.network <- function(tissue.network, dz.sig, neighbor, edge.weight){
  neighbor <- 1
  cut.graph <- tissue.network[((tissue.network$Connection >= 0.2) & (tissue.network$Connection < 0.51)),]
  new.sig <- dz.sig[(dz.sig %in% cut.graph$Gene1 | dz.sig %in% cut.graph$Gene2)]
  ## -- some of the dz genes are not in network -> trim them out -- ##
  dat.graph <- graph.data.frame(cut.graph, directed =  FALSE)
  ## couldn't just pass a list of nodes into the graph for some reason, so converted them into vertices.
  the.vert <- V(dat.graph)$name
  the.vert <- the.vert[the.vert %in% new.sig]
  mode.in <- "all" ## -- param for igraph -- ##
  neighborhood <- ego(dat.graph, order = neighbor, nodes = the.vert, mode = mode.in)
  print(neighborhood)
  nodes <- unlist(neighborhood)
  node.names <- unique(names(nodes))# this is the gene list you need from the neighborhood.
  node.names <- unique(c(node.names, dz.sig))
  print("Completed finding gene signature neighborhood, moving to Exome allele frequencies.")
  return(node.names)
}
get.genes.V1 <- function(tiss.net, dz.genes, the.node, the.edge, expression, exomeAF){
  network.genes <- get.network(tiss.net, dz.genes, the.node, the.edge)
  top.anc.dat <- order_frequencies(network.genes, exomeAF, expression)
  temp.anc <- unlist(top.anc.dat)
  ##only keep top n most variable acnestry genes 
  #temp.anc <- top.anc.dat[!(top.anc.dat %in% dz.genes)]
  temp.anc <- temp.anc[temp.anc %in% colnames(expression)] # make sure genes are in the expression matrix
  num.genes <- length(temp.anc)
  top.anc.dat <- unique(c(dz.genes, temp.anc))
  top.anc.dat <- c(top.anc.dat, "subtype")
  
  expr.variance <- apply(expression, 2, var)
  var.ordered <- expr.variance[order(expr.variance, decreasing = TRUE)]
  genes.keep <- var.ordered[1:num.genes]
  n.variable.genes<- names(genes.keep)
  n.variable.genes <- c(dz.genes, n.variable.genes, "subtype")
  genes <- list("phyloFrame" = top.anc.dat, "benchmark" = n.variable.genes)
  return(genes)
 
}

get.genes.V2 <- function(tiss.net, dz.genes, the.node, the.edge, expression, exomeAF){
  network.genes <- get.network(tiss.net, dz.genes, the.node, the.edge)
  top.anc.dat <- order_frequencies(network.genes, exomeAF, expression) # base genes will be model genes (the base sig)
  top.anc.dat <- unlist(top.anc.dat)
  # print(length(top.anc.dat))
  # if(length(top.anc.dat) == 0){
  #   top.anc.dat = pf_top_varying_genes(expression, 2000)
  # }
  # make sure genes are in the expression matrix
  top.anc.dat <- top.anc.dat[top.anc.dat %in% colnames(expression)] # make sure genes are in expression matrix
  top.anc.dat <- top.anc.dat[!(top.anc.dat %in% dz.genes)] 
  # do not increase expression of initial base genes
  print(length(top.anc.dat))
  genes <- list("phyloFrame" = top.anc.dat, "benchmark" = NULL) 
  return(genes)
}

ancestry.rescale<- function(expr.dat, ancestry.list){
  #do.call()
  #gene list of variable vs equitable do ttest two groups we are prediciting to see correlation  - double check we are actualyl gertting better expression. 
  for (i in 1:length(ancestry.list)) {
    gene.vector <- expr.dat[,ancestry.list[i]]
    cur.max <- max(gene.vector)
    cur.min <- min(gene.vector)
    new.max <- cur.max * 500
    rescaled.vector <- rescale(gene.vector, to = c(cur.min, new.max), from = range(gene.vector, na.rm = TRUE, finite = TRUE))
    expr.dat[,ancestry.list[i]] <- rescaled.vector
  }
  # index <- which(expr.dat== max(expr.dat), arr.ind = TRUE)
  # max.row <- index[1,1]
  # max.col <- index[1,2]
  # cur.max <- expr.dat[max.row, max.col]
  # new.max <- 2*cur.max
  # temp.expr <- expr.dat[,ancestry.genes]
  # lapply(temp[,ancestry.genes], gene.rescale)
  return(expr.dat)
}

## phyloFrame 2.0 added to accomodate for selected genes being null for inital elastic net phyloframe run 
phyloFrame <- function(selected.genes, expression.mat, directory, out.file, mixture){
  
  if(is.null(selected.genes) == FALSE){
    expr <- expression.mat[,colnames(expression.mat) %in% selected.genes]
  }else{
    expr <- expression.mat
  }
  #expr[1:5,1:5]
  expr$subtype <- as.factor(expr$subtype)
  the.model <- elasticnet.run(in.matrix = expr, directory = directory, out_file = out.file, mixture, 4831)
  return(the.model)
}



