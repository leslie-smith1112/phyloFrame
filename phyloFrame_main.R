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
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/new_exomeAF_nosex_ordering.R")
source("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/expression_elasticnet.R")


get.network <- function(tissue.network, dz.sig, neighbor, edge.weight){
  neighbor <- 2
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

get.genes.V2 <- function(tiss.net, dz.genes, the.node, the.edge, expression, exomeAF){
  network.genes <- get.network(tiss.net, dz.genes, the.node, the.edge)
  top.anc.dat <- order_frequencies(network.genes, exomeAF, expression) # base genes will be model genes (the base sig)
  top.anc.dat <- unlist(top.anc.dat)
  # make sure genes are in the expression matrix
  top.anc.dat <- top.anc.dat[top.anc.dat %in% colnames(expression)] # make sure genes are in expression matrix
  top.anc.dat <- top.anc.dat[!(top.anc.dat %in% dz.genes)] 
  print(length(top.anc.dat))
  genes <- list("phyloFrame" = top.anc.dat, "benchmark" = NULL) 
  return(genes)
}

## phyloFrame call to run elastic net on the given genes with correct penalty
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



