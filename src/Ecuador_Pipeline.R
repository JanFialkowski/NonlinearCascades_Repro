library(igraph)
library(data.table)
library(glue)

ReadRawGraph<-function(fileloc){
  frame <- read.table(fileloc,sep=",",header=TRUE)
  Graph <- graph_from_data_frame(frame)
  return(Graph)
}

ProcessRawGraph <- function(Graph){
  E(Graph)$weight <- E(Graph)$total_tax
  Graph <- t(simplify(Graph)) # t, because the direction in the file is reversed simplify removes self-loops and adds multi-edges
  return(Graph)
}

ReadIDData <- function(fileloc){
  output <- fread(fileloc)
  if("id" %in% names(output)){
    setkey(output,id)
  }
  return(output)
}

DataPrepESRI <- function(Graph,IDData,run_id, Targetfile=NULL, psi_mat=FALSE, track_h = FALSE){
  library(GLcascade)
  library(data.table)
  Dubbidu <- list()
  nace_conv_mat <- as.data.table(GLcascade::nace_conv_mat)
  
  Dubbidu$W <- as_adjacency_matrix(Graph, attr="weight")
  ISIC <- substring(IDData[.(as.numeric(V(Graph)$name)),on="id"]$ciiu_4n4,2,5)
  Dubbidu$p <- sapply(ISIC,function(x){nace_conv_mat[Reference.to.ISIC.Rev..4==x,nace4_num][1]})
  Dubbidu$p[names(Dubbidu$p)=="0513"]<-"0510"
  Dubbidu$p[is.na(Dubbidu$p)]<-"9999"
  Dubbidu$p <- as.numeric(Dubbidu$p)
  Dubbidu$p_market <- as.numeric(Dubbidu$p)
  
  ess_mat_n4_ihs <- GLcascade::ess_mat_n4_ihs
  colnames(ess_mat_n4_ihs) <- as.character(as.numeric(colnames(ess_mat_n4_ihs)))
  rownames(ess_mat_n4_ihs) <- as.character(as.numeric(rownames(ess_mat_n4_ihs)))
  firm_nace4_affiliation_subset <- as.character(sort(unique(as.numeric(Dubbidu$p))))
  Dubbidu$ess_mat_sec <- ess_mat_n4_ihs[firm_nace4_affiliation_subset , firm_nace4_affiliation_subset]
  
  Dubbidu$psi_mat <- psi_mat
  Dubbidu$track_h <- track_h
  Dubbidu$use_rcpp <- TRUE
  Dubbidu$run_id <- run_id
  if(!is.null(Targetfile)){saveRDS(Dubbidu,Targetfile)}
  return(Dubbidu)
}