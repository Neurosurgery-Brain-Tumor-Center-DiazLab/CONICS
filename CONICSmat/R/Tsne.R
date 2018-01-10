#' Identify the most variable genes in a matrix
#'
#' This function identifies the most variable genes from a matrix of raw gene counts to log2(CPM/10+1).
#' @param expmat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param ngenes Return the ngenes most variable genes.
#' @keywords Detect varibale genes
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

detectVarGenes = function(mat,ngenes){
  var.fit <- scran::trendVar(mat, parametric=TRUE, span=0.2)
  var.out <- scran::decomposeVar(mat, var.fit)
  hvg.out <- var.out[which(var.out$FDR <= 0.00001),]
  hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
  genes=rownames(hvg.out)[1:500]
  return(genes)
}

#' Calculate a 2D t-Distributed Stochastic Neighbor Embedding from a gene expression matrix
#'
#' This function calculates a 2D t-Distributed Stochastic Neighbor Embedding from a gene expression matrix in the space of a given set of genes
#' @param expmat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param genes A list of genes used for tSNE generation (HUGO format)
#' @keywords Plot tSNE
#' @export
#' @examples
#' calculateTsne(suva_exp,c("EGFR","PDGFRA","MET","SOX2"))

calculateTsne = function(mat,genes,dims=20,perplexity = 30){
  set.seed(10)
  rt=Rtsne::Rtsne(t(mat[genes,]),initial_dims=dims,perplexity=perplexity)
  return(rt)
}

#' Visualize the expression of a set of genes in tSNE space
#'
#' This function plots the 1st and 2nd tSNE coordinates of cells and visualized the expression of a given set of genes
#' @param rt An object returned by calculateTsne() holding tSNE coordinates
#' @param mat A log2(CPM+1/10) scaled expression matrix
#' @param genes A vector of genes to visualize
#' @keywords Plot tSNE
#' @export
#' @examples
#' plotTsneGene(rt,suva_exp,c("EGFR","PDGFRA","MET","SOX2"))

plotTsneGene  = function(rt,mat,genes,ncols=2){
  par(mfrow=c(ceiling((length(genes)/ncols)),ncols))
  par(mar=c(4,4,4,3))
  for (i in 1:length(genes)){
    map <- squash::makecmap(mat[genes[i],], colFn = squash::bluered)
    plot(rt$Y,pch=16,col=squash::cmap(mat[genes[i],], map = map),xlab="tSNE 1",ylab="tSNE 2",main=genes[i])
    squash::vkey(map, title="",stretch = 0.8)
  }
}


#' Visualize the posterior probabilities for one component of a gaussian mixture model for each cell in tSNE space
#'
#' This function plots the 1st and 2nd tSNE coordinates of cells and visualizes the posterior probabilities for one component of a gaussian mixture model for each cell
#' @param rt An object returned by calculateTsne() holding tSNE coordinates
#' @param mat A log2(CPM+1/10) scaled expression matrix
#' @param genes A vector of posterior probablities to visualize
#' @keywords Plot tSNE
#' @export
#' @examples
#' plotTsneProbabilities(rt,suva_exp,p1)

plotTsneProbabilities  = function(rt,mat,prob,pname,ncols=2){
  par(mfrow=c(1,1))
  par(mar=c(4,4,4,3))
  map <- squash::makecmap(prob, colFn = squash::bluered)
  plot(rt$Y,pch=16,col=squash::cmap(prob, map = map),xlab="tSNE 1",ylab="tSNE 2",main=pname)
  squash::vkey(map, title="",stretch = 0.8)
}

