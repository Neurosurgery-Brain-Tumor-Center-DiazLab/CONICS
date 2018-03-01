colMeans = function (expmat) {
  n = apply(expmat,2,mean)
  return(n)
}

rowMeans = function (expmat) {
  n = apply(expmat,1,mean)
  return(n)
}


#' Filter and expression matrix
#'
#' This function filters an expression matrix by removing genes which are expressed at low levels and in only very few samples.
#' @param mat1 A genes X samples expression matrix.
#' @param genelist Genes not in this list will be removed from the matrix.
#' @param minExp Only genes with rowSums > minExp are being kept. Default:10
#' @param minCells Only genes expressed in > minCells samples will be kept
#' @keywords Filter matrix
#' @export
#' @examples
#' filterMatrix(suva_exp, chr_genes)

filterMatrix=function(mat1,genelist,minExp=10,minCells=5){
  mat1=mat1[intersect(genelist,rownames(mat1)),]
  mat1=mat1[which(rowSums(mat1)>minExp),]
  percentgenes=apply(mat1,1,function(x){sum(x>0)})
  mat1=mat1[percentgenes>minCells,]
  return(mat1)
}

#' Calculate a normalization factor for each column in an expression matrix
#'
#' This function calculates an adjustment factor for each sample in a log2(CPM/10+1)-scaled expression matrix.
#' @param expmat A genes X samples expression matrix scaled to Log2(CPM/10+1).
#' @keywords Adjustment factor
#' @export
#' @examples
#' calcNormFactors(suva_exp)

calcNormFactors = function (expmat){
  n=colMeans(expmat)
  return(n)
}

#' Normalize a expression matrix of raw counts to log2 (CPM/10+1) values
#'
#' This function normlizes a matrix of raw gene counts to log2(CPM/10+1).
#' @param expmat A genes X samples expression matrix of raw (single cell) RNA-seq counts.
#' @keywords Normalize
#' @export
#' @examples
#' normMat(suva_exp)

normMat = function (expmat){
  sdepth=colSums(expmat)
  rmat=t( t(expmat) / sdepth*1000000 )
  rmat=log2(rmat/10+1)
  return(rmat)
}
