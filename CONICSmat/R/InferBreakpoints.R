#' Transforms data matrix from log2(CPM/10+1) to CPM
#'
#' This function transforms counts in a data matrix from log2(CPM/10+1) to CPM
#' @param expmat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @keywords Transform
#' @export
#' @examples
#' removeLogScale(suva_exp)

removeLogScale=function(expmat){
  res=((2^expmat)-1)*10
  return(res)
}

#' Visualize normal/tumor ratio of genes expressed from the same chromosome and sorted by position
#'
#' This function allows to visualize the smoothed normal/tumor ratio of genes expressed from the same chromosome and sorted by position. This give information on breakpoints as well as CNV status (gain/loss)
#' @param mat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param normal Indices of columns holding normal cells.
#' @param tumor Indices of columns holding tumor cells.
#' @param windowsize Size of the window smoothing the gene expression data
#' @param gene_pos Matrix with positional information for each gene
#' @param chr Chromosome to plot
#' @param patients Vector of patient identifiers for each cell
#' @param patient Id of the patient to visualize
#' @param breakpoints Chromosomal positions of chromosome arms
#' @export
#' @examples
#' detectBreakPoints (suva_exp,normal,tumor,101,gene_pos,1,patients,patient,breakpoints)

detectBreakPoints= function(mat,normal,tumor,windowsize,gene_pos,chr,patients,patient,breakpoints){
  normal_exp=rowMeans(removeLogScale(mat[,normal]))
  tumor_exp=rowMeans(removeLogScale(mat[,intersect(which(patients==patient),tumor)]))
  genes = which(normal_exp > 5 & tumor_exp > 5)
  normal_exp=log2(normal_exp+1)
  tumor_exp=log2(tumor_exp+1)
  #normal_exp[which(normal_exp>10)]=10
  #tumor_exp[which(tumor_exp>10)]=10
  ratio_nt=tumor_exp[names(genes)]-normal_exp[names(genes)]
  ratio_nt=ratio_nt-median(ratio_nt)
  gp=gene_pos[order(gene_pos[,"start_position"]),]
  gp=gp[which(gp[,3]==chr),]
  rownames(gp)=gp[,2]
  target_genes=intersect(gp[,2],names(genes))
  if (length(target_genes)>windowsize){
    rat=zoo::rollapply(ratio_nt[target_genes],windowsize,mean,1,align="center")
    pos=gp[target_genes,4][(windowsize/2):(length(target_genes)-(windowsize/2))]
    at=which(abs(pos-breakpoints[(chr*2)-1,"End"])==min(abs(pos-breakpoints[(chr*2)-1,"End"])))
    plot(rat,pch=16,ylim=c(-2,2),ylab="Log2 Tumor/Normal ratio",main=paste("Chromosome",chr))
    abline(h=0,lty=16,col="grey",lwd=2)
    abline(h=-1,lty=16,col="blue",lwd=2)
    abline(h=0.58,lty=16,col="red",lwd=2)
    #add line for centromer
    abline(v=at,lty=3)
    text((at/2)-0.5,2,"p")
    text((at+(((length(target_genes)-windowsize)-at)/2)+0.5),2,"q")
  }
  else{
	plot(c(1,1,1,1,1,1,1,1,1,1,1,1,1),pch=16,ylim=c(-2,2),ylab="Log2 Tumor/Normal ratio",main=paste("Chromosome",chr))
	text(6,0,"Too few genes \n on chromosome to \n infer CNV pattern")
  }
}

#' Applies detectBreakPoints() for all chromosomes and writes a pdf file with figures
#'
#' This function applies detectBreakPoints() for all chromosomes and writes a pdf file with figures
#' @param mat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param normal Indices of columns holding normal cells.
#' @param tumor Indices of columns holding tumor cells.
#' @param windowsize Size of the window smoothing the gene expression data
#' @param gene_pos Matrix with positional information for each gene
#' @param patients Vector of patient identifiers for each cell
#' @param patient Id of the patient to visualize
#' @param fname Name for the output file
#' @param breakpoints Chromosomal positions of chromosome arms
#' @export
#' @examples
#' detectBreakPoints (suva_exp,normal,tumor,101,gene_pos,1,patients,patient,breakpoints)

plotAllChromosomes= function (mat,normal,tumor,windowsize,gene_pos,patients,patient,fname,breakpoints){
  pdf(paste(fname,".pdf",sep=""))
  par(mfrow=c(2,2))
  for (i in 1:22){
    chr=i
    detectBreakPoints(mat,normal,tumor,windowsize,gene_pos,chr,patients,patient,breakpoints)
  }
  dev.off()
}
