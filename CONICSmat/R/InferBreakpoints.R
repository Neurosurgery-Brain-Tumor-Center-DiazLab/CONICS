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
	rat=c(1,1,1,1,1,1,1,1,1,1,1,1,1)
  }
  return(rat)
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

plotAllChromosomes= function (mat,normal,tumor,windowsize,gene_pos,patients,patient,fname,breakpoints,offs=40){
  pdf(paste(fname,".pdf",sep=""))
  par(mfrow=c(2,2))
  res=c()
  bps=c(0)
  for (i in 1:22){
    chr=i
    r=detectBreakPoints(mat,normal,tumor,windowsize,gene_pos,chr,patients,patient,breakpoints)
	res=c(res,r)
	bps=c(bps,bps[i]+length(r))
	if (i==1){
		#segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
	}
  }
  par(mfrow=c(1,1))
  map <- squash::makecmap(res, colFn = squash::bluered)
  plot(c(1:length(res)),rep(1,length(res)),col=squash::cmap(res, map = map),pch=15,xaxt='n',yaxt='n',ann=FALSE)
  squash::vkey(map, title="",stretch = 0.8)
  for (i in 1:21){
      bp1=bps[i+1]-offs
      abline(v=bp1,lty=16,col="lightgrey")
  }
  axis(1, at=bps[2:23]-((bps[2:23]-bps[1:22])/2)-offs, labels=1:22, las=2,col="grey")
  dev.off()
}

scaleTo=function(mat,thresh){
	res=apply(mat,2,function(x) {})
}
#' Generates a heatmap of gene expression across the genome for each cell 
#'
#' This function generates a heatmap of gene expression across the genome for each cell 
#' @param mat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param normal Indices of columns holding normal cells.
#' @param tumor Indices of columns holding tumor cells.
#' @param windowsize Size of the window smoothing the gene expression data
#' @param chr Optional: Order the heatmap based on the expression from this chromosome
#' @export
#' @examples
#' plotChromosomeHeatmap (suva_expr,normal,tumor,121,chr=11)

plotChromosomeHeatmap= function (mat,normal,tumor,windowsize,chr=NULL,expThresh=0.4,thresh=1){
	ref=rowMeans(mat[,colnames(mat)[normal]])
	gexp=mat[,tumor]
	gexp=gexp[which(ref>expThresh),]
	ref=ref[which(ref>expThresh)]
	gexp=gexp[which(rowMeans(gexp)>expThresh),]

	gp=getGenePositions(rownames(gexp))
	gp=gp[order(as.numeric(gp[,"chromosome_name"]),as.numeric(gp[,"start_position"])),]
	gp=gp[which(gp[,"chromosome_name"] %in% 1:22),]
	gexp=gexp[intersect(gp[,2],rownames(gexp)),]
	ref=ref[rownames(gexp)]

	ref=ref-mean(ref)
	gexp=gexp-colMeans(gexp)
	rat=gexp-ref

	d=apply (rat,2,function (x) zoo::rollapply(x,mean,width=windowsize,align="center"))
	chromo=chr
	# Improve this! Could be without renaming matrix rows
	rownames(gp)=gp[,"hgnc_symbol"]
	if (!is.null(chr)){
		hc = hclust(dist(t(rat[which(gp[rownames(d),"chromosome_name"]==chromo),])))
	}
	else{
		hc = hclust(dist(t(rat)))
	}
	cellOrder <- hc$order
	d=d[,cellOrder]
	
	for (i in 1:ncol(d)){
	  res=d[,i]-mean(d[,i])
	  maxd=max(d)
	  mind=min(d)
	  res[which(res>thresh)]=thresh
	  res[which(res<(-thresh))]=(-thresh)
	  #res=d[,i]
	  map <- squash::makecmap(c(max(-thresh,mind),min(thresh,maxd)), colFn = squash::darkbluered,n=256)
	  if (i==1){
		plot(c(1:length(res)),rep(1,length(res)),col=squash::cmap(res, map = map),pch=15,xaxt='n',yaxt='n',ann=FALSE,ylim=c(0,ncol(d)))
	  }
	  else{
		points(c(1:length(res)),rep(i,length(res)),col=squash::cmap(res, map = map),pch=15,xaxt='n',yaxt='n',ann=FALSE)
		if (i==ncol(d)){
			squash::vkey(map, title="")
		}
	  }
	}
	bps=c(0)
	for (i in 1:21){
	  bp1=which(gp[rownames(d),"chromosome_name"]==i+1)[1]-40
	  bps=c(bps,bp1)
	  abline(v=bp1,lty=16)
	}

	dim(d)
	axis(1, at=bps[2:23]-((bps[2:23]-bps[1:22])/2)-40, labels=1:22, las=2,col="grey")
	
}
