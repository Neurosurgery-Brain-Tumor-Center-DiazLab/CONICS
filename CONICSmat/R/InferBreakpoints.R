colMeans = function (expmat) {
  n = apply(expmat,2,mean)
  return(n)
}

rowMeans = function (expmat) {
  n = apply(expmat,1,mean)
  return(n)
}

#' Maps stings to colors
#'
#' This function transforms a vector of strings, eg. the patient ID for each cell, to colors
#' @param st A vector of strings.
#' @keywords Transform
#' @export
#' @examples
#' removeLogScale(suva_exp)

string.to.colors = function(st){
  st=as.character(st)
  n=length(unique(st))
  pt=unique(st)
  co=rainbow(n, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10, alpha = 1)
  for (i in 1:n){st=gsub(pt[i],co[i] , st)}
  return(st)
}

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
#' @param patients Optional: Vector of patient identifiers for each cell
#' @param patient Optional: Id of the patient to visualize
#' @param Optional: Breakpoints Chromosomal positions of chromosome arms
#' @export
#' @examples
#' detectBreakPoints (suva_exp,normal,tumor,101,gene_pos,1)

detectBreakPoints= function(mat,normal,tumor,windowsize,gene_pos,chr,patients=NULL,patient=NULL,breakpoints=NULL){
	#This can be refactored to speed up the function. The first part of the calculation is redundant for each function call
  normal_exp=rowMeans(removeLogScale(mat[,normal]))
  if (!is.null(patient)){
	if (length(intersect(which(patients==patient),tumor))>1){
		tumor_exp=rowMeans(removeLogScale(mat[,intersect(which(patients==patient),tumor)]))
	}
	else {
		tumor_exp=removeLogScale(mat[,intersect(which(patients==patient),tumor)])
	}
  }
  else{
	tumor_exp=rowMeans(removeLogScale(mat[,tumor]))
  }
  genes = which(normal_exp > 5 & tumor_exp > 5)
  normal_exp=log2(normal_exp+1)
  tumor_exp=log2(tumor_exp+1)
  ratio_nt=tumor_exp[names(genes)]-normal_exp[names(genes)]
  ratio_nt=ratio_nt-median(ratio_nt)
  gp=gene_pos[order(gene_pos[,"start_position"]),]
  gp=gp[which(gp[,3]==chr),]
  target_genes=intersect(gp[,2],names(genes))
  if (length(target_genes)>windowsize){
    rat=zoo::rollapply(ratio_nt[target_genes],windowsize,mean,1,align="center")
    plot(rat,pch=16,ylim=c(-2,2),ylab="Log2 Tumor/Normal ratio",main=paste("Chromosome",chr))
    abline(h=0,lty=16,col="grey",lwd=2)
    abline(h=-1,lty=16,col="blue",lwd=2)
    abline(h=0.58,lty=16,col="red",lwd=2)
    #add line for centromer
	if (!is.null(breakpoints)){
		rownames(gp)=gp[,2]
		pos=gp[target_genes,4][(windowsize/2):(length(target_genes)-(windowsize/2))]
		at=which(abs(pos-breakpoints[(chr*2)-1,"End"])==min(abs(pos-breakpoints[(chr*2)-1,"End"])))
		abline(v=at,lty=3)
		text((at/2)-0.5,2,"p")
		text((at+(((length(target_genes)-windowsize)-at)/2)+0.5),2,"q")
	}
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
#' @param fname Name for the output file
#' @param Optional: patients Vector of patient identifiers for each cell
#' @param Optional: patient Id of the patient to visualize
#' @param Optional: breakpoints Chromosomal positions of chromosome arms
#' @export
#' @examples
#' detectBreakPoints (suva_exp,normal,tumor,101,gene_pos,1,patients,patient,breakpoints)

plotAllChromosomes= function (mat,normal,tumor,windowsize,gene_pos,fname,patients=NULL,patient=NULL,breakpoints=NULL,offs=40){
  pdf(paste(fname,".pdf",sep=""))
  par(mfrow=c(2,2))
  res=c()
  bps=c(0)
  for (i in 1:22){
    chr=i
    r=detectBreakPoints(mat,normal,tumor,windowsize,gene_pos,chr,patients,patient,breakpoints)
	res=c(res,r)
	bps=c(bps,bps[i]+length(r))
	
	#For future breakpoint detection
	#if (i==1){
		#segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
	#}
  }
  par(mfrow=c(1,1))
  map = squash::makecmap(res, colFn = squash::bluered)
  plot(c(1:length(res)),rep(1,length(res)),col=squash::cmap(res, map = map),pch=15,xaxt='n',yaxt='n',ann=FALSE)
  squash::vkey(map, title="",stretch = 0.8)
  for (i in 1:21){
      bp1=bps[i+1]-offs
      abline(v=bp1,lty=16,col="lightgrey")
  }
  axis(1, at=bps[2:23]-((bps[2:23]-bps[1:22])/2)-offs, labels=1:22, las=2,col="grey")
  dev.off()
}

#' Fast calculation of a distance matrix from a gene expression matrix (R dist function is slow)
#'
#' This function calculates a distance matrix from a gene expression matrix (R dist function is slow)
#' @param mat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @export
#' @examples
#' pdist (suva_expr)

pdist= function(tmat){
  mtm=Matrix::tcrossprod(tmat)
  sq=rowSums(tmat^2)
  out0=outer(sq, sq, "+") - 2 * mtm
  out0[out0 < 0]=0
  return(sqrt(out0))
}

#' Generates a heatmap of gene expression across the genome for each cell 
#'
#' This function generates a heatmap of gene expression across the genome for each cell 
#' @param mat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param normal Vector of indices of columns holding normal cells.
#' @param plotcells Vector of indices of columns holding cells to plot.
#' @param gene_pos Matrix with positional information for each gene
#' @param windowsize Integer size of the window smoothing the gene expression data
#' @param chr Optional: Boolenan: Hierarchically cluster rows (cells) to identify clones. default: FALSE
#' @param expThresh Optional: Double threshold for the average expression of a gene across all cells to be considered in the calculation. 
#' @param thresh Optional: Double visualtization threshold: Rations above/below this threshold (x>t | x<(-t)) will be set to t | -t 
#' @param colo Optional: Color for each cell, e.g. for each patient. Will be displayed next to the heatmap
#' @param colo2 Optional: Second color for each cell, e.g. for each region. Will be displayed next to the heatmap
#' @param retMat Optional: return normalized matrix for custom plotting
#' @export
#' @examples
#' plotChromosomeHeatmap (suva_expr,normal,tumor,gene_pos,chr=c(11))


plotChromosomeHeatmap= function (mat,normal,plotcells,gene_pos,windowsize=121,chr=FALSE,expThresh=0.4,thresh=1,colo=NULL, colo2=NULL, retMat=FALSE,plotdendrogram=FALSE){
	
	#Create average of reference (for all non-tumor cells, even if not plotted)
	ref=rowMeans(mat[,colnames(mat)[normal]])
	#Define matrix to be plotted (remove not plotted normal and tumor cells)
	normal=intersect(normal,plotcells)
	tumor=setdiff(plotcells,normal)
	gexp=mat[,c(normal,tumor)]
	
	#Set color vectors
	if (!is.null(colo)){colo=colo[plotcells]}
	if (!is.null(colo2)){colo2=colo2[plotcells]}
	
	#Filter matrices
	ref=ref[which(ref>expThresh)]
	gexp=gexp[names(ref),]
	gexp=gexp[which(rowMeans(gexp)>expThresh),]
	ref=ref[rownames(gexp)]
	
	#Get sorted gene positions
	#gp=getGenePositions(rownames(gexp))
	gp=gene_pos
	gp=gp[which(gp[,"chromosome_name"] %in% 1:22),]
	gp=gp[order(as.numeric(gp[,"chromosome_name"]),as.numeric(gp[,"start_position"])),]
	
	#Reduce matrix to genes with known positions
	gexp=gexp[intersect(gp[,2],rownames(gexp)),]
	ref=ref[rownames(gexp)]
	
	#Calculate ratios
	rat=gexp-ref
	
	#Smoothed expression
	print("Generating distance Matrix")
	d=apply (rat,2,function (x) zoo::rollapply(x,mean,width=windowsize,align="center"))
	
	#Center in each cell and set to boundaries
	d=apply(d,2,function(x) x-mean(x))
	
	#Set min/max 
	d[which(d>thresh)]=thresh
	d[which(d<(-thresh))]=(-thresh)
	
	#Order by chromosome
	print(paste("Distance matrix ready with dimensions",nrow(d),ncol(d)))
	print("Clustering cells")
	if (chr){
		normal=1:length(normal)
		tumor=(length(normal)+1):ncol(d)
		dmat=as.dist(pdist(t(d[,tumor])))
		hc = hclust(dmat)
		cellOrder = c(normal,tumor[hc$order])
		d=d[,cellOrder]
		if (!is.null(colo)){
			colo=colo[cellOrder]
		}
		if (!is.null(colo2)){
			colo2=colo2[cellOrder]
		}
	}
	print("Cells clustered, starting to plot")
	
	if (plotdendrogram){
		fr=length(tumor)/(length(normal)+length(tumor))
		print(fr)
		mat = t(matrix(c(1,2,3,3), nrow=2, byrow=TRUE))
		layout(mat, widths = c(0.1,0.9), heights = c(fr, (1-fr), 1), FALSE)
		par(mar=c(2*fr,0,4,0))
		plot(as.dendrogram(hc),ann = F,axes = F,horiz=T,leaflab = "none")
		plot(1,1,col="white",axes=FALSE,xlab = "",ylab="")
		par(mar=c(5,0,4,0))
	}
	
	dsize=800/length(plotcells)
	i=0
	apply(d,2,function(x){
	  res=x
	  #Counter increase
	  i <<- i + 1
	  map = squash::makecmap(c(-thresh,thresh), colFn = squash::darkbluered,n=256)
	  #map = squash::makecmap(c(max(-thresh,mind),min(thresh,maxd)), colFn = squash::darkbluered,n=256)
	  if (i==1){
		plot(c(1:length(res)),rep(1,length(res)),col=squash::cmap(res, map = map),cex=dsize,pch=".",xaxt='n',yaxt='n',ann=FALSE,ylim=c(0,ncol(d)))
		if (!is.null(colo)){
			points((length(res)+1):(length(res)+20),rep(1,20),col="white",pch=".",cex=dsize)
			points((length(res)+21):(length(res)+110),rep(1,90),col=colo[i],pch=".",cex=dsize)
		}
		if (!is.null(colo2)){
			points((length(res)+111):(length(res)+160),rep(1,50),col="white",pch=".",cex=dsize)
			points((length(res)+161):(length(res)+250),rep(1,90),col=colo2[i],pch=".",cex=dsize)
		}
		squash::hkey(map, title="",stretch = 0.1)
	  }
	  else{
		points(c(1:length(res)),rep(i,length(res)),col=squash::cmap(res, map = map),pch=".",cex=dsize)
		if (!is.null(colo)){
			points((length(res)+1):(length(res)+20),rep(i,20),col="white",pch=".",cex=dsize)
			points((length(res)+21):(length(res)+110),rep(i,90),col=colo[i],pch=".",cex=dsize)
		}
		if (!is.null(colo2)){
			points((length(res)+111):(length(res)+160),rep(i,50),col="white",pch=".",cex=dsize)
			points((length(res)+161):(length(res)+250),rep(i,90),col=colo2[i],pch=".",cex=dsize)
		}
	  }
	})
	bps=c(0)
	for (i in 1:21){
	  gp=gp[which(gp[,2] %in% rownames(gexp)),]
	  bp1=which(gp[,"chromosome_name"]==i+1)[1]-windowsize/2
	  bps=c(bps,bp1)
	  abline(v=bp1,lty=16)
	}
	axis(1, at=bps[2:23]-((bps[2:23]-bps[1:22])/2), labels=1:22, las=2,col="grey")
	bps=round(bps,0)
	abline(h=length(normal),lty=16)
	if (retMat){
		chrpos=c()
		for (i in 2:22){chrpos=c(chrpos,rep(i-1,bps[i]-bps[i-1]))};chrpos=c(chrpos,rep(22,nrow(d)-bps[22]))
		return(rbind(chrpos,t(d)))
	}
	print("Done")
}
