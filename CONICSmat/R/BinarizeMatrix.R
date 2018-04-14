#' Binarize a vector of posterior probabilities.
#'
#' This function allows to binarize posterior probabilities of a GMM. Input is a vector of posterior probabilities for the component with the larger mean from a Gaussian Mixture Model.
#' @param mixmdl Vector of posterior probabilities between 0 and 1.
#' @param normal Vector of positions indicating the indices that identify posteriors assigned to normal cells.
#' @param tumor Vector of positions indicating the indices that identify posteriors assigned to tumor cells.
#' @param threshold Posterior threshold level. The presence of a CNV is assigned to a cell if its posterior exceeds the threshold.
#' @param withna Should posteriors that can't be assigned to component 1, but also not to component2, be set to 0 or NA (defualt is set to NA).
#' @keywords Binarizee vector
#' @export
#' @examples
#' binarizeCalls(mixmdl,normal,tumor,threshold,withna=T)

binarizeCalls = function (mixmdl,normal,tumor,threshold,withna=T){
  g1=length(which(mixmdl[normal]>threshold))
  g2=length(which((1-mixmdl[normal])>threshold))
  status=""
  if (g1>g2){
    status="del"
    resV=ifelse(mixmdl<(1-threshold),1,0)
    if (withna==T){
      resV[which(mixmdl<threshold & mixmdl>(1-threshold))]=NA
    }
  }
  else{
    status="amp"
    resV=ifelse(mixmdl>threshold,1,0)
    if (withna==T){
      resV[which(mixmdl<threshold & mixmdl>(1-threshold))]=NA
    }
  }
  res <- list("integer" = resV, "status" = status)
  return(res)
}

#' Binarize a matrix of posterior probabilities.
#'
#' This function allows to binarize a matrix of posterior probabilities. Input is a matrix containing the posterior probabilities for the component with the larger mean from a Gaussian Mixture Model across different regions. 
#' @param mixmat Matrix of posterior probabilities between 0 and 1.
#' @param normal Vector of positions indicating the indices that identify posteriors assigned to normal cells.
#' @param tumor Vector of positions indicating the indices that identify posteriors assigned to tumor cells.
#' @param threshold Optional: Posterior threshold level. The presence of a CNV is assigned to a cell if its posterior exceeds the threshold.
#' @param withna Optional: Should posteriors that can't be assigned to component 1, but also not to component2, be set to 0 or NA (defualt is set to NA).
#' @keywords Binarizee vector
#' @export
#' @examples
#' binarizeMatrix (mixmat,normal,tumor)

binarizeMatrix = function (mixmat,normal,tumor,threshold=0.8,withna=T){
  res=apply(mixmat,2,function (x) binarizeCalls(x,normal,tumor,threshold,withna=T)$integer)
  nms=apply(mixmat,2,function (x) binarizeCalls(x,normal,tumor,threshold,withna=T)$status)
  colnames(res)=paste(nms,colnames(res),sep="_")
  return(res)
}

#' Visualize a matrix of binary CNV assignments.
#'
#' This function visualizes a matrix of binary CNV assignment. A 1 indicates the presence, a 0 the absence of a CNV
#' @param mati A cells X regions matrix 
#' @param normal Vector of positions indicating the indices that identify posteriors assigned to normal cells.
#' @param tumor Vector of positions indicating the indices that identify posteriors assigned to tumor cells.
#' @param patients A vector of length(nrow(mati)) indicating the patient for each cell.
#' @param patient Optional: Which patient should the matrix be plotted for.
#' @keywords Binarizee vector
#' @export
#' @examples
#' plotBinaryMat(mati,patients,normal,tumor,patient="MGH96")

plotBinaryMat = function(mati,patients,normal,tumor,patient=NULL,k=3){
	celltypes=rep("Tumor",length(normal)+length(tumor));celltypes[normal]="Normal";names(celltypes)=c(names(normal),names(tumor))
	patientcolors =data.frame(celltypes)
	patientcolors=cbind(patientcolors,patients)
	rownames(patientcolors)=names(celltypes)
	rownames(mati)=names(celltypes)
	if (!is.null(patient)){
		p=pheatmap::pheatmap(t(mati[which(patients==patient),]),cluster_cols=T, cutree_cols = 3,annotation=patientcolors, col=c("lightgrey","black"),border_color = "grey60",show_colnames = F,clustering_distance_cols="euclidean")
	}
	else{
		p=pheatmap::pheatmap(t(mati),cluster_cols=T, cutree_cols = k,annotation=patientcolors, col=c("lightgrey","black"),border_color = "grey60",show_colnames = F,clustering_distance_cols="euclidean")
	}
	ord=unique(cutree(p$tree_col, k = k)[p$tree_col[["order"]]])
	numb=table(cutree(p$tree_col, k = k))[ord]
	n=length(numb)
	grid::grid.text(expression(bold("Cluster ID \n(left to right)")),x=rep(0.92),y=c(n*0.03+0.03),gp=grid::gpar(fontsize=8, col="grey"))
	grid::grid.text(ord,x=rep(0.92,length(numb)),y=seq(n*0.03, 0.03, -0.03),gp=grid::gpar(fontsize=8, col="grey"))
	return(cutree(p$tree_col, k = k))
}

#' Calculate p-values from a mixmdl object.
#'
#' This function generates p-values from a mixmdl object based on the mean and variance of the two components
#' @param mixmdl An mixmdl object returned by the mixtools package.
#' @param normal Vector of positions indicating the indices that identify posteriors assigned to normal cells.
#' @param tumor Vector of positions indicating the indices that identify posteriors assigned to tumor cells.
#' @param threshold Optional: A threshold for binarizing the matrix. Example: For gains all cells with a posterior probability >threshold for the component with the larger mean will be assigned a 1.
#' @keywords Binarize vector
#' @export
#' @examples
#' calcPvalue (mixmdl,normal,tumor)


calcPvalue = function (mixmdl,normal,tumor,threshold=0.8){
  largerMean=ifelse(mixmdl$mu[1]>mixmdl$mu[2],1,2)
  smallerMean=ifelse(mixmdl$mu[1]<mixmdl$mu[2],1,2)
  #Figure out if the normal cells are in the component with the larger or smaller mean
  g1=length(which(mixmdl$posterior[normal,largerMean]>threshold))
  g2=length(which((1-mixmdl$posterior[normal,largerMean])>threshold))
  status=""
  if (g1>g2){
    status="del"
	resV=c()
    for (i in 1:length(tumor)){resV=c(resV,pnorm(mean=mixmdl$mu[largerMean],sd=mixmdl$sigma[largerMean],mixmdl$x[tumor[i]],lower.tail=T))}
  }
  else{
    status="amp"
    resV=c()
    for (i in 1:length(tumor)){resV=c(resV,pnorm(mean=mixmdl$mu[smallerMean],sd=mixmdl$sigma[smallerMean],mixmdl$x[tumor[i]],lower.tail=F))}
  }
  resV=p.adjust(resV,method="BH")
  return(resV)
}

#' Calculate p-values from a mixmdl object.
#'
#' This function generates p-values from a log2(CPM/10+1) scaled gene expression matrix using the calcPvalue() function on all the regions given via the regions parameter
#' @param mat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param regions A matrix of regions to be evaluated 
#' @param normFactor A vector of normalization factors for each cell calculated with calcNormFactors()
#' @param normal Vector of positions indicating the indices that identify posteriors assigned to normal cells.
#' @param tumor Vector of positions indicating the indices that identify posteriors assigned to tumor cells.
#' @param gene_pos A matrix of positions for every gene is the expression matrix, calculated via getGenePositions()
#' @param threshold Optional: A threshold for binarizing the matrix. Example: For gains all cells with a posterior probability >threshold for the component with the larger mean will be assigned a 1.
#' @keywords Binarizee vector
#' @export
#' @examples
#' generatePvalMat(suva_expr,regions,normfactor,normal,tumor,gene_pos)

generatePvalMat = function (mat,regions,normfactor,normal,tumor,gene_pos,threshold=0.8){
	res=c()
	for (i in 1:nrow(regions)){
		one=plotChrEnichment (mat,regions[i,1],normFactor,gene_pos,nrow(regions),normal,tumor,regions[i,2],regions[i,3],vis=F)
		if (!is.null(one)){
			pv=calcPvalue(one,normal,tumor)
			res=cbind(res,pv);colnames(res)[i]=rownames(regions)[i]
		}
	}
	rownames(res)=colnames(mat[,tumor])
	return(res)
}

