#' Calculate the BIC (information-theoretic criteria) for a mixture model. 
#'
#' This function calculate the BIC (information-theoretic criteria) for a mixture model. 
#' @param l1 Log likelihood 1
#' @param l2 Log likelihood 1
#' @param n Number of total tests for Benfferroni correction
#' @keywords Detect varibale genes
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

BIC.mix = function(out){
	bc=-2*out$loglik+log(length(out$x))*(length(unlist(out[2:4]))-1)
	return(bc)
}

#' Perform likelihood ratio test given two log likelihoods and a integer numbe rindicating the number of tests performed for multiple testing correction 
#'
#' This function performs a likelihood ratio test given two log likelihoods and a integer numbe rindicating the number of tests performed for multiple testing correction 
#' @param l1 Log likelihood 1
#' @param l2 Log likelihood 1
#' @param n Number of total tests for Benfferroni correction
#' @keywords Detect varibale genes
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

likelihoodRatioTest = function(l1,l2,n){
	lrt= -2*l1 -(-2*l2)
	res=min(pchisq(lrt, df=2, lower.tail=FALSE)*n,1)
	print(n)
	return(res)
}

#' Estimate a two-component Gaussian Mixture Model based on a log2(CPM/10+1) expression matrix 
#'
#' This function estimates a two component Gaussian Mixture Model on log2(CPM/10+1) expression matrix.
#' @param expmat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param chr Number of the chromosome to be considered.
#' @param normFactor A vector of normalization factors for each cell calculated with calcNormFactors()
#' @param groups1 A vector of indices for normal cells
#' @param groups2 A vector of indices for tumor cells
#' @param gene_positions A matrix of positions for every gene in the expression matrix, calculated via getGenePositions()
#' @param start Optional: The start position of the region to be considered
#' @param end Optional: The end position of the region to be considered
#' @keywords Detect varibale genes
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

plotChrEnichment = function(expmat,chr,normFactor,gene_positions,n,groups1=NULL,groups2=NULL,start=NULL,end=NULL,k=2,vis=T){
  par(mfrow=c(2,2))
  if (!is.null(groups1)){
	cellcolor=rep("black",(length(groups1)+length(groups2)));cellcolor[groups2]="red"
  }
  else{
	cellcolor=NULL
  }
  if (!is.null(start)){
    chr_genes=gene_positions[which(gene_positions[,3]==chr & gene_positions[,4]>start &  gene_positions[,5]<end),2]
  }
  else{
    chr_genes=gene_positions[which(gene_positions[,3]==chr),2]
  }
  if(length(chr_genes)>100){
    chr_exp=scale(colMeans(expmat[intersect(chr_genes,row.names(expmat)),])-normFactor)
    wait = chr_exp
    mixmdl = mixtools::normalmixEM(wait,k=k,maxit = 1000,maxrestarts=10)
	out1 = list(x=wait,mu=mean(wait),sigma=sd(wait),lambda=1,loglik=sum(dnorm(wait,mean(wait),sd(wait),log=TRUE)))
	bics = c(BIC.mix(out1),BIC.mix(mixmdl))
	lrt= round(likelihoodRatioTest (out1$loglik,mixmdl$loglik,n),6)
	mixmdl$BIC=bics
	mixmdl$lrt=lrt
	if (mixmdl$lambda [1]<0.02 | mixmdl$lambda [2]<0.02){
			lrt="Not calculated. More than 98 perc of\n cells assigned to one component."
		}
	#print(lrt)
	if (vis==T){
		plot(mixmdl,which=2,breaks=50,col1=c("red","green"),main2=paste("Chr: ",chr,":",start,":",end,"\n","Log likelihood ",round(mixmdl$loglik,1),sep=""),lwd2=3,xlab2="Expression z-score")
	}
    if (length(cellcolor)>1){
      g1=length(which(mixmdl$posterior[groups1,1]>0.95))/length(groups1)*100
      g2=length(which(mixmdl$posterior[groups1,2]>0.95))/length(groups1)*100
      g3=length(which(mixmdl$posterior[groups1,2]<0.95 & mixmdl$posterior[groups1,1]<0.95))/length(groups1)*100
      g4=length(which(mixmdl$posterior[groups2,1]>0.95))/length(groups2)*100
      g5=length(which(mixmdl$posterior[groups2,2]>0.95))/length(groups2)*100
      g6=length(which(mixmdl$posterior[groups2,2]<0.95 & mixmdl$posterior[groups2,1]<0.95))/(length(groups2)+length(groups1))*100
	  if (vis==T){
		barplot(rbind(c(g1,g2,g3),c(g4,g5,g6)),ylim=c(0,100),beside=T,ylab="Percentage of cells",names=c("Cluster","Cluster","Ambigu"),legend = c("Non-malignant", "Malignant"),args.legend = list(title = "Pred. via transcript.", x = "topright", cex = .65),xlab="Predicted via transcriptomics")
		axis(1, at=c(0.5,1,2,3,3.3), line=2, tick=T, labels=rep("",5), lwd=3, lwd.ticks=0,col="red")
		axis(1, at=c(3.5,4,5,6,6.5), line=2, tick=T, labels=rep("",5), lwd=3, lwd.ticks=0,col="green")
		barplot(bics,names=c("1","2"),ylab="BIC",pch=16,xlab="Number of components",log="y")
		plot( runif(length(chr_exp), 0,100),chr_exp,pch=16,col=cellcolor,ylab="Expression z-score",ylim=c(min(chr_exp),(max(chr_exp)+2)),xlab="Cells")
		legend("topright", col=c("black","red"), c("Non-malignant","Malignant"), bty="o",  box.col="darkgreen", cex=.65,pch=16,title="Pred. via transcript.")
	  }
    }
    else{
	  if (vis==T){
		plot( runif(length(chr_exp), 0,100),chr_exp,pch=16,ylab="Expression z-score",ylim=c(min(chr_exp),(max(chr_exp)+2)),xlab="Cells")
		barplot(bics,names=c("1","2"),ylab="BIC",pch=16,xlab="Number of components",log="y")
		hist(mixmdl$posterior[,1],main="Posterior probablility distribution\n component 1",xlab="Posterior probability",breaks=20,xlim=c(0,1))
	  }
    }
    return(mixmdl)
  }
}

#' Estimate and visualize a two-component Gaussian Mixture Model based on a log2(CPM/10+1) expression matrix 
#'
#' This function estimates a two component Gaussian Mixture Model on log2(CPM/10+1) expression matrix and visualizes the results.
#' @param mat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param normFactor A vector of normalization factors for each cell calculated with calcNormFactors()
#' @param regions A matrix of regions to be evaluated 
#' @param gene_pos A matrix of positions for every gene is the expression matrix, calculated via getGenePositions()
#' @param fname A string indicating the name of the output pdf
#' @param normal Optional: A vector of indices for normal cells
#' @param tumor Optional: A vector of indices for tumor cells
#' @keywords Detect varibale genes
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

plotAll = function (mat,normFactor,regions,gene_pos,fname,normal=NULL,tumor=NULL){
  pdf(paste(fname,"_CNVs.pdf",sep=""))
  loglik=c()
  bic=c()
  lrt=c()
  restarts=c()
  for(i in 1:nrow(regions)){
    if (i==1){
      mixmdl=plotChrEnichment(mat,regions[i,1],normFactor,gene_pos,nrow(regions),normal,tumor,regions[i,2],regions[i,3])
      loglik=c(loglik,mixmdl$loglik)
	  bic=c(bic,mixmdl$BIC)
	  lrt=c(lrt,mixmdl$lrt)
	  restarts=c(restarts,mixmdl$restarts)
      if (mixmdl$mu[1]>mixmdl$mu[2]){
        l=mixmdl$posterior[,1]
      }
      else{
        l=mixmdl$posterior[,2]
      }
    }
    else{
      mixmdl=plotChrEnichment(mat,regions[i,1],normFactor,gene_pos,nrow(regions),normal,tumor,regions[i,2],regions[i,3])
      if (!is.null(mixmdl)){
        loglik=c(loglik,mixmdl$loglik)
		bic=c(bic,mixmdl$BIC)
		lrt=c(lrt,mixmdl$lrt)
		restarts=c(restarts,mixmdl$restarts)
        names(loglik)[1]=rownames(regions)[1]
        names(loglik)[length(loglik)]=rownames(regions)[i]
		names(bic)[1]=paste(rownames(regions)[1],"_1_comp",sep="_")
		names(bic)[2]=paste(rownames(regions)[1],"_2_comp",sep="_")
        names(bic)[length(bic)]=paste(rownames(regions)[i],"1_comp",sep="_")
		names(bic)[length(bic)-1]=paste(rownames(regions)[i],"2_comp",sep="_")
        if (mixmdl$mu[1]>mixmdl$mu[2]){
          r=mixmdl$posterior[,1]
        }
        else{
          r=mixmdl$posterior[,2]
        }
        l=cbind(l,r)
        colnames(l)[1]=rownames(regions)[1]
        colnames(l)[ncol(l)]=rownames(regions)[i]
      }
    }
  }
  par(mfrow=c(1,1))
  barplot(sort(loglik),names=names(sort(loglik)),cex.axis=0.8,cex.names=0.7,las=2,ylab="log-likelihood")
  dev.off()
  bicLRmat=matrix(ncol=4,nrow=length(loglik))
  #gmmConvergece=ifelse(restarts==1000,"N","Y")
  bicLRmat[,1]=bic[seq(1,(length(bic)-1),2)]
  bicLRmat[,2]=bic[seq(2,length(bic),2)]
  bicLRmat[,3]=bicLRmat[,1]-bicLRmat[,2]
  #bicLRmat[,4]=restarts
  #bicLRmat[,5]=gmmConvergece
  bicLRmat[,ncol(bicLRmat)]=lrt
  #colnames(bicLRmat)=c("BIC 1 component","BIC 2 components","BIC difference","GMM iterations","GMM converged","LRT adj. p-val")
  colnames(bicLRmat)=c("BIC 1 component","BIC 2 components","BIC difference","LRT adj. p-val")
  rownames(bicLRmat)=names(loglik)
  write.table(bicLRmat,paste(fname,"BIC_LR.txt",sep="_"),sep="\t")
  return(l)
}

color.map <- function(mol.biol) { 
  if (mol.biol=="malignant") "#FF0000" else "#0000FF" 
}


#' Visualize a heatmap of posterior probabilities of cells for component 1 of Gaussian Mixture Models fit for each user-defined region
#'
#' This function allows to plot a heatmap of posterior probabilities of cells for component 1 of Gaussian Mixture Models fit for each user-defined region.
#' @param pmat Matrix of posterior probabilites, as returned by plotAll().
#' @param expmat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param clusters The number of clusters to cut the dendrogram into
#' @param zscoreThreshold A threshold on the row-wise z-scores for visualiztion. Default: 4
#' @param gene_pos A matrix of positions for every gene is the expression matrix, calculated via getGenePositions()
#' @param patients Optional: A vector naming the patient for each cell
#' @param celltypes Optional: A vector naming a putative celltype for each cell
#' @keywords Heatmap of posteriors
#' @export
#' @examples
#' plotHistogram(pmat,suva_expr,2,zscoreThreshold=4)

plotHistogram = function (pmat,expmat,clusters,zscoreThreshold=4,patients=NULL,celltypes=NULL){
  t=zscoreThreshold
  pmat=scale(pmat)
  if (max(pmat)>t){
    pmat[which(pmat>t)]=t
    pmat[which(pmat<(-t))]=(-t)
  }
  else {
    mx=min(max(pmat),abs(min(pmat)))
    sc=t/mx
    pmat=pmat*sc
    pmat[which(pmat>t)]=t
    pmat[which(pmat<(-t))]=(-t)
  }
  if(is.null(patients) & is.null(celltypes)){
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,show_colnames = F,clustering_distance_cols="euclidean")
  }
  else if (is.null(celltypes)){
    patientcolors =data.frame(patients)
    rownames(patientcolors)=colnames(expmat)
    rownames(pmat)=colnames(expmat)
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F,silent=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,annotation=patientcolors,show_colnames = F,clustering_distance_cols="euclidean")
  }
  else if(is.null(patients)){
    patientcolors =data.frame(celltypes)
    rownames(patientcolors)=colnames(expmat)
    rownames(pmat)=colnames(expmat)
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,annotation=patientcolors,show_colnames = F,clustering_distance_cols="euclidean")
  }
  else {
    patientcolors =data.frame(celltypes)
    patientcolors=cbind(patientcolors,patients)
    rownames(patientcolors)=names(celltypes)
    rownames(pmat)=names(celltypes)
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,annotation=patientcolors,show_colnames = F,clustering_distance_cols="euclidean")
  }
print(p)
 return(heat.clusts = cutree(p$tree_col, k = clusters))
}

