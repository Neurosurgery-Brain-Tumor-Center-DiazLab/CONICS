library(scde)
library(PCIT)
library(boot)
library(parallel)
library(raster)
library(flashClust)

args=commandArgs(TRUE)

#smallest number representable by the machine
adj<-.Machine$double.eps

center_gene=args[1]
#Gene expression matrix, genes X cells
expression=args[2]
#Metadata File
#Number of cores for computation
nCores<-as.numeric(args[3])
#Starting threshold for neighbour correlationof gene of interest (is decreased until at least x neighbours are found)
cor_threshold=as.numeric(args[4])
# How many direct neighbours of gene of interest should be searched for
min_neighbours=as.numeric(args[5])

#How many raw reads should at least map to a gene to consider the gene
minRawReads=as.numeric(args[6])
#Percentage (0.15 =15%) of cells expressing a gene to consider the gene
percentCellsExpressing=as.numeric(args[7])
#How many genes should be expressed at least in a cell to consider the cell 
minGenesExpr=as.numeric(args[8])
#How deep should the fav gene analysis search
depth=as.numeric(args[9])
dir = args[10]
base = args[11]

##Start of calculations
#Read matrix
a=read.table(expression,sep="\t",row.names=1)

setwd(dir)


#Remove genes with low raw counts and ribosomal protein endcoding RNA (andoptional) 
a=a[rowSums(a)>minRawReads,]
#a=a[!stri_startswith_fixed(row.names(a),"RP"),]

#filter cells based on expression of favorite gene
a=a[,which(a[center_gene,]>5)]

#Filter for cells with at least x genes expressed
nexp=apply(a,2,function(x){sum(x>0)})
a=a[,which(nexp>minGenesExpr)]

#Filter for genes expressed in at least y cells
percentgenes=apply(a,1,function(x){sum(x>0)}/ncol(a))
a=a[percentgenes>percentCellsExpressing,]

#How many cells and fgenes go into analysis?
print(dim(a))

##fit model with SCDE
modelFit<-scde.error.models(counts=a,groups=NULL,zero.lambda=0.1,n.cores=nCores,threshold.segmentation=T,save.crossfit.plots=F, save.model.plots=F,verbose=1)
exprEstModel<-exp(scde.expression.magnitude(models=modelFit,counts=a))
failureProbs<-scde.failure.probability(models=modelFit,counts=a)


#calculate SCDE adjusted correlations
nonFailureProbs<-1-failureProbs
logData<-as.matrix(log2(exprEstModel+1))
nGenes<-nrow(exprEstModel)
genes<-rownames(exprEstModel)
weightMat<-matrix(NA,nGenes,nGenes)
diag(weightMat)<-3
rownames(weightMat)<-genes
colnames(weightMat)<-genes
corMat<-matrix(NA,nGenes,nGenes)
diag(corMat)<-1
rownames(corMat)<-genes
colnames(corMat)<-genes	
k<-0.95

tempTime<-system.time(for(i in 1:(nGenes-1)){	
tempCors<-mclapply((i+1):nGenes,function(j){
tempWeights<-k*sqrt(nonFailureProbs[i,]*nonFailureProbs[j,])+(1-k)
tempCorrs<-corr(t(rbind(logData[i,],logData[j,])),w=tempWeights)
return(list(corrs=tempCorrs,w=sum(tempWeights)))
},mc.cores=nCores)
tempW<-sapply(tempCors,function(temp){return(temp[["w"]])})
tempCors<-sapply(tempCors,function(temp){return(temp[["corrs"]])})
weightMat[i,(i+1):nGenes]<-tempW
weightMat[(i+1):nGenes,i]<-tempW	
corMat[i,(i+1):nGenes]<-tempCors
corMat[(i+1):nGenes,i]<-tempCors		
})[[3]]

#save correlation matrix as R object
save(file=paste(base, "_", center_gene,"_corMat.rd",sep=""),list="corMat")

#Determine correlation threshold
while (length(which(corMat[center_gene,]>cor_threshold))<min_neighbours){
cor_threshold=cor_threshold-0.01
}

print(paste("Correlation threshold determined:",cor_threshold))
friends=corMat[which(corMat[center_gene,]>cor_threshold),center_gene]

#Plot top correlations
pdf(paste(base,"_topCorrelations.pdf", sep=""))
par(mar=c(3,7,1,1))
barplot(sort(friends),hor=T,las=1,cex.names=0.8)
grid()
dev.off()

#Perform network analysis for gene of interest
cors=which(corMat[center_gene,]>cor_threshold)
geneNames=colnames(corMat[,cors])
all=geneNames
names=c()

#Find all genes that are in the network
for (j in 2:depth){
for (i in 1:length(geneNames)){
names=c(names,colnames(corMat[,which(corMat[geneNames[i],]>(cor_threshold+(0.1*j)))]))
}
names=unique(names)
geneNames=names
all=c(all,names)
names=c()
}
all=unique(all)

#Determine binary connection matrix and write it as csv
reduced=corMat[all,all]
reduced[which(reduced<cor_threshold)]=0
#reduced[which(reduced>=cor_threshold)]=1
el=getEdgeList(reduced, rm.zero=TRUE)
write.table(el, file=paste(base, "_", cor_threshold,"_",center_gene,"_cor.txt",sep=""), row.names=FALSE, col.names=TRUE, sep="\t")

