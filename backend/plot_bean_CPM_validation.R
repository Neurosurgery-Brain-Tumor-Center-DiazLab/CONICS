library(beanplot)

args=commandArgs(TRUE)

tumor=read.table(args[1],header=T,row.names=1,sep="\t")
control=read.table(args[2],header=T,row.names=1,sep="\t")
outFile=args[3]
fdrTresh=as.double(args[4])

tumor=log2(tumor+1)
control=log2(control+1)
result=tumor

pdf(paste(outFile,".pdf",sep=""))

for (i in 1:ncol(tumor)){
	beanplot(tumor[,i],control[,i],col="lightblue", bw="nrd0", names=c("Tumor","Control brain"),beanlinewd=0, what=c(0,1,1,1),ylab="log2 CPM")
	names=unlist(strsplit(colnames(tumor)[i], "\\.", fixed = FALSE, perl = FALSE, useBytes = FALSE))
	print(substr(names[1], 2, nchar(names[1])))
	title_string=paste("Chr",substr(names[1], 2, nchar(names[1]))," Pos ",names[2],":",names[3]," ", names[4], sep="")
	title(title_string)
	if (substr(unlist(names)[4],1,1) =="A"){
		cdf=ecdf(control[,i])
		pval=1-cdf(tumor[,i])
		fdr=p.adjust(pval,method="BH")
		discr= ifelse(fdr<fdrTresh,1,0)
		result[,i]=discr
		
		if (sum(fdr<fdrTresh)>0){
			cutoff=min(tumor[which(fdr<fdrTresh),i])-0.1
		}
		else{
			cutoff=max(control[,i])
		}
		abline(cutoff,0,lty=2,col="blue")
	}
	else {
		cdf=ecdf(control[,i])
		pval=cdf(tumor[,i])
		fdr=p.adjust(pval,method="BH")
		discr= ifelse(fdr<fdrTresh,1,0)
		result[,i]=discr
		
		if (sum(fdr<fdrTresh)>0){
			cutoff=max(tumor[which(fdr<fdrTresh),i])+0.1
			print(substr(names[1], 2, nchar(names[1])))
			print(cutoff)
		}
		else{
			cutoff=min(control[,i])
		}
		abline(cutoff,0,lty=2,col="orange")
	}
	grid()
}
dev.off()


write.table(result,paste(outFile,".incidenceMatrix.csv",sep=""),sep="\t")
