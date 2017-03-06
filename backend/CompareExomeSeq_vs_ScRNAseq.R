require(zoo)

args=commandArgs(TRUE)

a=read.table(args[1],header=T,row.names=1,sep="\t")
windowsize=as.numeric(args[2])
patient=args[3]


chrom=a[,1:2]
chrom[,1]=sub("\\;.*","\\",chrom[,1])
chrom[,2]=sub("\\;.*","\\",chrom[,2])
a=a[,c(3:ncol(a))]

#normalization
sdepth=apply(a,2,function(x){sum(x)})
norm=t( t(a) / sdepth*1000000)

#filter genes
norm=norm[which(rowSums(norm[,c(1,2)])>5),]
norm=norm[which(rowSums(norm[,c(3,4)])>5),]

#log transformation
norm=log2(norm+1)
chrom=chrom[row.names(norm),]

#generate matrix with one row to append results of each chromosome
windowNorm=matrix(ncol=ncol(norm))

#iterate over genes from first 20 chromosomes, not 22 or 21, too few genes
for (i in 1:19){
genes_on_chr=which(chrom[,1]==i)
print(length(genes_on_chr))
chromNorm=rollapply(norm[genes_on_chr,], width = windowsize, by = 1, FUN = mean,by.column=T)
windowNorm=rbind(windowNorm,chromNorm)
}

#remove initial NA colum
windowNorm=windowNorm[c(2:nrow(windowNorm)),]

#calculate fold change
r1=windowNorm[,2]-windowNorm[,1]
r2=windowNorm[,4]-windowNorm[,3]

#median center
r1=r1-median(r1)
r2=r2-median(r2)

pdf(paste("./output_", patient, "/", patient,"_",windowsize,".pdf",sep=""))

boxplot(r2[which(r1<(-.5))],r2[which(r1<.3&r1>(-.5))],r2[which(r1>(.3))],outline=F,names=c("Loss", "Neutral","Gain"),ylab="Log2 expression ratio tumor/control",boxwex=.3,main=paste(patient," PCC:",round(cor(r1,r2),3)," Windowsize: ",windowsize),xlab="Copy number status of locus")
grid()

plot(r1,r2,xlab="Log2 ratio Exome-seq",ylab="Log2 ratio RNA-seq")

dev.off()


