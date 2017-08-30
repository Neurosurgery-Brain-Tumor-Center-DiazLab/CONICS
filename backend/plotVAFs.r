library(ggplot2)
library(beanplot)

#setwd("/Users/soere/Documents/UCSF/Analysis/Phyngle/SF10345/")

#a=read.table("res_sf10345_af.txt",header=T,row.names=1,sep="\t",check.names=F)
#b=read.table("res_sf10345_bf.txt",header=T,row.names=1,sep="\t",check.names=F)

args=commandArgs(TRUE)
a=read.table(args[1],header=T,row.names=1,sep="\t",check.names=F)
b=read.table(args[2],header=T,row.names=1,sep="\t",check.names=F)

a =a[,colSums(is.na(a))<nrow(a)]
b=b[,colnames(a)]

regions=unlist(strsplit(colnames(a),"_",fixed=TRUE))[seq(2,(6*ncol(a)),6)]
normal_baf=as.numeric(unlist(strsplit(colnames(a),"_",fixed=TRUE))[seq(5,(6*ncol(a)),6)])
tumor_baf=as.numeric(unlist(strsplit(colnames(a),"_",fixed=TRUE))[seq(6,(6*ncol(a)),6)])
regions_idf=names(table(regions))[which(table(regions)>25)]

av_mat=matrix(ncol=length(regions_idf),nrow=nrow(a)+2)
colnames(av_mat)=regions_idf
rownames(av_mat)=c(rownames(a),"Exome_Norm","Exome_Tum")

for (i in 1:length(regions_idf)){
  af=rowSums(a[,which(regions==regions_idf[i])],na.rm = T)
  bf=rowSums(b[,which(regions==regions_idf[i])],na.rm = T)
  both=af+bf
  perc=length(which(both>20))/length(both);perc
  if (perc<0.4){
    av_mat[,i]=NA
  }
  else {
    r=bf/(af+bf)
    r[which(both<8)]=NA
    tum=mean(tumor_baf[which(regions==regions_idf[i])])
    norm=mean(normal_baf[which(regions==regions_idf[i])])
    av_mat[,i]=c(r,norm=norm,tum=tum)
  }
}
av_mat=av_mat[,colSums(is.na(av_mat))<nrow(av_mat)]
plot(av_mat[,1])

#Problem: Most TAMs do not have coverage at any SNVs, eg C7
# Why do you have to sum: Because it's random which of the two alleles will be measured in a single cell if locus heterozygous
colnames(a)[which(regions==regions_idf[i])][which(a["D6",which(regions==regions_idf[i])]>0)]
colnames(b)[which(b["C7",]>0)]
colnames(a)[which(a["B6",]>4)]


#Plot BAF per region
sel=!is.na(match(regions,regions_idf[c(1,3,4,5)]))
e=as.data.frame(tumor_baf[sel])
e$expr=tumor_baf[sel]
e$st=regions[sel]
p <- ggplot(e, aes(st, expr))
p + geom_boxplot() + geom_jitter(width = 0.2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ylab("BAF Exome-seq Tumor") +xlab("") +ylim(0,0.5)

#Plot heatmap of z-scores
par(mar=c(10,5,5,5))
heatmap.2(av_mat[which(rowSums(av_mat[,c(1,3,4,5)])>0),c(1,3,4,5)],na.color="blue",margins=c(13,8),srtCol=45)


