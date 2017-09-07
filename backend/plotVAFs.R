library(gplots)
library(ggplot2)

#setwd("/Users/soere/Documents/UCSF/Analysis/Phyngle/SF10281/")

#a=read.table("sf10281baf_af.txt",header=T,row.names=1,sep="\t",check.names=F)
#b=read.table("sf10281baf_bf.txt",header=T,row.names=1,sep="\t",check.names=F)

#Reads file with A allele frequency anf file with B allele frequency for every SNV in every region
args=commandArgs(TRUE)
a=read.table(args[1],header=T,row.names=1,sep="\t",check.names=F)
b=read.table(args[2],header=T,row.names=1,sep="\t",check.names=F)
ofi=args[3]

#Filter SNVs without coverage
a =a[,colSums(is.na(a))<nrow(a)]
b=b[,colnames(a)]

#Create vectors of 1: Regions of copy number alterations; 2 BAF from control exome-seq; 3 BAF from tumor exome-seq;
regions=unlist(strsplit(colnames(a),"_",fixed=TRUE))[seq(2,(6*ncol(a)),6)]
normal_baf=as.numeric(unlist(strsplit(colnames(a),"_",fixed=TRUE))[seq(5,(6*ncol(a)),6)])
tumor_baf=as.numeric(unlist(strsplit(colnames(a),"_",fixed=TRUE))[seq(6,(6*ncol(a)),6)])

# Only use regions with at least 25 SNVs 
regions_idf=names(table(regions))[which(table(regions)>25)]
av_mat=matrix(ncol=length(regions_idf),nrow=nrow(a)+2)
colnames(av_mat)=regions_idf
rownames(av_mat)=c(rownames(a),"Exome_Norm","Exome_Tum")

#Loop for every region
for (i in 1:length(regions_idf)){
	#Sum over all a alleles in each cell
  af=rowSums(a[,which(regions==regions_idf[i])],na.rm = T)
	#Sum over all B alleles in each cell
  bf=rowSums(b[,which(regions==regions_idf[i])],na.rm = T)
  both=af+bf
	#Do at least 40 percent of cells have more than 20 reads covering SNVs in the region?
  perc=length(which(both>20))/length(both);perc
  if (perc<0.4){
    av_mat[,i]=NA
  }
  else {
	#Calculate BAF for each cell
    r=bf/(af+bf)
	#NA if less than 8 reads, since BAF will be unreliable
    r[which(both<8)]=NA
	#Average tumor BAF in the region from exome-seq
    tum=mean(tumor_baf[which(regions==regions_idf[i])])
	#Average control BAF in the region from exome-seq
    norm=mean(normal_baf[which(regions==regions_idf[i])])
    av_mat[,i]=c(r,norm=norm,tum=tum)
  }
}
#Remove columns with only NA
av_mat=av_mat[,colSums(is.na(av_mat))<nrow(av_mat)]
#plot(av_mat[,1])

#Plot BAF per region
pdf(paste(ofi,"BAF_Tumor.pdf",sep=""))
sel=!is.na(match(regions,regions_idf[c(1,3,4,5)]))
e=as.data.frame(tumor_baf[sel])
e$expr=tumor_baf[sel]
e$st=regions[sel]
p <- ggplot(e, aes(st, expr))
p + geom_boxplot() + geom_jitter(width = 0.2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ylab("BAF Exome-seq Tumor") +xlab("") +ylim(0,0.5)
dev.off()

#Plot heatmap of z-scores
pdf(paste(ofi,"Heatmap_BAF.pdf",sep=""))
par(mar=c(10,5,5,5))
heatmap.2(av_mat[which(rowSums(av_mat[,])>0),],na.color="blue",margins=c(13,8),srtCol=45)
dev.off()

