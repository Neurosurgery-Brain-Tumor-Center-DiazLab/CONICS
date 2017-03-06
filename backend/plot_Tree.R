library(Rphylip)
require(methods)

args=commandArgs(TRUE)
input <- args[1]
k <- args[2]
phylip_dir <- args[3]
tag <- args[4]

output_cell <- paste(tag, "_cluster.pdf", sep="")
output_txt <- paste(tag, "_clusters.txt", sep="")



x <- read.table(input, sep="\t", row.names=1, header=T)
clusters <- as.integer(k)
dis=dist(x, 'binary')
hc=hclust(dis)
hc_label=cutree(hc, k=clusters)
clone=sapply(unique(hc_label), function(i) which(hc_label==i) )
n = length(clone)
hc_genotype = matrix(0, n, dim(x)[2])
cluster_name = rep("a", n)
cluster_member = rep("a", n)
colnames(hc_genotype)=colnames(x)
for (i in 1:n){
	cluster_name[i] = paste("cluster_", i, sep="")
	cluster_member[i] = paste(cluster_name[i], capture.output(cat(attributes(clone[[i]])$names, sep=",")), sep="  ")
	idx= clone[[i]]
	if ( length(idx)==1 ){
		hc_genotype[i,]= as.matrix(x[idx,])
	}
	else{
		hc_genotype[i,]= as.numeric( colMeans(as.matrix(x[idx,])) > 0.5 )
	}
}
rownames(hc_genotype) <- cluster_name
root=rep(0,ncol(hc_genotype))
x=rbind(hc_genotype,root)
row.names(x)[nrow(x)]="Normal cell"
dis=dist(x, "manhattan")
tree2=Rfitch(dis, path=phylip_dir, outgroup="Normal cell")
pdf(output_cell)
plot(unroot(tree2))

dev.off()
write(cluster_member, file=output_txt)
