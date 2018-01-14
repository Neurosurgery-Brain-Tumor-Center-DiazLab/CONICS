#' Receive a tab-separated expression matrix in gz format from the web.
#'
#' This function allows to receive a tab-separated expression matrix in gz format from the web, for example TCGA RNA-seq expression data from the ucsc xenabrowser.
#' @param weblink .
#' @param ensembl_version Version of the ENSEMBL database used to quantify gene expression data. Defaul: v87.
#' @keywords Chromosomal positions
#' @export
#' @examples
#' getGenePositions(gene_names=c("EGFR","PDGFRA"))

loadPublicData = function (weblink=NULL,cancer_id=NULL){
  one=c("GBM","LGG")
  two=c("https://gdc.xenahubs.net/download/TCGA-GBM/Xena_Matrices/TCGA-GBM.htseq_counts.tsv.gz","https://gdc.xenahubs.net/download/TCGA-LGG/Xena_Matrices/TCGA-LGG.htseq_counts.tsv.gz")
  res=cbind(one,two)
  if (!is.null(cancer_id)){
    weblink=res[which(one==cancer_id),2]
  }
  temp <- tempfile()
  download.file(weblink, temp)
  aap.file <- read.table(gzfile(temp), sep="\t",header=T,row.names=1)
  unlink(temp)
  return(aap.file)
}

