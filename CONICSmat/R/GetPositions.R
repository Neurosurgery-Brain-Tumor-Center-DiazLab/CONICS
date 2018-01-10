#' Receive genomic coordinates of a gene list
#'
#' This function allows to receive the genomic positions of a vector of genes in HUGO format.
#' @param gene_names A vector of gene names in HUGO format.
#' @param ensembl_version Version of the ENSEMBL database used to quantify gene expression data. Defaul: v87.
#' @keywords Chromosomal positions
#' @export
#' @examples
#' getGenePositions(gene_names=c("EGFR","PDGFRA"))

getGenePositions= function(gene_names,ensembl_version=87){
  ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version=ensembl_version)
  gene_positions <- biomaRt::getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters ='hgnc_symbol', values =gene_names, mart = ensembl)
  gene_positions=gene_positions[!duplicated(gene_positions[,2]),]
  length(gene_positions[,1])
  gene_positions[which(gene_positions[,3]=="X"),3]=23
  gene_positions[which(gene_positions[,3]=="Y"),3]=24
  gene_positions[which(gene_positions[,3]=="MT"),3]=0
  gene_positions[which(nchar(gene_positions[,3])>2),3]=0
  gene_positions=gene_positions[order(as.numeric(gene_positions[,3]),decreasing=F),]
  return(gene_positions)
}