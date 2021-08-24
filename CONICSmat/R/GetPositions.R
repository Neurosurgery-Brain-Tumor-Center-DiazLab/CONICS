#' Receive genomic coordinates of a gene list
#'
#' This function allows to receive the genomic positions of a vector of genes in HUGO format.
#' @param gene_names A vector of gene names in HUGO format.
#' @param ensembl_version Version of the ENSEMBL database used to quantify gene expression data. Defaul: v87.
#' @param ignoreAlt If set to TRUE: Ignore if multiple loci are reported for a gene, pick the one from the primary assembly.
#' @keywords Chromosomal positions
#' @export
#' @examples
#' getGenePositions(gene_names=c("EGFR","PDGFRA"))

getGenePositions= function(gene_names,ensembl_version="dec2016.archive.ensembl.org",species="human", ignoreAlt = F){
  if (species=="human"){
	ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
	gene_positions <- biomaRt::getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters ='hgnc_symbol', values =gene_names, mart = ensembl)
  }
  else {
	ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host=ensembl_version)
	gene_positions <- biomaRt::getBM(attributes=c('ensembl_gene_id','mgi_symbol','chromosome_name','start_position','end_position'), filters ='mgi_symbol', values =gene_names, mart = ensembl)
  }
  if (ignoreAlt == T){
  gene_positions=gene_positions[order(gene_positions[,3],decreasing=F),]
  	}
  gene_positions=gene_positions[!duplicated(gene_positions[,2]),]
  gene_positions[which(gene_positions[,3]=="X"),3]=23
  gene_positions[which(gene_positions[,3]=="Y"),3]=24
  gene_positions[which(gene_positions[,3]=="MT"),3]=0
  gene_positions[which(nchar(gene_positions[,3])>2),3]=0
  gene_positions=gene_positions[order(as.numeric(gene_positions[,3]),decreasing=F),]
  return(gene_positions)
}
