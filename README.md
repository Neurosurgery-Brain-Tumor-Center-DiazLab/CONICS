# phyngle
*PHY*logenies and clone-specific expression from si*NGLE*-cell RNA sequencing

## Table of contents
- [Identifying CNVs from scRNA-seq](#Calling_CNV)
- [Phylogenetic tree contruction](#Constructing_Tree)
- [Intra-clone co-expression networks](#CX_Net)
- [Assessing the correlation of CNV status with single-cell expression](#Corr)
- [False discovery rate estimation: Cross validation](#10x)
- [False discovery rate estimation: Empirical testing](#Empirical)


## <a id="Calling_CNV"></a> Calling presence/absence of CNVs
Identifying large-scale CNVs in scRNA-seq. 

### Requirements
  * [Python](https://www.python.org) and [Perl](https://www.perl.org)
  * [beanplot (R package)](https://www.jstatsoft.org/article/view/v028c01)
  * [samtools](http://www.htslib.org)
  * [bedtools](http://bedtools.readthedocs.io/en/latest)  
  * Two directories, the first containing the aligned scRNA-seq data to be classified by CNV status, and a second, containing aligned scRNA-seq data to be used as a control.
  * A file contianing the genomic coordinates of the CNVs in [BED](https://genome.ucsc.edu/FAQ/FAQformat#format1) format.

### Config file
Adjust __Phyngle.cfg__ to customize the following:
  * Path to python/samtools/bedtools/Rscript
  * Thresholds for mapping-quality and read-count
  * FDR for CNV calling

### Running

  ```
  bash run_Phyngle.sh [directory for tumor] [directory for normal] [.bed file for CNV segments] [base name]
  ```
  * __[directory for tumor]__: path to directory containing aligned bam files to be tested. Example glioblastoma data, used in the manuscript, can be obtained [here](https://www.ebi.ac.uk/ega/studies/EGAS00001002185).
    
  * __[directory for normal]__: path to directory containing aligned bam files to be used as a control. Example nonmalignant brain data, used in the manuscript, can be obtained [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835) was used as an examples for the journal 
   
  * __[.bed file for CNV segments]__: tab-delimited bed file of CNV segments to be quantified.
  
      ```
      [chromosome]	[start]	[end]	[chromosome:start:end:CNV]
      ```
      
      __Note: the 4th column of the file must have the exact format shown here:
        __ (__Amp__:amplification, __Del__:deletion)
    * example (SF10281c.cnv.merged_gt1500000_20percent.bed)
  ```
  7   19533   157408385   7:19533:157408385:Amp
  9   19116859    32405639    9:19116859:32405639:Del
  ```
  * __[base name]__ : base name for output directory
  

### Output
All output files will be located in the directory __output_[base name]__.
  1. __incidenceMatrix.csv__: matrix of presence/absence for all CNVs, in individual cells
  2. __pdf__: 
    1. Read-count distribution in CNV segments. (violin plot)
    2. Hierarchical clustering of the single cells by CNV status.

![violin](images/Phyngle_violin.jpg?raw=true "violin" )
![dendrogram](images/Phyngle_dendrogram.jpg?raw=true "dendrogram" )


## <a id="Constructing_Tree"></a> Phylogenetic tree contruction
Phyngle generates the phylogenetic tree from the presence/absence profiles of CNV and visualize it. First, Phyngle clusters cells for the user-defined numbers, and then generates phylogenetic tree. Each leaf correspond the different genotypes, CNV status. This requires additional R package (Rphylip) for execution.

### Requirement
  * Rphylip (R package) with Phylip 

### Config file

Adjusting __Tree.cfg__ can change the belows.
  * Path to Rscript
  * Path to Rphylip
  
### Running
  * Before running, the path to phylip need to be added in __Tree.cfg__ file.
```
bash run_Tree.sh [CNV presence/absence matrix][number of genotypes] [base name for output file]
```
  * __[CNV presence/absence matrix]__: .incidenceMatrix.csv files. Note that exclusion of non-malignent single cells from the matrix is recommended. 
  * __[number of genotypes]__: the number of genotypes which are used for constructing trees by the aid of hierarchical clustering from the previous output files.
  * __[base name]__ : base name for output directory

### Output
__[base name]_cluster.pdf__ (phylogenetic trees) and __[base name]_cluster.txt__ would be generated in __output_[base name]__ directory. Each leaf correspond the clusters of cells which pertains different genotypes (different CNV status). Cells which belong to the different clusteres are in __[base name]_cluster.txt__ . 



![tree](images/Trees_cluster.jpg?raw=true "tree" ) 

```
cluster_1  D12,E10,F9,G3,A12,C8,C9,A3,A5,A6,C3,C2,C1,C7,H12,C4,D8,D9,A9,E4,E7,E3,F1,E1,B5,B7,E9,B3,D7,D1
cluster_2  E8,G7,G9,A7,G2,B6,E2
cluster_3  H3,A2,A4,H8,G11,F2,F3,H1,H7
cluster_4  A10,B2
cluster_5  C5
cluster_6  F8,B1
```

## <a id="CX_Net"></a> Generating co-expression network
Co-expression network could be constructed by the correlations of expressions in single cells. First, Phygle calculate [SCDE](http://hms-dbmi.github.io/scde/) -adjusted expressional correlations between genes. Then, the network which is constituted of top-correlated genes with the user-supplied target genes will be constructed.  
  

### Requirement
  * scde, PCIT, boot, parallel, raster, flashClust (R package)
  
### Config file
Adjusting __CorrelationNetwork.cfg__ can change the belows.
  * Path to Rscript
  * ncore: Number of cores for computation  (default: 12)
  * cor_threshold: Starting threshold for neighbour correlationof gene of interest (is decreased until at least x neighbours are found)(default: 0.9)
  * min_neighbours: How many direct neighbours of gene of interest should be searched for (default: 20)
  * minRawReads: How many raw reads should at least map to a gene to consider the gene (default: 100)
  * percentCellsExpressing: Percentage (0.15 =15%) of cells expressing a gene to consider the gene (default: 0.15)
  * minGenesExpr: How many genes should be expressed at least in a cell to consider the cell (default: 800)
  * depth: How deep should the gene analysis search. (2=only direct neighbor genes would be considered) (default: 2)


### Running

  ```
  bash run_CorrelationNetwork.sh [input matrix] [centered gene] [base name]
  ```
  
  * __[input matrix]__: tab-delimited file of the number of mapped reads to each gene in single cells. 
  * __[centered gene]__: a target gene of which neighbor genes are analyzed.  
  * __[base name]__ : base name for output directory

### Output
The directory __output_[base name]__ would be generated and all the output files would be located in this directory.
  1. __[base name]_[correlstion_threshold]_[gene_name].txt__ : co-expression network
  2. __[base name]_[gene_name]_corMat.rd__: Rdata of entire correlation
  3. __[base name]_topCorrelations.pdf__: bar graph of top correlations. 
  
![CXnet](images/PTEN_topCorr.jpg?raw=true "CXnet" )



## <a id="Corr"></a> Correlation of scRNA-seq and exome-Seq
The correlation between Exome-seq and single-cell RNA-seq could be assessed. Briefly, the read counts resulting from the single-cell sequencing are compared to those from the exome-seq from the same patients.  

### Requirement
  * zoo (R package)
  
### Config file
Adjusting __CompareExomeSeq_vs_ScRNAseq.cfg__ can change the belows.
  * Path to Rscript
  * window size
  
### Running

  ```
  bash run_compareExomeSeq_vs_ScRNAseq.sh [matrix for read counts] [base name for output file]
  ```
  * __[matrix for read counts]__: tab-delimited file of the number of mapped reads to each gene in Exome-seq and single-cell RNA-seq
  
      ```
      [gene] [chromosome] [start] [#read in Exome-seq(normal)] [#read in Exome-seq(tumor)] [#read in scRNA-seq(normal)] [#read in scRNA-seq(tumor)]
      ```
    * example
  ```
  
DDX11L1   1   11874   538   199   5   0
WASH7P    1   14362   4263   6541   223   45
  ```
  * __[base name]__ : base name for output directory



### Output
__Compare_[window_size].pdf__ (Box plot) would be generated in __output_[base name]__ directory. 
![compare](images/Compare_200.jpg?raw=true "compare" )


## <a id="10x"></a> False discovery rate estimation: 10x Cross validation
Phyngle could estimate false discovery rate (FDR) by 10x cross validation. In cross validation, false positive CNV calls are derived from the [normal brain controls](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835) which randomly selected 10% were used as the test sets and the remaining 90% were used for training.

### Requirement
  * beanplot (R package)
  * samtools
  * bedtools  


### Config file
Adjusting __10X_cross_validation.cfg__ can change the belows.
  * Path to python/samtools/bedtools/Rscript
  * The threshold of mapping-qualities/read-count for alignment files
  * FDR for CNV calling

### Running

  ```
  bash run_10X_cross_validation.sh [directory for test] [.bed files for CNV segments] [base name]
  ```
  * __[directory for test]__: path to directory which aligned bam files of __testing__ for FDR estimation
   
  * __[.bed file for CNV segments], [base name]__ : same as above.
;
### Output
Box plot of 10 FDRs resulting from each pooled sample would be generated (__[base name]_boxplot.pdf__) in  __output_[base name]__ directories.
![10X](images/10X_boxplot.jpg?raw=true "10Xval_Test" )

## <a id="Empirical"></a> False discovery rate estimation: Empirical testing
The FDR could be estimated by empirical testing. The number of false positive CNV calls was calculated from non-malignant [fetal brain dataset](http://dx.doi.org/10.1016/j.cell.2015.09.004) which are independent on the [training set](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835)

### Requirement
  * beanplot (R package)
  * samtools
  * bedtools  
  
### Config file
Adjusting __Empirical_validation.cfg__ can change the belows.
  * Path to python/samtools/bedtools/Rscript
  * The threshold of mapping-qualities/read-count for alignment files
  * FDR for CNV calling

### Running

  ```
  bash run_empirical_validation.sh [directory for train] [directory for test] [.bed files for CNV segments] [base name]
  ```

  * __[directory for train]__: path to directory which aligned bam files for __training__ for FDR estimation
    
  * __[directory for test]__: path to directory which aligned bam files of __testing__ for FDR estimation
   
  * __[.bed file for CNV segments], [base name]__ : same as above.

### Output
Box plot of FDRs would be generated (__[base name]_boxplot.pdf__) in  __output_[base name]__ directories.
![empirical](images/Empirical_boxplot.jpg?raw=true "empirical_Test" )
