# __CONICS__
*CONICS*: *CO*py-*N*umber analysis *I*n single-*C*ell RNA-*S*equencing

The CONICS paper has been accepted for publication in Bioinformatics. Check it out [here](https://doi.org/10.1093/bioinformatics/bty316) !

## Table of contents
- [CONICSmat - Identifying CNVs from scRNA-seq with only a count table (R tutorial)](#Calling_CNV_R)
- [Identifying CNVs from scRNA-seq using aligned reads and a control dataset (advanced users)](#Calling_CNV)
- [Integrating the minor-allele frequencies of point mutations (advanced users)](#B_allele_freqs)
- [Phylogenetic tree contruction](#Constructing_Tree)
- [Intra-clone co-expression networks](#CX_Net)
- [Assessing the correlation of CNV status with single-cell expression](#Corr)
- [False discovery rate estimation: Cross validation](#10x)
- [False discovery rate estimation: Empirical testing](#Empirical)


## <a id="Calling_CNV_R"></a> CONICSmat - Identifying CNVs from scRNA-seq using a count table 

CONICSmat is an R package that can be used to identify CNVs in single cell RNA-seq data from a gene expression table, without the need of an explicit normal control dataset. A tutorial on how to use CONICSmat, and a Smart-Seq2 dataset, can be found [on the CONICSmat Wiki page [CLICK here]](https://github.com/diazlab/CONICS/wiki/Tutorial---CONICSmat;---Dataset:-SmartSeq2-scRNA-seq-of-Oligodendroglioma).

![overview](images/conics_overview.jpg "overview" )
<i> Visualizations of scRNA-seq data from [Oligodendroglioma](https://www.nature.com/articles/nature20123/) (Tirosh et al., 2016) generated with CONICSmat. </i>


<br>
<br>

## <a id="Calling_CNV"></a> CONICS - Identifying CNVs from scRNA-seq with alignment files 

### Requirements
  * [Python](https://www.python.org) and [Perl](https://www.perl.org)
  * [beanplot](https://www.jstatsoft.org/article/view/v028c01)
  * [samtools](http://www.htslib.org)
  * [bedtools](http://bedtools.readthedocs.io/en/latest)  IMPORTANT: Bedtools >2.2.5 is needed in order to correctly calculate the coverage using CONICS.
  * Two directories, the first containing the aligned scRNA-seq data to be classified by CNV status, and a second, containing aligned scRNA-seq data to be used as a control.
  * A file contianing the genomic coordinates of the CNVs in [BED](https://genome.ucsc.edu/FAQ/FAQformat#format1) format.

### Config file
Adjust __CONICS.cfg__ to customize the following:
  * Path to python/samtools/bedtools/Rscript
  * Thresholds for mapping-quality and read-count
  * FDR for CNV calling
  
### Test data 

A set of tumor cells from three glioblastoma patients and a normal brain control, as well as a file with genomic coordinate of large-scale CNVs are available [here](https://goo.gl/1dyafr). 

### Running

  ```
  bash run_CONICS.sh [directory for tumor] [directory for normal] [.bed file for CNV segments] [base name]
  ```
  * __[directory for tumor]__: path to directory containing aligned bam files to be tested. Example glioblastoma data, used in the manuscript, can be obtained [here](https://www.ebi.ac.uk/ega/studies/EGAS00001002185).
    
  * __[directory for normal]__: path to directory containing aligned bam files to be used as a control. Example nonmalignant brain data, used in the manuscript, can be obtained [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835) was used as an examples for the journal 
   
  * __[BED file for CNV segments]__: tab-delimited bed file of CNV segments to be quantified.
  
      ```
      [chromosome]	[start]	[end]	[chromosome:start:end:CNV]
      ```
      
      Note: the 4th column of the file must have the exact format shown here:(__Amp__: amplification, __Del__: deletion)
  ```
  7   19533   157408385   7:19533:157408385:Amp
  9   19116859    32405639    9:19116859:32405639:Del
  ```
  * __[base name]__ : base name for output directory
  

### Output
All output files will be located in the directory __output_[base name]__.
  1. __incidenceMatrix.csv__: matrix of presence/absence for all CNVs, in individual cells
  2. Read-count distribution in CNV segments. (violin plot)
  3. Hierarchical clustering of the single cells by CNV status.

![violin](images/Phyngle_violin.jpg?raw=true "violin" )
![dendrogram](images/Phyngle_dendrogram.jpg?raw=true "dendrogram" )

## <a id="B_allele_freqs"></a> Integrating estimates of point-mutation minor-allele frequencies
Regions of copy-number alteration will show a drop in the frequency of reads quantifying the minor allele. Averaged over large regions of copy-number alteration, this provides an additional metric to increase confidence in single-cell CNV-calls.

### Requirements
  * [Python](https://www.python.org) and [R](https://cran.r-project.org/mirrors.html)
  * [bam-readcount](https://github.com/genome/bam-readcount) 
  * [gplots](https://cran.r-project.org/web/packages/gplots/index.html) and [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
  * One directory containing the aligned tumor scRNA-seq data to be classified
  * Two variant [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) files from exome-seq of (blood) control and tumor tissue, eg generated with the [GATK](https://software.broadinstitute.org/gatk/) toolkit.
  * A file contianing the genomic coordinates of the CNVs in [BED](https://genome.ucsc.edu/FAQ/FAQformat#format1) format.
  
### Config file
Adjust __CONICS.cfg__ to customize the following:
  * Path to python/bam-readcount/Rscript
  * Path to genome which reads were aligned to (FASTA format)

### Running

 ```
  bash run_BAf_analysis.sh [directory for tumor] [VCF file for normal exome-seq] [VCF file for tumor exome-seq] [BED file for CNV segments] [base name]
  ```
  * __[directory for tumor]__: path to directory containing aligned bam files to be tested. Example glioblastoma data, used in the manuscript, can be obtained [here](https://www.ebi.ac.uk/ega/studies/EGAS00001002185).
    
  * __[VCF file for normal exome-seq]__: Vcf file containing mutations for a control exome-seq, e.g. from blood of the patient. This file can be generated with tools like GATK toolkit.
   
  * __[VCF file for tumor exome-seq]__: Vcf file containing mutations detected in exome-seq of the tumor. This file can be generated with tools like GATK toolkit.
  
 * __[BED file for CNV segments]__: tab-delimited bed file of CNV segments to be quantified.
  
  * __[base name]__ : base name for output directory
  

### Output
All output files will be located in the directory __output_[base name]__
  1. ___germline-snvs.bed__: BED file containing position and BAFs from Exome-seq, generated in step 1.
  2. ___af.txt__: TAB separated table containing the counts for the A allele at each locus in each cell, generated in step 2
  3. ___bf.txt__: TAB separated table containing the counts for the B allele at each locus in each cell, generated in step 2
  4. baf_hist.pdf Hierarchical clustering of the average B allele frequency in each of the loci altered by copy number for each cell, generated in step 3
  
![heatmap](images/heatmap_sf10345.jpg?raw=true "B-allele heatmap" )

## <a id="Constructing_Tree"></a> Phylogenetic tree contruction
CONICS can generate a phylogenetic tree from the CNV incidence matrix, using the Fitch-Margoliash algorithm. Other phylogenetic reconstruction algorithms can be applied, using the incidence matrix as a starting point.

### Requirements
  * [Rscript](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html)
  * [Rphylip](https://cran.r-project.org/web/packages/Rphylip/index.html)
  * [Phylip](http://evolution.genetics.washington.edu/phylip.html) 

### Config file

Adjust __Tree.cfg__ to change the following.
  * Path to Rscript
  * Path to Rphylip
  
### Running
  * Before running, set the path to Phylip in __Tree.cfg__ file.
```
bash run_Tree.sh [CNV presence/absence matrix][number of genotypes] [base name for output file]
```
  * __[CNV presence/absence matrix]__: .incidenceMatrix.csv files. 
  * __[number of genotypes]__: the number of genotypes to model
  * __[base name]__ : base name for output directory

### Output
__cluster.pdf__ (phylogenetic trees) and __cluster.txt__ will be generated in the output directory. Each leaf corresponds to a clusters of cells with a common genotype. Cluster assignments for each cell will be in __cluster.txt__. 



![tree](images/Trees_cluster.jpg?raw=true "tree" ) 

```
cluster_1  D12,E10,F9,G3,A12,C8,C9,A3,A5,A6,C3,C2,C1,C7,H12,C4,D8,D9,A9,E4,E7,E3,F1,E1,B5,B7,E9,B3,D7,D1
cluster_2  E8,G7,G9,A7,G2,B6,E2
cluster_3  H3,A2,A4,H8,G11,F2,F3,H1,H7
cluster_4  A10,B2
cluster_5  C5
cluster_6  F8,B1
```

## <a id="CX_Net"></a> Intra-clone co-expression networks
CONICS can construct the local co-expression network of a given gene, based on correlations across single cells.  
  

### Requirements
  * [scde](http://hms-dbmi.github.io/scde)
  * [PCIT](https://cran.r-project.org/web/packages/PCIT/index.html)
  * [boot](https://cran.r-project.org/web/packages/boot/)
  * [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)
  * [raster](https://cran.r-project.org/web/packages/raster/)
  * [flashClust](https://cran.r-project.org/web/packages/flashClust/index.html)
  
### Config file
Adjust __CorrelationNetwork.cfg__ to configure the following:
  * Path to Rscript
  * ncore: Number of cores (default: 12)
  * cor_threshold: Starting threshold to construct the co-expression network (default: 0.9)
  * min_neighbours: How many direct neighbours of gene of interest should be analyzed (default: 20)
  * minRawReads: How many raw reads should map to a gene for it to be included (default: 100)
  * percentCellsExpressing: Percentage (0.15 =15%) of cells expressing a gene for it to be included (default: 0.15)
  * minGenesExpr: How many genes should be expressed in a cell for it to be included (default: 800)
  * depth: How deep should the gene analysis search. (2=only direct neighbor genes would be considered) (default: 2)


### Running

  ```
  bash run_CorrelationNetwork.sh [input matrix] [centered gene] [base name]
  ```
  
  * __[input matrix]__: tab-delimited file of read counts for each gene (rows), for each cell (columns). 
  * __[centered gene]__: a target gene of which neighbor genes are analyzed.  
  * __[base name]__ : base name for output directory

### Output
All the output files will be located in __output__.
  1. __[correlstion_threshold]_[gene_name].txt__ : co-expression network
  2. __[gene_name]corMat.rd__: Rdata containing the adjusted correlation matrix
  3. __topCorrelations.pdf__: bar graph of top correlations. 
  
![CXnet](images/PTEN_topCorr.jpg?raw=true "CXnet" )



## <a id="Corr"></a> Assessing the correlation of CNV status with single-cell gene-expression 

### Requirements
  * [zoo](https://cran.r-project.org/web/packages/zoo/index.html)
  
### Config file
Adjust __CompareExomeSeq_vs_ScRNAseq.cfg__ to set the following:
  * Path to Rscript
  * window size for assessing CNV status
  
### Running

  ```
  bash run_compareExomeSeq_vs_ScRNAseq.sh [matrix for read counts] [base name for output file]
  ```
  * __[matrix for read counts]__: tab-delimited file of the number of mapped reads to each gene in the DNA sequencing and in scRNA-seq. Genes on each chromosome should be ordered by their chromosomal position.
  
      ```
      [gene] [chromosome] [start] [#read in DNA-seq(normal)] [#read in DNA-seq(tumor)] [#read in scRNA-seq(normal)] [#read in scRNA-seq(tumor)]
      ```
    * example
  ```
  
DDX11L1   1   11874   538   199   5   0
WASH7P    1   14362   4263   6541   223   45
  ```
  * __[base name]__ : base name for output directory



### Output
__Compare_[window_size].pdf__ (Box plot) will be generated in the output directory. 
![compare](images/Compare_200.jpg?raw=true "compare" )


## <a id="10x"></a> False discovery rate estimation: Cross validation
CONICS can estimate false discovery rate via 10-fold cross-validation, using the user-supplied control scRNA-seq dataset. For example, in the manuscript cross validation was performed using [normal brain controls](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835).

### Requirements
  * [beanplot](https://www.jstatsoft.org/article/view/v028c01)(R package)
  * [samtools](http://www.htslib.org)
  * [bedtools](http://bedtools.readthedocs.io/en/latest)  


### Config file
Adjust __10X_cross_validation.cfg__ to set the following:
  * Paths to python/samtools/bedtools/Rscript
  * Thresholds for mapping-quality and read-count.
  * FDR for CNV calling

### Running

  ```
  bash run_10X_cross_validation.sh [directory for control scRNA-seq] [.bed file containing CNV segments] [base name]
  ```
  * __[directory for test]__: path to directory containing the aligned BAM files of the scRNA-seq control data.
   
  * __[.bed file for CNV segments], [base name]__ : same as described in run_CONICS.sh
;
### Output
Box plot of 10 FDRs resulting from each pooled sample would be generated (__boxplot.pdf__) in the output directory.
![10X](images/10X_boxplot.jpg?raw=true "10Xval_Test" )

## <a id="Empirical"></a> False discovery rate estimation: Empirical testing
FDRs can also be estimated by empirical testing. In the manuscript, the number of false positive CNV calls was calculated using a non-malignant [fetal brain dataset](http://dx.doi.org/10.1016/j.cell.2015.09.004). These data are independent from the [training set](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835)

### Requirements
  * [beanplot](https://www.jstatsoft.org/article/view/v028c01)
  * [samtools](http://www.htslib.org)
  * [bedtools](http://bedtools.readthedocs.io/en/latest)
  
### Config file
Adjust __Empirical_validation.cfg__ to change the following:
  * Paths to python/samtools/bedtools/Rscript
  * Thresholds for mapping-quality and read count
  * FDR for CNV calling

### Running

  ```
  bash run_empirical_validation.sh [directory for train] [directory for test] [.bed files for CNV segments] [base name]
  ```

  * __[directory for train]__: path to directory containing aligned bam files of scRNA-seq data used as a control to call CNVs
    
  * __[directory for test]__: path to directory containing aligned bam files of scRNA-seq data known not to have CNVs, used as a gold standard.
   
  * __[BED file for CNV segments]__ : same as described in run_CONICS.sh

### Output
Box plot of FDRs will be generated (__boxplot.pdf__) in the output directory.
![empirical](images/Empirical_boxplot.jpg?raw=true "empirical_Test" )
