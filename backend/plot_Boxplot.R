
args=commandArgs(TRUE)

dir <- args[1]
base <- args[2]
method <- args[3]

pdf_file <- paste(paste(dir, base, sep="/"), "_boxplot.pdf", sep="")
input_file <- paste(paste(dir, base, sep="/"), "_FDRs_for_Boxplot.txt", sep="")

pdf(pdf_file)
a <- read.table(input_file, header=FALSE)
boxplot(a,  main=paste(base, method, sep="_"), ylab="FDR")
grid()
dev.off()

