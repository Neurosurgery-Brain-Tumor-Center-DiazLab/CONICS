#!/bin/bash

source GoodMethod.cfg

tumor_bam=${1}
normal_bam=${2}
cnv_bed=${3}
output_base=${4}

if [ $# -ne 4 ]
then
	echo "usage: [script.sh] [directory_tumor_bam_files] [directory_normal_bam_files] [CNV list] [base name for output file]"
	exit 1
fi

out_dir="output_"${output_base}
if [ -d ${out_dir} ]
then
	echo "[error] output directory \"${out_dir}\" already exists"
	exit 1
fi

mkdir ${out_dir}


echo "========= <generating CPM matrix> ================="
echo "mapping threshold = ${mappingCutoff}"
echo "read number threshold = ${readCutoff}"
echo "=================================================="
echo "   =============================="
echo "   > Tumor"
echo "   =============================="
${path_to_python} ./backend/CalculateCoverageMatrix.py ${tumor_bam} ${cnv_bed} ${mappingCutoff} ${readCutoff} ${path_to_samtools}  ${path_to_bedtools} ./${out_dir}/${output_base}"_tumor"
echo "   =============================="
echo "   > Normal"
echo "   =============================="
${path_to_python} ./backend/CalculateCoverageMatrix.py ${normal_bam} ${cnv_bed} ${mappingCutoff} ${readCutoff} ${path_to_samtools}  ${path_to_bedtools} ./${out_dir}/${output_base}"_normal"


echo "============== <CNV calling> ====================="
${path_to_rscript} ./backend/plot_bean_CPM.R  ${out_dir}/${output_base}_tumor_cpmMatrix.csv ${out_dir}/${output_base}_normal_cpmMatrix.csv ${out_dir}/${output_base} ${fdrCutoff}


rm ${out_dir}/${output_base}_normal_rawMatrix.csv
rm ${out_dir}/${output_base}_tumor_rawMatrix.csv
rm ${out_dir}/${output_base}_normal_cpmMatrix.csv
rm ${out_dir}/${output_base}_tumor_cpmMatrix.csv


