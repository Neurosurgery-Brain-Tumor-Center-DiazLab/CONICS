#!/bin/bash

source 10X_cross_validation.cfg

normal_bam=${1}
cnv_bed=${2}
output_base=${3}

if [ $# -ne 3 ]
then
	echo "usage: [script.sh] [directory_normal_bam_files] [CNV list] [base name for output file]"
	exit 1
fi

out_dir="output_"${output_base}
if [ -d ${out_dir} ]
then
	echo "[error] output directory \"${out_dir}\" already exists"
	exit 1
fi

mkdir ${out_dir}
bam_list=${out_dir}/${output_base}"_bam_list.txt"
mkdir ${out_dir}/sample_list
mkdir ${out_dir}/cpmMatrix
mkdir ${out_dir}/incidenceMatrix
touch ${bam_list}

for f in `ls  ${normal_bam}`
do
	echo "${f}" >> ${bam_list}
done

echo "=================================================="
echo "sampling ..."
	${path_to_perl} ./backend/mk_10_random_sample.pl ${bam_list} ${out_dir}

echo "=================================================="
echo "calculating CPM ..."
for i in {0..9}
do
	echo "   =============================="
	echo "   > sample "${i}
	echo "   =============================="
	echo "      train======================"
	${path_to_python} ./backend/CalculateCoverageMatrix_10X.py ${normal_bam} ${cnv_bed} ${mappingCutoff} ${readCutoff} ${path_to_samtools}  ${path_to_bedtools} ./${out_dir}/"cpmMatrix/Train_"${i} ./${out_dir}"/sample_list/Train_"${i}
	echo "      test======================="
	${path_to_python} ./backend/CalculateCoverageMatrix_10X.py ${normal_bam} ${cnv_bed} ${mappingCutoff} ${readCutoff} ${path_to_samtools}  ${path_to_bedtools} ./${out_dir}/"cpmMatrix/Test_"${i} ./${out_dir}"/sample_list/Test_"${i}

done

echo "=================================================="
echo "calculating FDR ..."
	${path_to_perl} ./backend/process_10x_crossVal_result.pl ${out_dir} ${path_to_rscript} 0.05 ${out_dir}/${output_base}"_FDRs_for_Boxplot.txt"

rm -rf ${out_dir}"/sample_list"
rm -rf ${out_dir}"/cpmMatrix"
rm -rf ${out_dir}"/incidenceMatrix/Result_"*".pdf"


echo "=================================================="
echo "drawing box plot ..."
	${path_to_rscript} ./backend/plot_Boxplot.R ${out_dir} ${output_base}



