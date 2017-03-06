#!/bin/bash

source Empirical_validation.cfg

train_bam=${1}
test_bam=${2}
cnv_bed=${3}
output_base=${4}

if [ $# -ne 4 ]
then
	echo "usage: [script.sh] [directory for TRAINING bam files] [directory for sets of TESTING bam files] [CNV list] [base name for output file]"
	echo "*** Bam files for each testing set need to be grouped as a directory"
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
mkdir ${out_dir}/cpmMatrix
mkdir ${out_dir}/incidenceMatrix


echo ${test_bam}
echo "=================================================="
echo "calculating CPM ..."
echo "   =============================="
${path_to_python} ./backend/CalculateCoverageMatrix.py ${train_bam} ${cnv_bed} ${mappingCutoff} ${readCutoff} ${path_to_samtools}  ${path_to_bedtools} ./${out_dir}/"cpmMatrix/train"
for i in `ls ${test_bam}`
do
	echo "   =============================="
	echo "   > sample "${i}
	${path_to_python} ./backend/CalculateCoverageMatrix_empirical.py ${i} ${cnv_bed} ${mappingCutoff} ${readCutoff} ${path_to_samtools}  ${path_to_bedtools} ./${out_dir}/"cpmMatrix/"${i} ${test_bam}
	${path_to_rscript} ./backend/plot_bean_CPM_validation.R ./${out_dir}/"cpmMatrix/"${i}"_cpmMatrix.csv" ./${out_dir}/"cpmMatrix/train_cpmMatrix.csv" ./${out_dir}/"incidenceMatrix/"${i}"_result" ${fdrCutoff}
done


echo "=================================================="
echo "calculating FDR ..."

rm -rf ${out_dir}"/incidenceMatrix/"*"_result.pdf"
cd ${out_dir}"/incidenceMatrix"
touch ${output_base}"_result.txt"
ls >> ${output_base}"_result.txt"
cd ..
cd ..
	${path_to_perl} ./backend/process_empirical_Val_result.pl ${output_base} ${out_dir}"/incidenceMatrix/"${output_base}"_result.txt" ${out_dir}"/"${output_base}"_FDRs_for_Boxplot.txt"


echo "=================================================="
echo "drawing box plot ..."
${path_to_rscript} ./backend/plot_Boxplot.R ${out_dir} ${output_base} "empirical validation"
rm -rf ${out_dir}"/cpmMatrix"
rm -rf ${out_dir}"/incidenceMatrix"



