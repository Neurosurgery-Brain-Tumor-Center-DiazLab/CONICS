#!/bin/bash

source ./CONICS.cfg

tumor_bam=${1}
snvs_normal=${2}
snvs_tumor=${3}
cnv_bed=${4}
output_base=${5}

if [ $# -ne 5 ]
then
	echo "usage: [script.sh] [directory_tumor_bam_files] [SNVs Control] [SNVs Tumor] [CNV list] [base name for output file]"
	exit 1
fi

out_dir="output_"${output_base}
if [ -d ${out_dir} ]
then
	echo "[error] output directory \"${out_dir}\" already exists"
	#exit 1
fi

mkdir ${out_dir}

echo "========= <generating BED file of heterozygous germline SNVs> ================="
#${path_to_python} ./backend/FilterSNVsFromVCF.py ${snvs_normal} ${snvs_tumor} ${cnv_bed} ./${out_dir}/${output_base}
echo "   =============================="
echo "   BED file generated"
echo "   =============================="
echo "========= <Summarizing all SNVs in each cell> ================="

#${path_to_python} ./backend/summarizeSNVs.py ${tumor_bam} ./${out_dir}/${output_base} ./${out_dir}/${output_base}_germline_snvs.bed ${genome}

echo "========= <Summary done> ================="

echo "============== <Plotting results> ====================="
${path_to_rscript} ./backend/plotVAFs.R  ${out_dir}/${output_base}_af.txt ${out_dir}/${output_base}_bf.txt ${out_dir}/${output_base}
