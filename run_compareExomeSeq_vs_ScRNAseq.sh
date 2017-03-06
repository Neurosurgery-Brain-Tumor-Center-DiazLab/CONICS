#!/bin/bash

source CompareExomeSeq_vs_ScRNAseq.cfg

input=${1}
output_base=${2}

if [ $# -ne 2 ]
then
	echo "usage: [script.sh] [input_matrix] [base name for output file]"
	exit 1
fi

out_dir="output_"${output_base}
if [ -d ${out_dir} ]
then
	echo "[error] output directory \"${out_dir}\" already exists"
	exit 1
fi

mkdir ${out_dir}


echo "=================================================="
echo "window size = ${window_size}"
echo "calculating..."
${path_to_rscript} ./backend/CompareExomeSeq_vs_ScRNAseq.R ${input} ${window_size} ${output_base}
echo "done";
echo "=================================================="





