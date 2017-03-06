#!/bin/bash

source CorrelationNetwork.cfg

expression=${1}
gene=${2}
output_base=${3}

if [ $# -ne 3 ]
then
	echo "usage: [script.sh] [input matrix] [gene] [base name for output file]"
	exit 1
fi

out_dir="output_"${output_base}
if [ -d ${out_dir} ]
then
	echo "[error] output directory \"${out_dir}\" already exists"
	exit 1
fi

mkdir ${out_dir}


echo "========= <calculating co-expression> ================="

${path_to_rscript} ./backend/CorrelationNetwork.R ${gene} ${expression} ${ncore} ${cor_threshold} ${min_neighbours} ${minRawReads} ${percentCellsExpressing} ${minGenesExpr} ${depth} ${out_dir} ${output_base}
