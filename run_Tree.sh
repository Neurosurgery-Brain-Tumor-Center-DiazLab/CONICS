#!/bin/bash

source Tree.cfg

input_mat=${1}
cluster=${2}
output_base=${3}

if [ $# -ne 3 ]
then
	echo "usage: [script.sh] [CNV presence/absence matrix][number of genotypes] [base name for output file]"
	exit 1
fi

out_dir="output_"${output_base}
if [ ! -d ${out_dir} ]
then
	mkdir ${out_dir}
fi


${path_to_rscript} ./backend/plot_Tree.R  ${input_mat} ${cluster} ${path_to_phylip} ${out_dir}/${output_base}






