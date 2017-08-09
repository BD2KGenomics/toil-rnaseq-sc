#!/bin/bash

if [ "$#" != 4 ]; then
    echo "Usage: after docker command, pass kallisto_index output_folder num_threads fastq_input_folder"
    echo "\tfastq folder must already exist"
else
    for id in $(ls -1 $4 | sed s/_.*// | uniq); do
        mkdir ${2}/${id}
        /dep/kallisto quant -i $1 -o ${2}/${id} -t $3 ${4}/${id}_1.fastq ${4}/${id}_1.fastq
    done
fi
