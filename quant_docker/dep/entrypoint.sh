#!/bin/bash

if [ "$#" != 4 ]; then
    echo "Usage: after docker command, pass kallisto_index output_folder num_threads fastq_input_folder"
    echo "\tfastq folder must already exist"
else
    /dep/kallisto quant -i $1 -o $2 -t $3 ${4}/1.fastq ${4}/2.fastq
fi
