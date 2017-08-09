#!/bin/bash

if [ "$#" != 5 ]; then
    echo "Usage: after docker command, --args kallisto_index output_folder num_threads fastq_input_folder"
    echo "fastq folder must already exist"
fi

/dep/kallisto quant -i $1 -o $2 -t $3 ${4}/1.fastq ${5}/2.fastq