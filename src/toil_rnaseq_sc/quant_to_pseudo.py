import os
import sys
import argparse
from toil.job import Job

def quant_to_pseudo(job, input_dir, output_dir):
    """
    Convert output of kallisto_sc_quant docker image (individual abundance.tsv's for each cell)
        to the same output format as kallisto pseudo (ec, cells, tsv) as described here:
            https://groups.google.com/forum/#!topic/kallisto-sleuth-users/PRTG9yPVBzg
    In the converted output, the ec file is one-to-one since each 'equivalence class' is in fact a transcript.
    The transcript numbers are based on the order of transcripts in the abundance.tsv's and are one-indexed.
    The one-indexing means that it is possible to access transcripts directly without removing the header of an abundance tsv.
    It is assumed that all the abundance.tsv's used the same index and have the same transcripts in the same order.
    The actual transcript names are not exported; use any of the original abundance.tsv's to obtain the transcript names.
    """
    ec_path = os.path.join(output_dir, "matrix.ec")
    cells_path = os.path.join(output_dir, "matrix.cells")
    tsv_path = os.path.join(output_dir, "matrix.tsv")
    input_files = sorted([os.path.join(input_dir, file) for file in os.listdir(input_dir)])
    
    if len(input_files) == 0:
        raise IOError("Empty input directory")
    
    # Write the equivalence class mapping. Every "equivalence class" is actually a transcript...
    # ...so it maps to the index of the transcript in the first TSV.
    # This assumes that all the TSVs have the same transcripts in the same order,
    with open(ec_path, "w") as ec_map, open(input_files[0], "r") as input_file:
        for index, line in enumerate(input_file):
            ec_map.write(str(index)+"\t"+str(index)+"\n")

    # Write the cell names. It is assumed that each TSV is associated with a single cell,
    # and the name of the TSV is the ID of that cell.
    with open(cells_path, "w") as cell_list:
        for input_file in input_files:
            file_id = os.path.splitext(os.path.basename(input_file))[0]
            cell_list.write(file_id+"\n")

    # Write the sparse matrix. It is assumed that the count data is in column 4 of the TSV.
    with open(tsv_path, "w") as sparse_tsv:
        for file_num, input_file in enumerate(input_files):
            with open(input_file, "r") as input:
                for col_num, line in enumerate(input):
                    if col_num == 0: continue # first line has labels only
                    counts = line.split("\t")[3]
                    if counts == "0": continue # sparse matrix -- record only nonzero values
                    icounts = int(counts)
                    if (job is not None): job.fileStore.logToMaster(str(counts) + " vs " + str(icounts))
                    sparse_tsv.write(str(col_num)+"\t")
                    sparse_tsv.write(str(file_num)+"\t")
                    sparse_tsv.write(str(icounts)+"\n")
