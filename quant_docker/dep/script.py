# docker build -t kallisto_sc_quant . ; docker run kallisto_sc_quant

import StringIO
import os
import sys
import subprocess

def quant(index, output, threads, f1, f2):
    call=["/dep/kallisto", "quant", "-i", index, "-o", output, "--threads", str(threads), f1, f2]
    string = open("/dep/_tmp", "w+")
    try:
        subprocess.call(call, stdout = string, stderr = subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print "CalledProcessError occured"
    string.seek(0, 0)
    print string.read()
    string.close()

quant("/quake/kallisto/index.idx", "/quake/kallisto/out/", 4, "/quake/kallisto/srr_1.fastq", "/quake/kallisto/srr_2.fastq")
