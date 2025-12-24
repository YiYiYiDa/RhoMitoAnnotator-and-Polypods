# Date: 20251128
# Author: Linchun Shi, Yanda Zhu
# Copyright: (c) 2025 Linchun Shi, Yanda Zhu. All rights reserved.

#################################################
import sys
import os
import subprocess
import fnmatch
from Bio import SeqIO

###############################################
print ("Usage: program dataDir(input) dataMerge(output)")
print ("Step02:\tmerge same sample to subdir\n##########")

'''
Purpose:
------------
        * Merge (concatenate) all R1 and R2 FASTQ files that belong to the same sample into a single R1 file and a single R2 file, respectively.
        ** Sample identity is taken from the sub-directory name that was created in the previous organisation step (step01).

Usage:
------------
        python step02_merge_fastq.py /path/to/DataDir(from step01) /path/to/DataMergeDir

Workflow:
------------
        1. Validate input and output directories.
        2. For every sample sub-folder found in the input directory:
                (1) create the same sub-folder in the output directory,
                (2) cat all *R1*/*R2* files  -> <out>/<sample>/<sample>.R1/R2.fastq.gz,
        3. Compare cumulative disk usage of input vs. merged files to confirm

Notes:
------------
        * The script skips re-merging if the target files already exist and are non-empty, allowing safe re-runs.
        * Shell command `cat` is used for concatenation, ensure the format of files are "xxx.R1/R2.xxx".
        * Disk-usage check relies on `du -h -c`; numbers are human-readable.
'''

########################################

# def make_dir(outdir):
#         if not os.path.exists(outdir):
#                 os.mkdir(outdir)


def merge_same_sample_fastq(workdir, outdir):
        for subdir in [it for it in os.listdir(workdir) if os.path.isdir(workdir+"/"+it)]:
                if not os.path.exists(outdir+"/"+subdir):
                        os.mkdir(outdir+"/"+subdir)

                outR1File = outdir+"/"+subdir+"/"+subdir+".R1.fastq.gz"
                outR2File = outdir+"/"+subdir+"/"+subdir+".R2.fastq.gz"
                if not os.path.exists(outR1File) or os.path.getsize(outR1File) == 0:
                        cmdR1 = f"cat {workdir}/{subdir}/*R1* > {outR1File}"
                        cmdR2 = f"cat {workdir}/{subdir}/*R2* > {outR2File}"
                        P1 = subprocess.Popen(cmdR1, shell=True)
                        P1.wait()
                        P2 = subprocess.Popen(cmdR2, shell=True)
                        P2.wait()
                        print (f"{subdir} merge success")
                else:
                        print (f"{subdir} already merged")


def count_data_size(workdir, outdir):
        for subdir in [it for it in os.listdir(workdir) if os.path.isdir(workdir + "/" + it)]:
                outR1File = outdir + "/" + subdir + "/" + subdir + ".R1.fastq.gz"
                outR2File = outdir + "/" + subdir + "/" + subdir + ".R2.fastq.gz"

                inputR1FileSizeCmd = f"du -h -c {workdir}/{subdir}/*R1* |grep total |cut -f1"
                P3 = subprocess.Popen(inputR1FileSizeCmd, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.STDOUT,  shell=True)
                #L = [P3.communicate(stdin), P3.communicate(stdout), P3.communicate(stderr)]
                P3.wait()
                inputR1FileSizeCmd = P3.communicate()[0].rstrip()

                outputR1FileSizeCmd = f"du -h -c  {outR1File} |grep total |cut -f 1"
                P4 = subprocess.Popen(outputR1FileSizeCmd, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.STDOUT,  shell=1)
                P4.wait()
                outputR1FileSizeCmd = P4.communicate()[0].rstrip()

                inputR2FileSizeCmd = f"du -h -c  {workdir}/{subdir}/*R2* |grep total |cut -f 1"
                P5 = subprocess.Popen(inputR2FileSizeCmd, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.STDOUT,  shell=1)
                P5.wait()
                inputR2FileSizeCmd = P5.communicate()[0].rstrip()

                outputR2FileSizeCmd = f"du -h -c  {outR2File} |grep total |cut -f 1"
                P6 = subprocess.Popen(outputR2FileSizeCmd, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.STDOUT,  shell=1)
                P6.wait()
                outputR2FileSizeCmd = P6.communicate()[0].rstrip()

                print (f"{subdir}\tinput R1 file size {inputR1FileSizeCmd}\t output R1 file size{outputR1FileSizeCmd}")
                print (f"{subdir}\tinput R2 file size {inputR2FileSizeCmd}\t output R2 file size{outputR2FileSizeCmd}")

if __name__ == "__main__":
        workdir = sys.argv[1]
        outdir = sys.argv[2]
        if len(sys.argv) != 3:
                sys.exit("Error: input!!! Usage: python step02.py inputDir outputDir")
        if not os.path.exists(workdir):
                sys.exit("Error: input!!! input directory does not exist")

        os.makedirs(outdir, exist_ok=True)
        workdir = os.path.abspath(workdir)
        outdir = os.path.abspath(outdir)
        make_dir(outdir)
        merge_same_sample_fastq(workdir, outdir)
        count_data_size(workdir, outdir)
