# Date: 20251128
# Author: Linchun Shi, Yanda Zhu
# Copyright: (c) 2025 Linchun Shi, Yanda Zhu. All rights reserved.

####################################################
import sys
import os
import shutil
from Bio import SeqIO
import fnmatch
import subprocess
import multiprocessing
import argparse


#####################################################
print ("Usage: program dataMerge(input) Trimmomatic(output) -t/--threads N (optional, default: 6)")
print("Step03:\trun trimmomatic on dataMerge to control data quality\n##########")

"""
Purpose: 
------------
    !!! If you have already run Trimmomatic or controlled quality using other tools, you can skip this step.
    * Perform quality controlling and adapter trimming on paired-end Illumina reads using Trimmomatic.
    ** Process every sample subdirectory generated in the previous merge step (step02).

Usage:
------------
    python step03_run_Trimmomatic.py /path/to/DataMergeDir /path/to/TrimmomaticOutDir -t/--thread <threads>

Workflow:
------------
    1. Set input parameters;
    2. Verify Java and Trimmomatic availability,  build Trimmomatic jar and adapter path; 
    3. For each sample sub-folder in <DataMergeDir>:
        (1) create the same sub-folder in <TrimmomaticOutDir>,
        (2) build Trimmomatic PE command with adapters, quality sliding window, leading/trailing crop, minimum length, etc.,
        (3) append the full command to a job list.
    4. Run all commands in parallel (5 worker processes by default) and wait until completion, print each command before execution for logging/debugging.

Notes:
------------
    * step03_run_Trimmomatic.py -h/--help for help.
    * Default Trimmomatic verison is 0.39, if you want to use other versions, replace 'trimmomatic-0.39.jar' in this script.
    * Output files are written as <sample>.R1/R2.paired.fastq and <sample>.R1/R2.unpaired.fastq alongside a <sample>.summary log.
    * The script skips re-merging if the target files already exist and are non-empty, allowing safe re-runs.
    * Adjust ‘processes=5’ in the multiprocessing.Pool call to match local CPU/memory resources.
"""
#####################################################
def args_parse():
    prog = os.path.basename(sys.argv[0])
    parser = argparse.ArgumentParser(description='run trimmomatic on dataMerge to control data quality',
                                    usage=f'{prog} <dataMerge> <trimmomaticDir> [-t/--threads N (default: 6)]\n'
                                          f'{prog} -h/--help\n'
                                    )
    parser.add_argument('dataMerge', help='path to DataMergeDir')
    parser.add_argument('trimmomaticDir', help='path to TrimmomaticOutDir')
    parser.add_argument('-t', '--threads', type=int, default=6, metavar='N', help='number of threads to use for Trimmomatic (default: 6)')
    return parser.parse_args()

def main():
    args = args_parse()
    if args.threads <= 0:
        print("Error: threads must be a positive integer")
        sys.exit(1)
    return args.dataMerge, args.trimmomaticDir, int(args.threads)


def run_trimmomatic_on_dataMerge(workdir, outdir, threadsTrimmomatic):
    if not shutil.which('trimmomatic-0.39.jar'):
        print("Trimmomatic-0.39.jar not found in PATH")
    else:
        trimmomaticPath = '/'.join(shutil.which('trimmomatic-0.39.jar').split('/')[:-1])
        adapterDir = os.path.join(trimmomaticPath, "adapters")
        if not os.path.exists(adapterDir):
            print(f"Adapter directory not in : {adapterDir}")
    if not shutil.which('java'):
        print("Java not found in PATH")

    L_cmd = []
    for sub_workdir in os.listdir(workdir):
            if os.path.isdir(workdir+"/"+sub_workdir):
                os.makedirs(outdir+'/'+sub_workdir, exist_ok=True)
                PE1 = sub_workdir+".R1.fastq.gz"
                PE2 = sub_workdir+".R2.fastq.gz"

                PE1_paired = outdir+"/"+sub_workdir+"/"+PE1.replace(".fastq", ".paired.fastq")
                PE1_unpaired = outdir+"/"+sub_workdir+"/"+PE1.replace(".fastq", ".unpaired.fastq")

                PE2_paired = outdir+"/"+sub_workdir+"/"+PE2.replace(".fastq", ".paired.fastq")
                PE2_unpaired = outdir+"/"+sub_workdir+"/"+PE2.replace(".fastq", ".unpaired.fastq")

                PE1_input = workdir+"/"+sub_workdir+"/"+PE1
                PE2_input = workdir+"/"+sub_workdir+"/"+PE2

                summaryFile = outdir+"/"+sub_workdir+"/"+sub_workdir+".summary"

                if not os.path.exists(PE1_paired) or os.path.getsize(PE1_paired) == 0:
                    cmd = f"java -jar {trimmomaticPath}/trimmomatic-0.39.jar PE -threads {threadsTrimmomatic} -phred33 -summary {summaryFile} {PE1_input} {PE2_input} {PE1_paired} {PE1_unpaired} {PE2_paired} {PE2_unpaired} ILLUMINACLIP:{adapterDir}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
                    L_cmd.append(cmd)
                else:
                    print(f"File {PE1_paired} already exists, skipping Trimmomatic PE for {sub_workdir}")
    return L_cmd

def run_task(name):
        my_Cmd = name
        print(my_Cmd)
        print ("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        subprocess.call(my_Cmd, shell=True)


if __name__ == '__main__':
    workdir, outdir, threadsTrimmomatic = main()
    if not os.path.exists(workdir):
        sys.exit("Error: input!!! input directory does not exist")

    os.makedirs(outdir, exist_ok=True)
    workdir = os.path.abspath(workdir)
    outdir = os.path.abspath(outdir)
    L_cmd = run_trimmomatic_on_dataMerge(workdir, outdir, threadsTrimmomatic)

    if L_cmd:
        print(f'current process {os.getpid()}')
        p = multiprocessing.Pool(processes=5)
        for name in L_cmd:
            p.apply_async(run_task, args=(name,))
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All processes done!')


