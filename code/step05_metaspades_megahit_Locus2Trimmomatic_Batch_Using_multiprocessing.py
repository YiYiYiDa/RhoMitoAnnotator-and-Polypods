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
import time
import argparse
from step04_Locus_Blast2_Trimmomatic_Batch_Using_multiprocessing_V2 import create_database_dict_and_get_L_locus as get_locus


####################################################
print("Usage:step05_metaspades_megahit_Locus2Trimmomatic_Batch_Using_multiprocessing.py /path/to/blastnFilteredReads /path/to/locus_database [-g/--gene <gene1> <gene2> ...] [-m/--megahit] [-st int(default=32)] [-mt int(default=32)]")
print("Step05:\tAssemble gene contigs from Blast reads via SPAdes and MEGAHIT(optional)\n##########")

"""
Purpose:
------------
    * De-novo assemble target-gene reads (filtered by step04 BLAST) using either SPAdes and MEGAHIT(optional).
    ** Output contigs are renamed as <sample>.<locus>.(spades|megahit).contigs.fasta for downstream mapping or annotation.

Usage:
------------
    python step05_metaspades_megahit_Locus2Trimmomatic_Batch_Using_multiprocessing.py /path/to/blastnFilteredReads /path/to/locus_database [-g/--gene <gene1> <gene2> ...] [-m/--megahit] [-st int(default=32)] [-mt int(default=32)]

Workflow:
------------
    1. Set input parameters;
    2. Create dict {gene:path_to_gene_fasta} from <databaseDir> (reuse step04 function);
    3. For each sample/<gene>/ folder:
        (1) check existence of R1/R2.paired.<locus>.fastq,
        (2) generate corresponding SPAdes and MEGAHIT(optional) command line;
    4. Run assemblies in parallel (4 worker processes by default);
    5. Rename contigs file to <sample>.<locus>.(spades|megahit).contigs.fasta and place in original assembly folder.

Notes:
------------
    * -m/--megahit flag is BOOLEAN; if present MEGAHIT will be executed, else only SPAdes.
    * The script will skip if contigs.fasta already exists and is non-empty, allowing safe re-runs.
    * Adjust ‘processes=’ in the multiprocessing.Pool call to match local CPU/memory resources.
    * SPAdes uses k-mer list 21,33,45...117 and 128 Gb RAM by default; MEGAHIT uses k-min 21, k-max 141, k-step 12 and 0.4 of the machine's total memory.
"""


####################################################
def args_parse():
	parser = argparse.ArgumentParser(description='This script is used to run spades.py to assemble gene from blastn-filtered reads',formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('workDir', help='path to blastn-filtered reads')
	parser.add_argument('databaseDir', help='path to databaseDir')
	parser.add_argument('-g', '--gene', type=str, default='all', nargs='+', metavar='str',
					   help='the name of gene to be assembled, support multi gene input(default: "all") :\n'
							'1.gene1;\n2.gene1 gene2 ... geneN;\n3.all')
	parser.add_argument('-m', '--megahit', help='use megahit to assemble gene', action='store_true')
	parser.add_argument('-st', type=int, default=32, metavar='int', help='the number of threads to use for spades.py (default: 32)')
	parser.add_argument('-mt', type=int, default=32, metavar='int', help='the number of threads to use for megahit (default: 32)')
	return parser.parse_args()

def main():
	args = args_parse()
	return args.workDir, args.databaseDir, args.gene, args.megahit, args.st, args.mt

def step1_run_spades_and_megahit(workdir, L_locus, megahitMaker, st, mt):
	if megahitMaker:
		if not shutil.which('megahit'):
			sys.exit("megahit not found in PATH")
	if not shutil.which('spades.py'):
		sys.exit("spades.py not found in PATH")
		
	L_Process_spades = []
	L_Process_megahit = []

	for sub_workdir in os.listdir(workdir):
		for locus in L_locus:
			locus_dir = workdir + "/" + sub_workdir + "/" + locus
			locus_R1_fqFile = locus_dir + "/" + sub_workdir + ".R1.paired." + locus + ".fastq"
			locus_R2_fqFile = locus_dir + "/" + sub_workdir + ".R2.paired." + locus + ".fastq"

			# obtain megahit command line
			if megahitMaker:
				outdir_megahit = locus_dir + "/megahit"
				if os.path.exists(outdir_megahit):
					if os.path.exists(outdir_megahit+f"/{sub_workdir}.{locus}.contigs.fa") and os.path.getsize(outdir_megahit+f"/{sub_workdir}.{locus}.contigs.fa"):
						print(f"-----megahit is OK\t{outdir_megahit}/contigs.fa")
					else:
						shutil.rmtree(outdir_megahit)
						cmd_megahit = f"megahit -o {outdir_megahit} --out-prefix {sub_workdir}.{locus} -t {mt} -m 0.4 --k-min 21 --k-max 141 --k-step 12 -1 {locus_R1_fqFile} -2 {locus_R2_fqFile}"
						L_Process_megahit.append(cmd_megahit)
				else:
					cmd_megahit = f"megahit -o {outdir_megahit} --out-prefix {sub_workdir}.{locus} -t {mt} -m 0.4 --k-min 21 --k-max 141 --k-step 12 -1 {locus_R1_fqFile} -2 {locus_R2_fqFile}"
					L_Process_megahit.append(cmd_megahit)


			# obtain spades command line
			outdir_spades = locus_dir + "/spades"
			os.makedirs(outdir_spades, exist_ok=True)
			if os.path.exists(outdir_spades+"/contigs.fasta") and os.path.getsize(outdir_spades+"/contigs.fasta"):
				print(f"-----spades is OK\t{outdir_spades}/contigs.fasta")
			else:
				cmd_spades = f"spades.py -1 {locus_R1_fqFile} -2 {locus_R2_fqFile} -t {st} --memory 128 -k 21,33,45,57,69,81,93,105,117 -o {outdir_spades}"
				L_Process_spades.append(cmd_spades)
	return L_Process_spades, L_Process_megahit

def run_spades_or_megahit(name):
	my_Cmd = name
	print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print(my_Cmd)
	subprocess.call(my_Cmd, shell=True)

def rename_contigs(workdir, L_locus, megahitMaker):
	for sub_workdir in os.listdir(workdir):
		for locus in L_locus:
			locus_dir = workdir + "/" + sub_workdir + "/" + locus

			# modify the output name of megahit
			if megahitMaker:
				outdir_megahit = locus_dir + "/megahit"
				output_megahit_name = outdir_megahit + "/" + sub_workdir + "." + locus + ".megahit.contigs.fasta"

				if os.path.exists(outdir_megahit):
					if os.path.exists(outdir_megahit + "/" + sub_workdir + "." + locus + ".contigs.fa") and os.path.getsize(outdir_megahit + "/" + sub_workdir + "." + locus + ".contigs.fa") != 0:
						shutil.copyfile(outdir_megahit + "/" + sub_workdir + "." + locus + ".contigs.fa",outdir_megahit + "/" + sub_workdir + "." + locus + ".megahit.contigs.fasta")
						print(f'rename\t{outdir_megahit}/{sub_workdir}.{locus}.contigs.fa ==> {outdir_megahit + "/" + sub_workdir + "." + locus + ".megahit.contigs.fasta"}')
					else:
						print("------ERROR: Does not exist or is empty",outdir_megahit + "/" + sub_workdir + "." + locus + ".megahit.contigs.fasta")
				else:
					print("------ERROR: Does not exist:", outdir_megahit)

			# modify the output name of spades
			outdir_spades = locus_dir + "/spades"
			output_spades_name = outdir_spades + "/" + sub_workdir + "." + locus + ".spades.contigs.fasta"

			if os.path.exists(outdir_spades):
				if os.path.exists(outdir_spades + "/contigs.fasta") and os.path.getsize(outdir_spades + "/contigs.fasta") != 0:
					shutil.copyfile(outdir_spades + "/contigs.fasta",outdir_spades + "/" + sub_workdir + "." + locus + ".spades.contigs.fasta")
					print(f'rename\t{outdir_spades}/contigs.fasta ==> {outdir_spades + "/" + sub_workdir + "." + locus + ".spades.contigs.fasta"}')
				else:
					print("------ERROR: Does not exist or is empty", output_spades_name)
			else:
				print("------ERROR: Does not exist:", outdir_spades)



if __name__ == '__main__':
	workdir,databaseDir, marker, megahitMaker, st, mt = main()
	if not os.path.exists(workdir):
		sys.exit("Error: input!!! input directory does not exist")

	D_queryFile, L_locus = get_locus(databaseDir, marker)
	L_Process_spades, L_Process_megahit = step1_run_spades_and_megahit(workdir, L_locus, megahitMaker, st, mt)

	if L_Process_spades:
		print('current process {0}'.format(os.getpid()))
		p = multiprocessing.Pool(processes=4)  # number of multiprocessing defalut:4
		for name in L_Process_spades:
			p.apply_async(run_spades_or_megahit, args=(name,))
		print('Waiting for all subprocesses done...')
		p.close()
		p.join()
		print('All processes done!')
	if L_Process_megahit:
		print('current process {0}'.format(os.getpid()))
		p = multiprocessing.Pool(processes=4)  # number of multiprocessing defalut:4
		for name in L_Process_megahit:
			p.apply_async(run_spades_or_megahit, args=(name,))
		print('Waiting for all subprocesses done...')
		p.close()
		p.join()
		print('All processes done!')



	rename_contigs(workdir, L_locus, megahitMaker)

