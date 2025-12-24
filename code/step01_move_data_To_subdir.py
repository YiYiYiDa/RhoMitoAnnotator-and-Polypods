# Date: 20251128
# Author: Linchun Shi, Yanda Zhu
# Copyright: (c) 2025 Linchun Shi, Yanda Zhu. All rights reserved.


#####################################################
import sys
import os
import subprocess
import multiprocessing

#####################################################
print ("Usage: DataDir(input)")
print ("Step01:\tmake subdir\n##########")

# workdir = sys.argv[1] # Data dir

'''
Purpose:
------------
	* Automatically create one sub-directory per biological sample and move all FASTQ files that belong to the same sample into its corresponding folder.
	** The script recognises samples solely by the characters that precede the first underscore (“_”) in the file name.  
	** Any file that does not carry one of the accepted extensions (fastq, fastq.gz, fq, fq.gz) is reported but left untouched.

Usage:
------------
	python step01_organize_fastq.py /absolute/or/relative/path/to/DataDir

Workflow:
------------
	1. Validate the user-supplied directory.
	2. Walk through the directory tree (only the top level) and collect all files whose names end with an accepted extension.
	3. Extract the sample prefix (text before the first “_”) and build aunique list of sample identifiers.
	4. Create one sub-directory for every sample identifier.
	5. Move all files that start with “<sample>” into the corresponding sub-directory by issuing a single shell “mv” command per sample.

Notes:
------------
	* The directory walk is intentionally on the top level (dirnames[:] = []), i.e. only the content directly below DataDir is scanned.
	* purposed data file must end with [fastq.gz, fastq, fq.gz, fq]
	* The script uses subprocess.Popen(shell=True) to perform the move operation. Ensure that file names do not contain shell-meta-characters.
'''

###################################################

def get_sample_name(workdir):
	L_sampleID = []
	N_sampleID = []
	valid_suffix = ('fastq.gz', 'fastq', 'fq.gz', 'fq')
	for dirpath, dirnames, filenames in os.walk(workdir):
		dirnames[:] = []
		for workfile in filenames:
			if workfile.lower().endswith(valid_suffix):
				handle = workfile.split("_")[0]
				if handle not in L_sampleID:
					L_sampleID.append(handle)
			else:
				N_sampleID.append(workfile)
	if not L_sampleID:
		print('Not have the file meets requirement:fastq.gz, fastq, fq.gz, fq')
	if N_sampleID:
		N_sampleID_info = "\n".join(N_sampleID)
		print(f'There are some files not meets requirement:\n{N_sampleID_info}')
	if L_sampleID:
		return L_sampleID


def move_same_sample_data_to_subdir(L_sampleID, workdir):
	L_sampleID.sort()
	for it in L_sampleID:
		# print(it)
		subdir = workdir + "/" + it
		if not os.path.exists(subdir):
			os.mkdir(subdir)
		cmd = f"mv {workdir}/{it}_* {subdir}"
		p = subprocess.Popen(cmd, shell=True)
		p.wait()
	print("all data have been moved to subdir")

if __name__ == '__main__':
	workdir = sys.argv[1]
	if len(sys.argv) != 2:
		sys.exit("Error: input!!! Usage: python step01.py DataDir(input)")
	if not workdir:
		sys.exit("Error: input!!! Don't have that dir")

	workdir = os.path.abspath(workdir)
	L_sampleID = get_sample_name(workdir)
	if L_sampleID:
		move_same_sample_data_to_subdir(L_sampleID, workdir)



