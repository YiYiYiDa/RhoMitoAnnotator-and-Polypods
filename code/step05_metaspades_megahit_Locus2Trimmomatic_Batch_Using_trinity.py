#encoding=utf-8
####################################################
import sys
import os
import shutil
from Bio import SeqIO
import fnmatch
import subprocess
import multiprocessing
import time

homepath = os.path.dirname(os.path.abspath(__file__))
modules_path = os.path.abspath(homepath+"/../modules/")
#miniconda_bin_path = "/home/mmshi/miniconda2/bin"
#miniconda_bin_path = "

blast_path = os.path.abspath(homepath+"/../modules/ncbi-blast-2.16.0+/bin")
megahit_path = os.path.abspath(homepath+"/../modules/megahit/bin")
db_path = os.path.abspath(homepath+"/../database")


# L_locus = ['atp1', 'atp4', 'atp6', 'atp8', 'atp9', 'ccmB', 'ccmC', 'ccmFc', 'ccmFn', 'cob', 'cox1', 'cox2', 'cox3', 'matR', 'mttB', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4L', 'nad5', 'nad6', 'nad7', 'nad9', 'rpl5', 'rpl10', 'rpl16', 'rps7', 'rps12', 'rps13', 'sdh4']

####################################################
print("Usage:%s Locus_Blast2Trimmomatic")

workdir = sys.argv[1] #Locus_Blast2Trimmomatic
#outdir = sys.argv[2]
marker = sys.argv[2].upper() #输入要组装的基因的名称

if marker == "NAD7":
	L_locus = ["nad7"]

elif marker == "NAD6":
	L_locus = ["nad6"]
elif marker == 'NAD2':
	L_locus = ["nad2"]
elif marker == 'MTTB':
	L_locus = ["mttB"]
elif marker == 'CCMFC':
	L_locus = ["ccmFc"]
elif marker == 'ALL':
	L_locus = ['atp1', 'atp4', 'atp6', 'atp8', 'atp9', 'ccmB', 'ccmC', 'ccmFc', 'ccmFn', 'cob', 'cox1', 'cox2', 'cox3', 'matR', 'mttB', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4L', 'nad5', 'nad6', 'nad7', 'nad9', 'rpl5', 'rpl10', 'rpl16', 'rps7', 'rps12', 'rps13', 'sdh4']
elif marker == 'MISS':
	L_locus = ['atp8', 'sdh4']
else:
	print("Error:", marker)
	sys.exit()

####################################################
# step 1: run Trinity and megahit

L_Process_trinity = []
L_Process_megahit = []

for sub_workdir in os.listdir(workdir):
	for locus in L_locus:
		locus_dir = workdir+"/"+sub_workdir+"/"+locus
		locus_R1_fqFile = locus_dir+"/"+sub_workdir+".R1.paired."+locus+".fastq"
		locus_R2_fqFile = locus_dir+"/"+sub_workdir+".R2.paired."+locus+".fastq"

		#obtain megahit command line
		outdir_megahit = locus_dir+"/megahit"
		if not os.path.exists(outdir_megahit):
			os.mkdir(outdir_megahit)
		if os.path.exists(outdir_megahit):
			shutil.rmtree(outdir_megahit)
			cmd_megahit = "%s/megahit -o %s --out-prefix %s -t 32 -m 0.4 --k-min 21 --k-max 141 --k-step 12 -1 %s -2 %s" % (megahit_path,outdir_megahit, sub_workdir+"."+locus, locus_R1_fqFile, locus_R2_fqFile)
			L_Process_megahit.append(cmd_megahit)

		"""
				if os.path.exists(outdir_megahit):
						if os.path.exists(outdir_megahit+"/"+sub_workdir+"."+locus+".contigs.fa") and os.path.getsize(outdir_megahit+"/"+sub_workdir+"."+locus+".contigs.fa"
) != 0:
								print "++++++\t", outdir_megahit+"/"+sub_workdir+"."+locus+".contigs.fa", " is OK!!!"
						else:
								shutil.rmtree(outdir_megahit)
								cmd_megahit = "megahit -o %s --out-prefix %s -t 16 -m 0.4 --k-min 21 --k-max 141 --k-step 12 -1 %s -2 %s" % (outdir_megahit, sub_workdir+"."+locus, locus_R1_fqFile, locus_R2_fqFile)
				else:
						#megahit -o matK.Saussurea_costus.megahit2 --out-prefix 111 -t 16 -m 0.4 --k-min 21 --k-max 141 --k-step 12 -1 matK.Saussurea_costus.R1.paired.fastq
						cmd_megahit = "megahit -o %s --out-prefix %s -t 16 -m 0.4 --k-min 21 --k-max 141 --k-step 12 -1 %s -2 %s" % (outdir_megahit,
sub_workdir+"."+locus, locus_R1_fqFile, locus_R2_fqFile)
						L_Process_megahit.append(cmd_megahit)
						#print cmd_megahit
		"""

		# obtain trinity command line
		outdir_trinity = locus_dir + "/trinity"
		if not os.path.exists(outdir_trinity):
			os.mkdir(outdir_trinity)
		if os.path.exists(outdir_trinity):
			shutil.rmtree(outdir_trinity)
			cmd_trinity = "Trinity --seqType fq --left %s --right %s --CPU 32 --max_memory 100G --output %s" % (locus_R1_fqFile, locus_R2_fqFile, outdir_trinity) # k-mer一般使用33、57、93（默认应该是33）[原始设置的k-mer数：21,33,45,57,69,81,93,105,117]
			L_Process_trinity.append(cmd_trinity)

cmd = ""

def run_task(name):
	my_Cmd = name
	print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	print(my_Cmd)

	subprocess.call(my_Cmd, shell=1)


if __name__ == '__main__':

	print('current process {0}'.format(os.getpid()))

	p = multiprocessing.Pool(processes=2) # number of multiprocessing defalut:4

	for name in L_Process_trinity:

		p.apply_async(run_task, args=(name,))

	print('Waiting for all subprocesses done...')

	p.close()
	p.join()
	print('All processes done!')

#sys.exit()

if __name__ == '__main__':

	print('current process {0}'.format(os.getpid()))

	p = multiprocessing.Pool(processes=2) #number of multiprocessing default:4

	for name in L_Process_megahit:

		p.apply_async(run_task, args=(name,))

	print('Waiting for all subprocesses done...')

	p.close()
	p.join()
	print('All processes done!')

for sub_workdir in os.listdir(workdir):
		for locus in L_locus:
				locus_dir = workdir+"/"+sub_workdir+"/"+locus
				#modify the output name of megahit
				outdir_megahit = locus_dir+"/megahit"
				output_megahit_name = outdir_megahit+"/"+sub_workdir+"."+locus+".megahit.contigs.fasta"

				if os.path.exists(outdir_megahit):
						if os.path.exists(outdir_megahit+"/"+sub_workdir+"."+locus+".contigs.fa") and os.path.getsize(outdir_megahit+"/"+sub_workdir+"."+locus+".contigs.fa") != 0:
							shutil.copyfile(outdir_megahit+"/"+sub_workdir+"."+locus+".contigs.fa", outdir_megahit+"/"+sub_workdir+"."+locus+".megahit.contigs.fasta")
						else:
							print("------ERROR: Does not exist or is empty", outdir_megahit+"/"+sub_workdir+"."+locus+".megahit.contigs.fasta")
				else:
					print("------ERROR: Does not exist:", outdir_megahit)

				#modify the output name of trinity

				outdir_trinity = locus_dir+"/trinity"
				output_trinity_name = outdir_trinity+"/"+"Trinity.fasta"

				if os.path.exists(outdir_trinity):

						if os.path.exists(outdir_trinity+"/Trinity.fasta") and os.path.getsize(outdir_trinity+"/Trinity.fasta") !=0:
								shutil.copyfile(outdir_trinity+"/Trinity.fasta", outdir_trinity+"/"+sub_workdir+"."+locus+".trinity.contigs.fasta")

						else:
								print("------ERROR: Does not exist or is empty", output_trinity_name)
				else:
					print("------ERROR: Does not exist:", outdir_trinity)


"""
This CAP3 program was temporarily closed #20200424

print "=======================\trun CAP3"

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

L_Process_Cap3 = []
#L_locus = ["COI", "ITS", "ITS2", "matK", "rbcL"]

for sub_workdir in os.listdir(workdir):

		for locus in L_locus:
				locus_dir = workdir+"/"+sub_workdir+"/"+locus

		CAP3_outdir = locus_dir+"/CAP3"
				if not os.path.exists(CAP3_outdir):
						os.mkdir(CAP3_outdir)

		#print "++++++++ I am here"
		if not os.path.exists(CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta") or os.path.getsize(CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta") ==0:

						outFasta = open(CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta", "w")
						outQual = open(CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta.qual", "w")

			#print CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta"
			#print "++++++++----- I am here"

						notCombined_1File = locus_dir+"/"+sub_workdir+".R1.paired."+locus+".fastq"
						notCombined_2File = locus_dir+"/"+sub_workdir+".R2.paired."+locus+".fastq"

						for record in SeqIO.parse(notCombined_1File, "fastq"):
								newRecord = SeqRecord(id = record.id+".x1", seq=record.seq, letter_annotations = record.letter_annotations, description="")
								outFasta.write(newRecord.format("fasta"))
								outQual.write(newRecord.format("qual"))
			#print "++++++++----- I am here"

						for record in SeqIO.parse(notCombined_2File, "fastq"):
								newRecord = SeqRecord(id = record.id+".y1", seq=record.seq, letter_annotations = record.letter_annotations, description="")
								outFasta.write(newRecord.format("fasta"))
								outQual.write(newRecord.format("qual"))

			#print "++++++++----- I am here"
						outFasta.close()
						outQual.close()

		if not os.path.exists(CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta.con") or os.path.getsize(CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta.con") ==0:
						cmd = "%s/CAP3/formcon %s 1 500" % (modules_path, CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta")
						P1 = subprocess.Popen(cmd, shell=1)
						P1.wait()

		#if not os.path.exists(CAP3_outdir+"/"+sub_workdir+"."+locus+".CAP3.contigs.fasta") or os.path.getsize(CAP3_outdir+"/"+sub_workdir+"."+locus+".CAP3.contigs.fasta"):
		#.cap.contigs
		#GCH01.COI.fasta.cap.contigs

		if not os.path.exists(CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta.cap.contigs") or os.path.getsize(CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta.cap.contigs")==0:

			cap3_cmd = "%s/CAP3/cap3 %s -o 25 -p 90 >%s" % (modules_path, CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta", locus_dir+"/"+sub_workdir+"."+locus+".fasta.cap.log")
					L_Process_Cap3.append(cap3_cmd)
		else:
			print "++++++++++++\t", CAP3_outdir+"/"+sub_workdir+"."+locus+".fasta.cap.contigs", "is OK!!!"

def run_task(name):
		my_Cmd = name
		print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print my_Cmd

		subprocess.call(my_Cmd, shell=1)


if __name__ == '__main__':

	print('current process {0}'.format(os.getpid()))

	p = multiprocessing.Pool(processes=10)

	for name in L_Process_Cap3:

		p.apply_async(run_task, args=(name,))

	print('Waiting for all subprocesses done...')

	p.close()
	p.join()
	print('All processes done!')


for sub_workdir in os.listdir(workdir):

		for locus in L_locus:

				locus_dir = workdir+"/"+sub_workdir+"/"+locus

				#modify the output name of megahit
				outdir_Cap3 = locus_dir+"/CAP3"
				output_Cap3_FastaName = outdir_Cap3+"/"+sub_workdir+"."+locus+".fasta.cap.contigs"
		output_Cap3_QualName = outdir_Cap3+"/"+sub_workdir+"."+locus+".fasta.cap.contigs.qual"

		if os.path.exists(outdir_Cap3):
			if os.path.exists(output_Cap3_FastaName) and os.path.getsize(output_Cap3_FastaName) !=0:
				shutil.copyfile(output_Cap3_FastaName, outdir_Cap3+"/"+sub_workdir+"."+locus+".Cap3.contigs.fasta")
				shutil.copyfile(output_Cap3_QualName, outdir_Cap3+"/"+sub_workdir+"."+locus+".Cap3.contigs.fasta.qual")
			else:
				print ("------ERROR: Does not exist or is empty", output_Cap3_FastaName)
		else:
			print ("------ERROR: Does not exist:", outdir_Cap3)
"""				
