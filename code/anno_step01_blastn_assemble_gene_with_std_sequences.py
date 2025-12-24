# Date: 20251128
# Author: Linchun Shi, Yanda Zhu
# Copyright: (c) 2025 Linchun Shi, Yanda Zhu. All rights reserved.

####################################################

import sys
import os
import fnmatch
import subprocess
import argparse
from Bio import SeqIO
from distutils.core import extension_keywords

####################################################

# 需要修改的地方：（1）原报错文件中同nad1、nad2的插入不报错【所有正常blast比对中的gap、insert不报错】？（2）加入非经典起始密码子、终止密码子、碱基非3的倍数报错？（3）nad6只要注释长度在600以上就不报错？
def args_parse():
    parser = argparse.ArgumentParser(description='This script is used to run blastn to create annoation files')
    parser.add_argument('assembleDir', help='path to assemble results')
    parser.add_argument('databaseDir', help='path to databaseDir')
    parser.add_argument('outDir', help='path to ouputDir')
    return parser.parse_args()

def main():
    args = args_parse()
    return args.assembleDir, args.databaseDir, args.outDir

def find_exon(txt, first_end_site, gene, multiple_exon_annotation_file, n=0):
    with open(multiple_exon_annotation_file, 'a') as ma:
        exonLengthDict = {
            'ccmfc':[767, 547],
            'cox2':[382, 318, 83],
            'nad1':[385, 83, 192, 59, 259],
            'nad2':[153, 392, 161, 573, 188],
            'nad4':[461, 515, 423, 89],
            'nad5':[230, 1216, 22, 395, 150],
            'nad7':[143, 69, 467, 244, 262]
        }
        n += 1

        if n + 1 <= len(exonLengthDict[gene]):
            exon_length = int(exonLengthDict[gene][n])

            exon_marker = False
            if exon_length == 22 and gene == 'nad5':
                mismatch = 3
                for line in txt:
                    if not line.strip():
                        continue
                    start = int(line.split('\t')[6])
                    end = int(line.split('\t')[7])
                    if int(end-start+1) in range(20, 24):
                        if start in range(first_end_site - 2, first_end_site + 2) and int(line.strip().split('\t')[4]) < mismatch:
                            site_start = start
                            site_end = end
                            readsName = line.split('\t')[1]
                            if int(line.split('\t')[8]) < int(line.split('\t')[9]):
                                exon_start = int(line.split('\t')[8]) + int(first_end_site - start) + 1
                                exon_end = int(line.split('\t')[9])
                                length = int(abs(exon_end - exon_start)) + 1
                            else:
                                exon_start = int(line.split('\t')[8]) - int(first_end_site - start) - 1
                                exon_end = int(line.split('\t')[9])
                                length = int(abs(exon_end - exon_start)) + 1
                            ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\n')
                            print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
                            exon_marker = True
                            break
                if not exon_marker:
                    for line in txt:
                        if not line.strip():
                            continue
                        start = int(line.split('\t')[6])
                        end = int(line.split('\t')[7])
                        if start > first_end_site:
                            if exon_length - 2 <= end - start <= exon_length + 2:
                                site_start = start
                                site_end = end
                                readsName = line.split('\t')[1]
                                exon_start = int(line.split('\t')[8])
                                exon_end = int(line.split('\t')[9])
                                length = int(abs(exon_end - exon_start)) + 1
                                ma.write(
                                    f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n} and {n + 1} has gap {start - first_end_site}\n')
                                print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
                                exon_marker = True
                if not exon_marker:
                    ma.write(f'{gene} blastn out has no exon {n + 1} and may annotation wrong \n')
                    return print(f'{gene} blastn out has no exon {n + 1} and may annotation wrong')
                if exon_marker:
                    return find_exon(txt, site_end, gene, multiple_exon_annotation_file, n)

            else:
                for line in txt:
                    if not line.strip():
                        continue
                    start = int(line.split('\t')[6])
                    end = int(line.split('\t')[7])
                    if start == first_end_site + 1:
                        if exon_length - 20 <= end - start <= exon_length + 20:

                            not_end_maker = False
                            if n + 1 == len(exonLengthDict[gene]):
                                std_length = 0
                                for i in exonLengthDict[gene]:
                                    std_length += int(i)
                                if end != std_length:
                                    not_end_maker = True
                                    end_gap = std_length - end

                            site_start = start
                            site_end = end
                            readsName = line.split('\t')[1]
                            gap = int(line.split('\t')[5])
                            exon_start = int(line.split('\t')[8])
                            exon_end = int(line.split('\t')[9])
                            length = int(abs(exon_end - exon_start)) + 1
                            if not_end_maker:
                                if gap == 0:
                                    ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has gap {end_gap} with end\n')
                                    print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has gap {end_gap} with end')
                                else:
                                    ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n + 1} has gap {end_gap} with end(2)exon {n + 1} has insert {gap} sites\n')
                                    print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n + 1} has gap {end_gap} with end(2)exon {n + 1} has insert {gap} sites')
                            else:
                                if gap == 0:
                                    ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\n')
                                    print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
                                else:
                                    ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has insert {gap} sites\n')
                                    print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has insert {gap} sites')
                            exon_marker = True
                            break
                if not exon_marker:
                    for line in txt:
                        if not line.strip():
                            continue
                        gap = int(line.split('\t')[5])
                        start = int(line.split('\t')[6])
                        end = int(line.split('\t')[7])
                        if start <= first_end_site:
                            if exon_length - 20 <= end - start <= exon_length + 20 and (first_end_site - start <= 20):
                                if gene == 'nad1' and gap:
                                    align_length = int(line.split('\t')[3])
                                    exon_start = int(line.split('\t')[8])
                                    exon_end = int(line.split('\t')[9])
                                    if exon_end - exon_start + 1 + gap == align_length:
                                        site_start = first_end_site + 1
                                        site_end = end
                                        readsName = line.split('\t')[1]
                                        if int(line.split('\t')[8]) < int(line.split('\t')[9]):
                                            exon_start = int(line.split('\t')[8]) + int(first_end_site - start) + 1 - gap
                                            exon_end = int(line.split('\t')[9])
                                            length = int(abs(exon_end - exon_start)) + 1
                                        else:
                                            exon_start = int(line.split('\t')[8]) - int(first_end_site - start) - 1 + gap
                                            exon_end = int(line.split('\t')[9])
                                            length = int(abs(exon_end - exon_start)) + 1
                                        ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\n')
                                        print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
                                        exon_marker = True
                                        break
                                else:

                                    not_end_maker = False
                                    if n + 1 == len(exonLengthDict[gene]):
                                        std_length = 0
                                        for i in exonLengthDict[gene]:
                                            std_length += int(i)
                                        if end != std_length:
                                            not_end_maker = True
                                            end_gap = std_length - end

                                    site_start = first_end_site + 1
                                    site_end = end
                                    readsName = line.split('\t')[1]
                                    if int(line.split('\t')[8]) < int(line.split('\t')[9]):
                                        gap = int(line.split('\t')[5])
                                        exon_start = int(line.split('\t')[8]) + int(first_end_site - start) + 1
                                        exon_end = int(line.split('\t')[9])
                                        length = int(abs(exon_end - exon_start)) + 1
                                    else:
                                        gap = int(line.split('\t')[5])
                                        exon_start = int(line.split('\t')[8]) - int(first_end_site - start) - 1
                                        exon_end = int(line.split('\t')[9])
                                        length = int(abs(exon_end - exon_start)) + 1
                                    if not_end_maker:
                                        if gap == 0:
                                            ma.write( f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has gap {end_gap} with end\n')
                                            print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has gap {end_gap} with end')
                                        else:
                                            ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n + 1} has gap {end_gap} with end(2)exon {n + 1} has insert {gap} sites\n')
                                            print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n + 1} has gap {end_gap} with end(2)exon {n + 1} has insert {gap} sites')
                                    else:
                                        if gap == 0:
                                            ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\n')
                                            print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
                                        else:
                                            ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has insert {gap} sites\n')
                                            print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has insert {gap} sites')
                                    exon_marker = True
                                    break
                if not exon_marker:
                    for line in txt:
                        if not line.strip():
                            continue
                        start = int(line.split('\t')[6])
                        end = int(line.split('\t')[7])
                        if start > first_end_site:
                            if exon_length - 20 <= end - start <= exon_length + 20:

                                not_end_maker = False
                                if n + 1 == len(exonLengthDict[gene]):
                                    std_length = 0
                                    for i in exonLengthDict[gene]:
                                        std_length += int(i)
                                    if end != std_length:
                                        not_end_maker = True
                                        end_gap = std_length - end

                                site_start = start
                                site_end = end
                                readsName = line.split('\t')[1]
                                gap = int(line.split('\t')[5])
                                exon_start = int(line.split('\t')[8])
                                exon_end = int(line.split('\t')[9])
                                length = int(abs(exon_end - exon_start)) + 1
                                if not_end_maker:
                                    if gap == 0:
                                        ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n + 1} has gap {start - first_end_site}(2)exon {n + 1} has gap {end_gap} with end\n')
                                        print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n + 1} has gap {start - first_end_site}(2)exon {n + 1} has gap {end_gap} with end')
                                    else:
                                        ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n + 1} has gap {start - first_end_site}(2)exon {n + 1} has gap {end_gap} with end(3)exon {n + 1} has insert {gap} sites\n')
                                        print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n + 1} has gap {start - first_end_site}(2)exon {n + 1} has gap {end_gap} with end(3)exon {n + 1} has insert {gap} sites')
                                else:
                                    if gap == 0:
                                        ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n} and {n+1} has gap {start - first_end_site}\n')
                                        print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n} and {n+1} has gap {start - first_end_site}')
                                    else:
                                        ma.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n+1} has gap {start - first_end_site}(2)exon {n + 1} has insert {gap} sites\n')
                                        print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n+1} has gap {start - first_end_site}(2)exon {n + 1} has insert {gap} sites')
                                exon_marker = True
                                break
                if not exon_marker:
                    ma.write(f'{gene} blastn out has no exon {n+1}\n')
                    return print(f'{gene} blastn out has no exon {n+1}')
                if exon_marker:
                    return find_exon(txt, site_end, gene, multiple_exon_annotation_file, n)
        else:
            return None


def run_blastn_for_assemble_results(assembleDir,databaseDir):
    for sample in os.listdir(assembleDir):
        if os.path.isdir(os.path.join(assembleDir, sample)):
            sampleDir = os.path.join(assembleDir, sample)

            blastOutDir = os.path.join(outDir, sample, 'blastOut')
            blastDatabaseDir = os.path.join(outDir, sample, 'database')

            os.makedirs(blastOutDir, exist_ok=True)
            os.makedirs(blastDatabaseDir, exist_ok=True)

            for gene in os.listdir(sampleDir):
                ## find sequence assembled by spades
                spadesDir = os.path.join(sampleDir, gene, 'spades')
                matchAssFile =  fnmatch.filter(os.listdir(spadesDir), '*.contigs.fasta')
                if not matchAssFile:
                    continue
                assemble_sequence = os.path.join(spadesDir, matchAssFile[0])

                ## find reference sequence from database
                for file in fnmatch.filter(os.listdir(databaseDir), '*.fasta'):
                    if gene.lower() == file.lower().split('_')[0]:
                        std_sequence_file = os.path.join(databaseDir, file)

                ## blastn iterly
                if os.path.exists(std_sequence_file):
                    if gene.lower() == 'nad5':
                        if not os.path.exists(os.path.join(blastDatabaseDir, f'[{matchAssFile[0]}].ndb')):
                            makeblastdb_cmd = f'makeblastdb -in {assemble_sequence} -dbtype nucl -out {blastDatabaseDir}/{matchAssFile}'
                            subprocess.call(makeblastdb_cmd, shell=True)
                        else:
                            print(f"-----{gene} database is OK\t{blastDatabaseDir}/{matchAssFile}")
                        if not os.path.exists(os.path.join(blastOutDir, f'std_{gene}_with_assemble.out')):
                            blastn_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 6 -word_size 7 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble.out'
                            subprocess.call(blastn_cmd, shell=True)
                        else:
                            print(f"-----{gene} blastout is OK\t{blastOutDir}/std_{gene}_with_assemble.out")
                        if not os.path.exists(os.path.join(blastOutDir, f'std_{gene}_with_assemble_0.out')):
                            blastn0_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 0 -word_size 7 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble_0.out'
                            subprocess.call(blastn0_cmd, shell=True)
                        else:
                            print(f"-----{gene} blastout0 is OK\t{blastOutDir}/std_{gene}_with_assemble_0.out")

                    else:
                        if not os.path.exists(os.path.join(blastDatabaseDir, f'[{matchAssFile[0]}].ndb')):
                            makeblastdb_cmd = f'makeblastdb -in {assemble_sequence} -dbtype nucl -out {blastDatabaseDir}/{matchAssFile}'
                            subprocess.call(makeblastdb_cmd, shell=True)
                        else:
                            print(f"-----{gene} database is OK\t{blastDatabaseDir}/{matchAssFile}")
                        if not os.path.exists(os.path.join(blastOutDir, f'std_{gene}_with_assemble.out')):
                            blastn_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 6 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble.out'
                            subprocess.call(blastn_cmd, shell=True)
                        else:
                            print(f"-----{gene} blastout is OK\t{blastOutDir}/std_{gene}_with_assemble.out")
                        if not os.path.exists(os.path.join(blastOutDir, f'std_{gene}_with_assemble_0.out')):
                            blastn0_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 0 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble_0.out'
                            subprocess.call(blastn0_cmd, shell=True)
                        else:
                            print(f"-----{gene} blastout0 is OK\t{blastOutDir}/std_{gene}_with_assemble_0.out")

def annotation_for_single_exon_gene(assembleDir, outDir, singleExonGene):
    for sample in os.listdir(outDir):
        if os.path.isdir(os.path.join(outDir, sample)):

            ## get blastout file path
            blastOutDir = os.path.join(outDir, sample, 'blastOut')
            ## get annotation file path
            annotationDir = os.path.join(outDir, sample, 'annotation')
            os.makedirs(annotationDir, exist_ok=True)
            singe_exon_annotation_file = os.path.join(annotationDir, 'single_exon_annotation.txt')
            with open(singe_exon_annotation_file, 'w') as sa:
                sa.write('gene\treadsName\tassemble_sequence_path\tlen_max\ts_start\ts_end\tinfo\n')

                ## check unassemble single exon gene
                for gene in singleExonGene:
                    geneMaker = False
                    for file in os.listdir(blastOutDir):
                        if gene.lower() == file.split('_')[1].lower():
                            geneMaker = True
                    if not geneMaker:
                        sa.write(f'don not have {gene}\n')
                        print(f'don not have {gene}')

                ## !!! annotation from blastout
                for outFIle in fnmatch.filter(os.listdir(blastOutDir), '*assemble.out'):
                    gene = outFIle.split('_')[1]

                    single_gene_marker = False
                    for i in singleExonGene:
                        if gene.lower() == i.lower():
                            single_gene_marker = True
                    if single_gene_marker:

                        len_max = 0
                        singleBlastInfo = ''
                        with open(os.path.join(blastOutDir, outFIle), 'r') as blastFile:
                            for line in blastFile:
                                singleBlastInfo += line
                                length = int(line.split('\t')[3])
                                mismatch = int(line.split('\t')[4])
                                if length > len_max:
                                    len_max = length
                                    queryName = line.split('\t')[0]
                                    readsName = line.split('\t')[1]
                                    mismatch = int(line.split('\t')[4])
                                    insert = int(line.split('\t')[5])
                                    q_start = int(line.split('\t')[6])
                                    q_end = int(line.split('\t')[7])
                                    s_start = int(line.split('\t')[8])
                                    s_end = int(line.split('\t')[9])

                            if s_start and s_end and len_max:
                                ### get path to assemble sequence
                                spadesDir = os.path.join(assembleDir, sample, gene, 'spades')
                                matchAssFile = fnmatch.filter(os.listdir(spadesDir), '*.contigs.fasta')
                                assemble_sequence = os.path.join(spadesDir, matchAssFile[0])

                                #### get std sequence length
                                stdSequenceLen = {}
                                for file in os.listdir(databaseDir):
                                    if gene.lower() == file.split('_')[0].lower():
                                        matchStdFile = os.path.join(databaseDir, file)
                                with open(matchStdFile, 'r') as std:
                                    for record in SeqIO.parse(std, "fasta"):
                                        stdSequenceLen[record.id] = len(record)

                                s_maker = False
                                e_maker = False
                                ### write annotation file
                                if q_start ==1 and q_end == stdSequenceLen[queryName]:
                                    sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\n')
                                    print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}')
                                else:
                                    if q_start !=1:
                                        s_maker = False
                                        for line in singleBlastInfo.split('\n'):
                                            if line:
                                                if line.strip().split('\t')[0] == queryName and int(line.strip().split('\t')[6]) == 1:
                                                    s_maker = True
                                    if q_end != stdSequenceLen[queryName]:
                                        e_maker = False
                                        for line in singleBlastInfo.split('\n'):
                                            if line:
                                                if line.strip().split('\t')[0] == queryName and int(line.strip().split('\t')[7]) == stdSequenceLen[queryName]:
                                                    e_maker = True

                                    if e_maker and s_maker:
                                        if gene.lower() == 'nad6' and len_max > 619:
                                            sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tend site lose {stdSequenceLen[queryName] - q_end} and may annotation wrong\n')
                                            print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tend site lose {stdSequenceLen[queryName] - q_end} and may annotation wrong')
                                        else:
                                            sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {stdSequenceLen[queryName] - q_end} and may annotation wrong(2)start site lose {q_start - 1} sites and may annotation wrong\n')
                                            print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {stdSequenceLen[queryName] - q_end} and may annotation wrong(2)start site lose {q_start - 1} sites and may annotation wrong')
                                    elif not e_maker and s_maker:
                                        sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tstart site lose {q_start - 1} sites and may annotation wrong\n')
                                        print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tstart site lose {q_start - 1} sites and may annotation wrong')
                                    elif e_maker and not s_maker:
                                        if gene.lower() == 'nad6' and len_max > 619:
                                            sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\n')
                                            print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}')
                                        else:
                                            sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tend site lose {stdSequenceLen[queryName] - q_end} and may annotation wrong\n')
                                            print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tend site lose {stdSequenceLen[queryName] - q_end} and may annotation wrong')
                                    else:
                                        sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\n')
                                        print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}')

def annotation_for_multi_exon_gene(assembleDir, outDir, multiExonGene, exonLengthDict):
    for sample in os.listdir(outDir):
        if os.path.isdir(os.path.join(outDir, sample)):

            blastOutDir = os.path.join(outDir, sample, 'blastOut')
            annotationDir = os.path.join(outDir, sample, 'annotation')
            os.makedirs(annotationDir, exist_ok=True)
            multiple_exon_annotation_file = os.path.join(annotationDir, 'multiple_exon_annotation.txt')
            with open(multiple_exon_annotation_file, 'w') as ma:
                ma.write('gene\treadsName\tassemble_sequence_path\texon_length\texon_start\texon_end\n')

                ## check unassemble multi exon gene
            for gene in multiExonGene:
                geneMaker = False
                for file in os.listdir(blastOutDir):
                    if gene.lower() == file.split('_')[1].lower():
                        geneMaker = True
                if not geneMaker:
                    sa.write(f'don not have {gene}\n')
                    print(f'don not have {gene}')

            ## !!! annotation from blastout
            for outFIle in fnmatch.filter(os.listdir(blastOutDir), '*assemble.out'):
                gene = outFIle.split('_')[1]

                multi_gene_marker = False
                for i in multiExonGene:
                    if gene.lower() == i.lower():
                        multi_gene_marker = True
                if multi_gene_marker:

            # for outFIle in fnmatch.filter(os.listdir(blastOutDir), '*assemble.out'):
            #     gene = outFIle.split('_')[1]

                    #### get path to assemble sequence
                    spadesDir = os.path.join(assembleDir, sample, gene, 'spades')
                    matchAssFile = fnmatch.filter(os.listdir(spadesDir), '*.contigs.fasta')
                    assemble_sequence = os.path.join(spadesDir, matchAssFile[0])

                    #### get std sequence length
                    stdSequenceLen = {}
                    for file in os.listdir(databaseDir):
                        if gene.lower() == file.split('_')[0].lower():
                            matchStdFile = os.path.join(databaseDir, file)
                    with open(matchStdFile, 'r') as std:
                        for record in SeqIO.parse(std, "fasta"):
                            stdSequenceLen[record.id] = len(record)


                    # matchStdFile = fnmatch.filter(os.listdir(databaseDir), f'{gene.upper()}*.fasta')
                    # with open(os.path.join(databaseDir, matchStdFile[0]), 'r') as std:
                    #     sequence_std = ''
                    #     for line_s in std:
                    #         if not line_s.startswith('>'):
                    #             sequence_std += line_s.strip()
                    #     length_std = int(len(sequence_std))



                    #### 1.先检查blast比对的结果能否覆盖起始位点（坐标1）和终止位点（标准序列长度的坐标），并且写入第一个exon
                    first_end_site = 0
                    if gene.lower() in exonLengthDict:
                        exonLength = exonLengthDict[gene.lower()][0]
                        with open(os.path.join(blastOutDir, outFIle), 'r') as f, open(multiple_exon_annotation_file, 'a') as ma:
                            start_marker = False
                            end_marker = False
                            for line in f:
                                queryName = line.split('\t')[0]
                                q_start = int(line.split('\t')[6])
                                q_end = int(line.split('\t')[7])
                                # if q_end == length_std:
                                #     end_marker = True
                                if q_start == 1:
                                    if (q_end - q_start) in range(exonLength - 15, exonLength + 15):
                                        start_marker = True
                                        readsName = line.split('\t')[1]
                                        length = int(line.split('\t')[3])
                                        q_end = int(line.split('\t')[7])
                                        first_end_site = q_end
                                        s_start = int(line.split('\t')[8])
                                        s_end = int(line.split('\t')[9])
                                        ma.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{s_start}\t{s_end}\n')
                                        print(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{s_start}\t{s_end}')
                                        break
                            if not start_marker:
                                f.seek(0)
                                for line in f:
                                    q_start = int(line.split('\t')[6])
                                    q_end = int(line.split('\t')[7])
                                    if q_start in range(0, 20):
                                        if exonLength - 20 <= (q_end - q_start) <= exonLength + 20:
                                            readsName = line.split('\t')[1]
                                            length = int(line.split('\t')[3])
                                            q_end = int(line.split('\t')[7])
                                            first_end_site = q_end
                                            s_start = int(line.split('\t')[8])
                                            s_end = int(line.split('\t')[9])
                                            ma.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{s_start}\t{s_end}\texon1 blast site lose{q_start - 1}\n')
                                            print(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{s_start}\t{s_end}\texon1 blast site lose{q_start - 1}')
                            if not first_end_site:
                                ma.write(
                                    f'{gene}\t----------\t----------\t----------\t----------\t----------\tfirst exon has problem and may annotation wrong\n')

                        #### 2.写入其他exon
                        if first_end_site:
                            with open(os.path.join(blastOutDir, outFIle), 'r') as f:
                                blast_txt = ''
                                for line in f:
                                    blast_txt += line
                                # with open(multiple_exon_annotation_file, 'a') as ma:
                                find_exon(blast_txt.split('\n'), first_end_site, gene.lower(), multiple_exon_annotation_file, n=0)









if __name__== '__main__':
    singleExonGene = {'atp1', 'atp4', 'atp6', 'atp8', 'atp9', 'ccmB', 'ccmC', 'ccmFn', 'cob', 'cox1', 'cox3', 'matR', 'mttB', 'nad3', 'nad4L', 'nad6', 'nad9', 'rpl5', 'rpl10', 'rpl16', 'rps7', 'rps12', 'rps13', 'sdh4'}
    multiExonGene = {'ccmFc', 'cox2', 'nad1', 'nad2', 'nad4', 'nad5', 'nad7'}
    exonLengthDict = {
        'ccmfc': [767, 547],
        'cox2': [382, 318, 83],
        'nad1': [385, 83, 192, 59, 259],
        'nad2': [153, 392, 161, 573, 188],
        'nad4': [461, 515, 423, 89],
        'nad5': [230, 1216, 22, 395, 150],
        'nad7': [143, 69, 467, 244, 262]
    }

    assembleDir, databaseDir, outDir = main()

    if not os.path.exists(assembleDir):
        sys.exit("Error: input!!! assemble directory does not exist")

    os.makedirs(outDir, exist_ok=True)

    run_blastn_for_assemble_results(assembleDir,databaseDir)
    annotation_for_single_exon_gene(assembleDir, outDir, singleExonGene)
    annotation_for_multi_exon_gene(assembleDir, outDir, multiExonGene, exonLengthDict)

# # 前两个参数最好使用绝对路径
# assembleDir = sys.argv[1] #assembleDir为批量处理的变量，指的是所有样本组装结果所在的目录（assembleDir是单个样本的目录）
# databaseDir = sys.argv[2] #输入标准序列所在的目录，后续考虑可以更改为数据库序列所在的目录
# outDir = sys.argv[3] #想要将结果输出的目录
#
# ### 测试运行的条件
# # 目录:/datb/workdir_zlzhang/yida_project/project/transcriptomeLnc_hjt/mitoAnalysis/geneAssemble/annotation/blast
# # 命令行:python blastn_assemble_gene_with_std_sequences.py /datb/workdir_zlzhang/yida_project/project/transcriptomeLnc_hjt/mitoAnalysis/geneAssemble/HJT120/Blast/FNQWHJTG01 /datb/workdir_zlzhang/yida_project/project/transcriptomeLnc_hjt/mitoAnalysis/pipline/database ./
#
#
#### 注释的逻辑：（1)先看后一个的开头是否等于、或者小于上一个的结尾；（2）比对的长度是否为参考外显子长度的±16bp左右
# def find_exon(txt, first_end_site, gene, n=0):
#     exonLengthDict = {
#         'ccmfc':[767, 547],
#         'cox2':[382, 318, 83],
#         'nad1':[385, 83, 192, 59, 259],
#         'nad2':[153, 392, 161, 573, 188],
#         'nad4':[461, 515, 423, 89],
#         'nad5':[230, 1216, 22, 395, 150],
#         'nad7':[143, 69, 467, 244, 262]
#     }
#     n += 1
#
#     if n + 1 <= len(exonLengthDict[gene]):
#         exon_length = int(exonLengthDict[gene][n])
#         # print(n, exon_length)
#         exon_marker = False
#
#         if exon_length == 22 and gene == 'nad5':
#             mismatch = 3
#             for line in txt:
#                 if not line.strip():
#                     continue
#                 start = int(line.split('\t')[6])
#                 end = int(line.split('\t')[7])
#                 if int(end-start+1) in range(20, 24):
#                     if start in range(first_end_site - 2, first_end_site + 2) and int(line.strip().split('\t')[4]) < mismatch:
#                         site_start = start
#                         site_end = end
#                         readsName = line.split('\t')[1]
#                         if int(line.split('\t')[8]) < int(line.split('\t')[9]):
#                             exon_start = int(line.split('\t')[8]) + int(first_end_site - start) + 1
#                             exon_end = int(line.split('\t')[9])
#                             length = int(abs(exon_end - exon_start)) + 1
#                         else:
#                             exon_start = int(line.split('\t')[8]) - int(first_end_site - start) - 1
#                             exon_end = int(line.split('\t')[9])
#                             length = int(abs(exon_end - exon_start)) + 1
#                         af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\n')
#                         print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
#                         exon_marker = True
#                         break
#             if not exon_marker:
#                 for line in txt:
#                     if not line.strip():
#                         continue
#                     start = int(line.split('\t')[6])
#                     end = int(line.split('\t')[7])
#                     if start > first_end_site:
#                         if exon_length - 2 <= end - start <= exon_length + 2:
#                             site_start = start
#                             site_end = end
#                             readsName = line.split('\t')[1]
#                             exon_start = int(line.split('\t')[8])
#                             exon_end = int(line.split('\t')[9])
#                             length = int(abs(exon_end - exon_start)) + 1
#                             af.write(
#                                 f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n} and {n + 1} has gap {start - first_end_site}\n')
#                             print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
#                             exon_marker = True
#             if not exon_marker:
#                 af.write(f'{gene} blastn out has no exon {n + 1} and may annotation wrong \n')
#                 return print(f'{gene} blastn out has no exon {n + 1} and may annotation wrong')
#             if exon_marker:
#                 return find_exon(txt, site_end, gene, n)
#         else:
#             for line in txt:
#                 if not line.strip():
#                     continue
#                 start = int(line.split('\t')[6])
#                 end = int(line.split('\t')[7])
#                 if start == first_end_site + 1:
#                     if exon_length - 20 <= end - start <= exon_length + 20:
#
#                         not_end_maker = False
#                         if n + 1 == len(exonLengthDict[gene]):
#                             std_length = 0
#                             for i in exonLengthDict[gene]:
#                                 std_length += int(i)
#                             if end != std_length:
#                                 not_end_maker = True
#                                 end_gap = std_length - end
#
#                         site_start = start
#                         site_end = end
#                         readsName = line.split('\t')[1]
#                         gap = int(line.split('\t')[5])
#                         exon_start = int(line.split('\t')[8])
#                         exon_end = int(line.split('\t')[9])
#                         length = int(abs(exon_end - exon_start)) + 1
#                         if not_end_maker:
#                             if gap == 0:
#                                 af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has gap {end_gap} with end\n')
#                                 print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has gap {end_gap} with end')
#                             else:
#                                 af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n + 1} has gap {end_gap} with end(2)exon {n + 1} has insert {gap} sites\n')
#                                 print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n + 1} has gap {end_gap} with end(2)exon {n + 1} has insert {gap} sites')
#                         else:
#                             if gap == 0:
#                                 af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\n')
#                                 print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
#                             else:
#                                 af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has insert {gap} sites\n')
#                                 print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has insert {gap} sites')
#
#                         exon_marker = True
#                         break
#             if not exon_marker:
#                 for line in txt:
#                     if not line.strip():
#                         continue
#                     gap = int(line.split('\t')[5])
#                     start = int(line.split('\t')[6])
#                     end = int(line.split('\t')[7])
#                     if start <= first_end_site:
#                         if exon_length - 20 <= end - start <= exon_length + 20 and (first_end_site - start <= 20):
#                             if gene == 'nad1' and gap:
#                                 align_length = int(line.split('\t')[3])
#                                 exon_start = int(line.split('\t')[8])
#                                 exon_end = int(line.split('\t')[9])
#                                 if exon_end - exon_start + 1 + gap == align_length:
#                                     site_start = first_end_site + 1
#                                     site_end = end
#                                     readsName = line.split('\t')[1]
#                                     if int(line.split('\t')[8]) < int(line.split('\t')[9]):
#                                         exon_start = int(line.split('\t')[8]) + int(first_end_site - start) + 1 - gap
#                                         exon_end = int(line.split('\t')[9])
#                                         length = int(abs(exon_end - exon_start)) + 1
#                                     else:
#                                         exon_start = int(line.split('\t')[8]) - int(first_end_site - start) - 1 + gap
#                                         exon_end = int(line.split('\t')[9])
#                                         length = int(abs(exon_end - exon_start)) + 1
#                                     af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\n')
#                                     print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
#                                     exon_marker = True
#                                     break
#                             else:
#
#                                 not_end_maker = False
#                                 if n + 1 == len(exonLengthDict[gene]):
#                                     std_length = 0
#                                     for i in exonLengthDict[gene]:
#                                         std_length += int(i)
#                                     if end != std_length:
#                                         not_end_maker = True
#                                         end_gap = std_length - end
#
#                                 site_start = first_end_site + 1
#                                 site_end = end
#                                 readsName = line.split('\t')[1]
#                                 if int(line.split('\t')[8]) < int(line.split('\t')[9]):
#                                     gap = int(line.split('\t')[5])
#                                     exon_start = int(line.split('\t')[8]) + int(first_end_site - start) + 1
#                                     exon_end = int(line.split('\t')[9])
#                                     length = int(abs(exon_end - exon_start)) + 1
#                                 else:
#                                     gap = int(line.split('\t')[5])
#                                     exon_start = int(line.split('\t')[8]) - int(first_end_site - start) - 1
#                                     exon_end = int(line.split('\t')[9])
#                                     length = int(abs(exon_end - exon_start)) + 1
#                                 if not_end_maker:
#                                     if gap == 0:
#                                         af.write( f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has gap {end_gap} with end\n')
#                                         print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has gap {end_gap} with end')
#                                     else:
#                                         af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n + 1} has gap {end_gap} with end(2)exon {n + 1} has insert {gap} sites\n')
#                                         print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n + 1} has gap {end_gap} with end(2)exon {n + 1} has insert {gap} sites')
#                                 else:
#                                     if gap == 0:
#                                         af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\n')
#                                         print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}')
#                                     else:
#                                         af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has insert {gap} sites\n')
#                                         print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n + 1} has insert {gap} sites')
#                                 exon_marker = True
#                                 break
#             if not exon_marker:
#                 for line in txt:
#                     if not line.strip():
#                         continue
#                     start = int(line.split('\t')[6])
#                     end = int(line.split('\t')[7])
#                     if start > first_end_site:
#                         if exon_length - 20 <= end - start <= exon_length + 20:
#
#                             not_end_maker = False
#                             if n + 1 == len(exonLengthDict[gene]):
#                                 std_length = 0
#                                 for i in exonLengthDict[gene]:
#                                     std_length += int(i)
#                                 if end != std_length:
#                                     not_end_maker = True
#                                     end_gap = std_length - end
#
#                             site_start = start
#                             site_end = end
#                             readsName = line.split('\t')[1]
#                             gap = int(line.split('\t')[5])
#                             exon_start = int(line.split('\t')[8])
#                             exon_end = int(line.split('\t')[9])
#                             length = int(abs(exon_end - exon_start)) + 1
#                             if not_end_maker:
#                                 if gap == 0:
#                                     af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n + 1} has gap {start - first_end_site}(2)exon {n + 1} has gap {end_gap} with end\n')
#                                     print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n + 1} has gap {start - first_end_site}(2)exon {n + 1} has gap {end_gap} with end')
#                                 else:
#                                     af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n + 1} has gap {start - first_end_site}(2)exon {n + 1} has gap {end_gap} with end(3)exon {n + 1} has insert {gap} sites\n')
#                                     print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n + 1} has gap {start - first_end_site}(2)exon {n + 1} has gap {end_gap} with end(3)exon {n + 1} has insert {gap} sites')
#                             else:
#                                 if gap == 0:
#                                     af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n} and {n+1} has gap {start - first_end_site}\n')
#                                     print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\texon {n} and {n+1} has gap {start - first_end_site}')
#                                 else:
#                                     af.write(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n+1} has gap {start - first_end_site}(2)exon {n + 1} has insert {gap} sites\n')
#                                     print(f'{gene}\t{readsName}\t----------\t{length}\t{exon_start}\t{exon_end}\t(1)exon {n} and {n+1} has gap {start - first_end_site}(2)exon {n + 1} has insert {gap} sites')
#                             exon_marker = True
#             if not exon_marker:
#                 af.write(f'{gene} blastn out has no exon {n+1}\n')
#                 return print(f'{gene} blastn out has no exon {n+1}')
#             if exon_marker:
#                 return find_exon(txt, site_end, gene, n)
#     else:
#         return None
#
#
#
# if __name__ == '__main__':
#
#     ### 批量处理assembleDir
#     for sample in os.listdir(assembleDir):
#         if os.path.isdir(os.path.join(assembleDir, sample)):
#             assembleDir = os.path.join(assembleDir, sample)
#
#             ### 多外显子基因元组
#             multiExonGene = ('ccmfc', 'cox2', 'nad1', 'nad2', 'nad4', 'nad5', 'nad7')
#             singleExonGene = {'atp1', 'atp4', 'atp6', 'atp8', 'atp9', 'ccmb', 'ccmc', 'ccmfn', 'cob', 'cox1', 'cox3', 'matr', 'mttb', 'nad3', 'nad4l', 'nad6', 'nad9', 'rpl10', 'rpl16', 'rpl5', 'rps12', 'rps13', 'rps7', 'sdh4'}
#
#
#             speName = assembleDir.split('/')[-1] # 注意输入的路径末尾不要有'/'
#             blastOutDir = os.path.join(outDir, sample, 'blastOut')
#             blastDatabaseDir = os.path.join(outDir, sample, 'database')
#             if not os.path.exists(blastOutDir):
#                 os.makedirs(blastOutDir, exist_ok=True)
#                 os.makedirs(blastDatabaseDir, exist_ok=True)
#
#             ### 批量blastn
#             for gene in os.listdir(assembleDir):
#                 ## 寻找spades组装的序列
#                 spadesDir = os.path.join(assembleDir, gene, 'spades')
#                 matchAssFile =  fnmatch.filter(os.listdir(spadesDir), '*.contigs.fasta')
#                 if not matchAssFile:
#                     continue
#                 assemble_sequence = os.path.join(spadesDir, matchAssFile[0])
#
#                 ## 对应寻找标准序列
#                 matchStdFile = fnmatch.filter(os.listdir(databaseDir), f'{gene.upper()}*.fasta')
#                 std_sequence_file = os.path.join(databaseDir, matchStdFile[0])
#                 if os.path.exists(std_sequence_file):
#                     if gene.lower() == 'nad5':
#                         # if sample in {'QHSLHJT01', 'QHTGHJT01', 'QHXYHJT01', 'WCQWHJTG02', 'WCQWHJTM02', 'YBCBHJT01', 'YBKYHJT02'}:
#                         #     blastn_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 6 -word_size 7 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble.out'
#                         #     subprocess.call(blastn_cmd, shell=True)
#
#                         if not os.path.exists(os.path.join(blastDatabaseDir, f'[{matchAssFile[0]}].ndb')):
#                             makeblastdb_cmd = f'makeblastdb -in {assemble_sequence} -dbtype nucl -out {blastDatabaseDir}/{matchAssFile}'
#                             subprocess.call(makeblastdb_cmd, shell=True)
#                             blastn_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 6 -word_size 7 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble.out'
#                             subprocess.call(blastn_cmd, shell=True)
#                             blastn0_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 0 -word_size 7 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble_0.out'
#                             subprocess.call(blastn0_cmd, shell=True)
#                         else:
#                             if not os.path.exists(os.path.join(blastOutDir, f'std_{gene}_with_assemble.out')):
#                                 blastn_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 6 -word_size 7 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble.out'
#                                 subprocess.call(blastn_cmd, shell=True)
#                             if not os.path.exists(os.path.join(blastOutDir, f'std_{gene}_with_assemble_0.out')):
#                                 blastn0_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 0 -word_size 7 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble_0.out'
#                                 subprocess.call(blastn0_cmd, shell=True)
#                     else:
#                         if not os.path.exists(os.path.join(blastDatabaseDir, f'[{matchAssFile[0]}].ndb')):
#                             makeblastdb_cmd = f'makeblastdb -in {assemble_sequence} -dbtype nucl -out {blastDatabaseDir}/{matchAssFile}'
#                             subprocess.call(makeblastdb_cmd, shell=True)
#                             blastn_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 6 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble.out'
#                             subprocess.call(blastn_cmd, shell=True)
#                             blastn0_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 0 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble_0.out'
#                             subprocess.call(blastn0_cmd, shell=True)
#                         else:
#                             if not os.path.exists(os.path.join(blastOutDir, f'std_{gene}_with_assemble.out')):
#                                 print(f'{assemble_sequence} has been in {blastDatabaseDir}')
#                                 blastn_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 6 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble.out'
#                                 subprocess.call(blastn_cmd, shell=True)
#                             if not os.path.exists(os.path.join(blastOutDir, f'std_{gene}_with_assemble_0.out')):
#                                 blastn0_cmd = f'blastn -query {std_sequence_file} -db {blastDatabaseDir}/{matchAssFile} -outfmt 0 -num_threads 32 -out {blastOutDir}/std_{gene}_with_assemble_0.out'
#                                 subprocess.call(blastn0_cmd, shell=True)
#
#             # ### 单外显子基因的注释(在blastout中，q_start永远小于q_end，但是s_start不一定小于s_end)
#             # #### 基本逻辑：1.找比对长度最长的条目；2.如果q_start起始不为1，且q_end的终止不为标准序列的长度的时候，提示可能组装序列存在问题
#             annotationDir = os.path.join(outDir, sample, 'annotation')
#             if not os.path.exists(annotationDir):
#                 os.makedirs(annotationDir, exist_ok=True)
#             singe_exon_annotation_file = os.path.join(annotationDir, 'single_exon_annotation.txt')
#
#             sa = open(singe_exon_annotation_file, 'w')
#             sa.write('gene\treadsName\tassemble_sequence_path\tlen_max\ts_start\ts_end\tinfo\n')
#
#
#             for gene in singleExonGene:
#                 geneMaker = False
#                 for file in os.listdir(blastOutDir):
#                     if gene in file.split('_')[1].lower():
#                         geneMaker = True
#                 if not geneMaker:
#                     sa.write(f'don not have {gene}\n')
#                     print(f'don not have {gene}')
#
#             for outFIle in fnmatch.filter(os.listdir(blastOutDir), '*assemble.out'):
#                 gene = outFIle.split('_')[1]
#                 if gene.lower() not in multiExonGene:
#
#                     #### 获取标准序列的长度
#                     matchStdFile = fnmatch.filter(os.listdir(databaseDir), f'{gene.upper()}*.fasta')
#                     with open(os.path.join(databaseDir, matchStdFile[0]), 'r') as std:
#                         sequence_std = ''
#                         for line_s in std:
#                             if not line_s.startswith('>'):
#                                 sequence_std += line_s.strip()
#                         length_std = int(len(sequence_std))
#
#
#                     blastInfo = {}
#                     blastInfoList = []
#                     len_max = 0
#                     mismatch = 100
#                     singleBlastInfo = ''
#                     with open(os.path.join(blastOutDir, outFIle), 'r') as f:
#                         for line in f:
#                             singleBlastInfo += line
#                             length = int(line.split('\t')[3])
#                             # if (length - len_max in range(-20, 20)) and (int(line.split('\t')[4]) < mismatch):
#                             #     len_max = length
#                             #     readsName = line.split('\t')[1]
#                             #     mismatch = int(line.split('\t')[4])
#                             #     q_start = int(line.split('\t')[6])
#                             #     q_end = int(line.split('\t')[7])
#                             #     s_start = int(line.split('\t')[8])
#                             #     s_end = int(line.split('\t')[9])
#                             if length > len_max:
#                                 if int(line.split('\t')[4]) < mismatch:
#                                     len_max = length
#                                     readsName = line.split('\t')[1]
#                                     mismatch = int(line.split('\t')[4])
#                                     insert = int(line.split('\t')[5])
#                                     q_start = int(line.split('\t')[6])
#                                     q_end = int(line.split('\t')[7])
#                                     s_start = int(line.split('\t')[8])
#                                     s_end = int(line.split('\t')[9])
#                             else:
#                                 if (length - length_std in range(-20,20)) and (int(line.split('\t')[4]) < mismatch):
#                                     len_max = length
#                                     readsName = line.split('\t')[1]
#                                     mismatch = int(line.split('\t')[4])
#                                     insert = int(line.split('\t')[5])
#                                     q_start = int(line.split('\t')[6])
#                                     q_end = int(line.split('\t')[7])
#                                     s_start = int(line.split('\t')[8])
#                                     s_end = int(line.split('\t')[9])
#                         if s_start and s_end and len_max:
#
#                             spadesDir = os.path.join(assembleDir, gene, 'spades')
#                             matchAssFile = fnmatch.filter(os.listdir(spadesDir), '*.contigs.fasta')
#                             assemble_sequence = os.path.join(spadesDir, matchAssFile[0])
#
#                             matchStdFile = fnmatch.filter(os.listdir(databaseDir), f'{gene.upper()}*.fasta')
#                             with open(os.path.join(databaseDir, matchStdFile[0]), 'r') as std:
#                                 sequence_std = ''
#                                 for line_s in std:
#                                     if not line_s.startswith('>'):
#                                         sequence_std += line_s.strip()
#                                 length_std = int(len(sequence_std))
#
#                             if length_std == len_max:
#                                 if insert and gene.lower() != 'matr':
#                                     sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tsequences has insert {insert} sites and may annotation wrong\n')
#                                     print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tsequences has insert {insert} sites and may annotation wrong')
#                                 else:
#                                     sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\n')
#                                     print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}')
#                             elif q_start != 1 and q_end != length_std:
#                                 q_maker = False
#                                 s_maker = False
#                                 for line in singleBlastInfo.split('\n'):
#                                     if line:
#                                         if int(line.strip().split('\t')[6]) == 1:
#                                             q_maker = True
#                                         if int(line.strip().split('\t')[7]) == length_std:
#                                             s_maker = True
#                                 if q_maker and not s_maker:
#                                     if insert:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {length_std - q_end}(2)start site lose {q_start - 1} sites and may annotation wrong(3)sequences has insert {insert} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {length_std - q_end}(2)start site lose {q_start - 1} sites and may annotation wrong(3)sequences has insert {insert} sites and may annotation wrong')
#                                     else:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {length_std - q_end}(2)start site lose {q_start - 1} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {length_std - q_end}(2)start site lose {q_start - 1} sites and may annotation wrong')
#                                 elif s_maker and not q_maker:
#                                     if insert:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites(2)end site lose {length_std - q_end} sites and may annotation wrong(3)sequences has insert {insert} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites(2)end site lose {length_std - q_end} sites and may annotation wrong(3)sequences has insert {insert} sites and may annotation wrong')
#                                 elif s_maker and q_maker:
#                                     if insert:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites and may annotation wrong(2)end site lose {length_std - q_end} sites and may annotation wrong(3)sequences has insert {insert} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites and may annotation wrong(2)end site lose {length_std - q_end} sites and may annotation wrong(3)sequences has insert {insert} sites and may annotation wrong')
#                                     else:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites(2)end site lose {length_std - q_end} sites\n')
#                                 if not q_maker and not s_maker:
#                                     if insert:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites(2)end site lose {length_std - q_end} sites(3)sequences has insert {insert} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites(2)end site lose {length_std - q_end} sites(3)sequences has insert {insert} sites and may annotation wrong')
#                                     else:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites(2)end site lose {length_std - q_end} sites\n')
#                                     print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites(2)end site lose {length_std - q_end} sites')
#                             elif q_start != 1:
#                                 maker = False
#                                 for line in singleBlastInfo.split('\n'):
#                                     if line:
#                                         if int(line.strip().split('\t')[6]) == 1:
#                                             maker = True
#                                 if not maker:
#                                     if insert:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites(2)sequences has insert {insert} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites(2)sequences has insert {insert} sites and may annotation wrong')
#                                     else:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tstart site lose {q_start - 1} sites\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tstart site lose {q_start - 1} sites')
#                                 else:
#                                     if insert:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites and may annotation wrong(2)sequences has insert {insert} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)start site lose {q_start - 1} sites and may annotation wrong(2)sequences has insert {insert} sites and may annotation wrong')
#                                     else:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tstart site lose {q_start - 1} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tstart site lose {q_start - 1} sites and may annotation wrong')
#                             elif q_end != length_std:
#                                 maker = False
#                                 for line in singleBlastInfo.split('\n'):
#                                     if line:
#                                         if int(line.strip().split('\t')[7]) == length_std:
#                                             maker = True
#                                 if not maker:
#                                     if insert:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {length_std - q_end} sites(2)sequences has insert {insert} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {length_std - q_end} sites(2)sequences has insert {insert} sites and may annotation wrong')
#                                     else:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tend site lose {length_std - q_end} sites\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tend site lose {length_std - q_end} sites')
#                                 else:
#                                     if insert:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {length_std - q_end} sites and may annotation wrong(2)sequences has insert {insert} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\t(1)end site lose {length_std - q_end} sites and may annotation wrong(2)sequences has insert {insert} sites and may annotation wrong')
#                                     else:
#                                         sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tend site lose {length_std - q_end} sites and may annotation wrong\n')
#                                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tend site lose {length_std - q_end} sites and may annotation wrong')
#                             else:
#                                 if len_max :
#                                     sa.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tsequence have {len_max - length_std} gap\n')
#                                     print(f'{gene}\t{readsName}\t{assemble_sequence}\t{len_max}\t{s_start}\t{s_end}\tsequence have {len_max - length_std} gap')
#
#
#
#
            # ### 多外显子基因的注释
            # #### 思考1：指定每个外显子的长度范围，如果不符合要求，则报错；(exon的顺序应该如何排列，如果遇到重复应该如何做)思考2：选取首尾可以相连的exon
            # exonLengthDict = {
            #     'ccmfc':[767, 547],
            #     'cox2':[382, 318, 83],
            #     'nad1':[385, 83, 192, 59, 259],
            #     'nad2':[153, 392, 161, 573, 188],
            #     'nad4':[461, 515, 423, 89],
            #     'nad5':[230, 1216, 22, 395, 150],
            #     'nad7':[143, 69, 467, 244, 262]
            # }
            #
            #
            #
            # annotationDir = os.path.join(outDir, sample, 'annotation')
            # if not os.path.exists(annotationDir):
            #     os.makedirs(annotationDir, exist_ok=True)
            # multiple_exon_annotation_file = os.path.join(annotationDir, 'multiple_exon_annotation.txt')
            #
            # with open(multiple_exon_annotation_file, 'w') as ma:
            #     ma.write('gene\treadsName\tassemble_sequence_path\texon_length\texon_start\texon_end\n')
            #
            # for outFIle in fnmatch.filter(os.listdir(blastOutDir), '*assemble.out'):
            #     gene = outFIle.split('_')[1]
            #     first_end_site = 0
            #
            # #### 获取标准序列的长度
            #     matchStdFile = fnmatch.filter(os.listdir(databaseDir), f'{gene.upper()}*.fasta')
            #     with open(os.path.join(databaseDir, matchStdFile[0]), 'r') as std:
            #         sequence_std = ''
            #         for line_s in std:
            #             if not line_s.startswith('>'):
            #                 sequence_std += line_s.strip()
            #         length_std = int(len(sequence_std))
            #
            # #### 获得组装序列的路径
            #     spadesDir = os.path.join(assembleDir, gene, 'spades')
            #     matchAssFile = fnmatch.filter(os.listdir(spadesDir), '*.contigs.fasta')
            #     assemble_sequence = os.path.join(spadesDir, matchAssFile[0])
            #
            # #### 1.先检查blast比对的结果能否覆盖起始位点（坐标1）和终止位点（标准序列长度的坐标），并且写入第一个exon
            #
            #     if gene.lower() in exonLengthDict:
            #         exonLength = exonLengthDict[gene.lower()][0]
            #         with open(os.path.join(blastOutDir, outFIle), 'r') as f, open(multiple_exon_annotation_file, 'a') as ma:
            #             start_marker = False
            #             end_marker = False
            #             for line in f:
            #                 q_start = int(line.split('\t')[6])
            #                 q_end = int(line.split('\t')[7])
            #                 if q_end == length_std:
            #                     end_marker = True
            #                 if q_start == 1:
            #                     if (q_end-q_start) in range(exonLength-15, exonLength+15):
            #                         if first_end_site:
            #                             continue
            #                         start_marker = True
            #                         readsName = line.split('\t')[1]
            #                         length = int(line.split('\t')[3])
            #                         q_end = int(line.split('\t')[7])
            #                         first_end_site = q_end
            #                         s_start = int(line.split('\t')[8])
            #                         s_end = int(line.split('\t')[9])
            #                         ma.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{s_start}\t{s_end}\n')
            #                         print(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{s_start}\t{s_end}')
            #             if not start_marker:
            #                 f.seek(0)
            #                 for line in f:
            #                     q_start = int(line.split('\t')[6])
            #                     q_end = int(line.split('\t')[7])
            #                     if q_start in range(0,20):
            #                         if exonLength-20 <= (q_end-q_start) <= exonLength+20:
            #                             readsName = line.split('\t')[1]
            #                             length = int(line.split('\t')[3])
            #                             q_end = int(line.split('\t')[7])
            #                             first_end_site = q_end
            #                             s_start = int(line.split('\t')[8])
            #                             s_end = int(line.split('\t')[9])
            #                             ma.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{s_start}\t{s_end}\texon1 blast site lose{q_start-1}\n')
            #                             print(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{s_start}\t{s_end}\texon1 blast site lose{q_start-1}')
            #             if not first_end_site:
            #                 ma.write(f'{gene}\t----------\t----------\t----------\t----------\t----------\tfirst exon has problem and may annotation wrong\n')
            #             # if not end_marker:
            #             #     ma.write(f'{gene}\tstart site has problem\n')
            #             #     print(f'{gene}\tblast has no end site')
            #
            # #### 2.写入其他exon
            #         if first_end_site:
            #             with open(os.path.join(blastOutDir, outFIle), 'r') as f:
            #                 blast_txt = ''
            #                 for line in f:
            #                     blast_txt += line
            #                 with open(multiple_exon_annotation_file, 'a') as af:
            #                     find_exon(blast_txt.split('\n'), first_end_site, gene.lower(), n=0)
            #
            #
            #
            #
            #                 # site_start = 1
            #                 # site_end = first_end_site
            #                 # n = 1
            #                 # while site_end != length_std:
            #                 #     f.seek(0)
            #                 #     for line in f:
            #                 #         q_start = int(line.split('\t')[6])
            #                 #         q_end = int(line.split('\t')[7])
            #                 #         if q_start == site_start + 1:
            #                 #             n += 1
            #                 #             readsName = line.split('\t')[1]
            #                 #             length = int(line.split('\t')[3])
            #                 #             exon_start = int(line.split('\t')[8])
            #                 #             exon_end = int(line.split('\t')[9])
            #                 #             site_start = q_start
            #                 #             site_end = q_end
            #                 #             ma.write(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{exon_start}\t{exon_end}\n')
            #                 #             print(f'{gene}\t{readsName}\t{assemble_sequence}\t{length}\t{exon_start}\t{exon_end}')
            #
            #
            #
            # # txt = '''
            # # 213\t679\t1\t1
            # # 1\t143\t1\t1
            # # 144\t212\t1\t1
            # # 923\t1185\t1\t1
            # # 679\t923\t1\t1
            # # 333\t475\t1\t1
            # # 601\t661\t1\t1
            # # '''
            # #
            # #
            # #
            # #
            # #
            # #
            # #
            # #
            # # exon_info = find_exon(txt.split('\n'), 143, 'nad7', multiple_exon_annotation_file)
            #
            #
            #
            #
