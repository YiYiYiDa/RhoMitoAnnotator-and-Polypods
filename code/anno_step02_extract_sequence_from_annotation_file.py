# Date: 20251128
# Author: Linchun Shi, Yanda Zhu
# Copyright: (c) 2025 Linchun Shi, Yanda Zhu. All rights reserved.

####################################################

import sys
import os
import argparse
from Bio import SeqIO
from pandas.core.groupby.numba_ import generate_numba_agg_func

####################################################

def arg_parse():
    parser = argparse.ArgumentParser(description='Extract sequence from annotation file')
    parser.add_argument('annotationResult', help='path to annotationResult')
    parser.add_argument('errorRecord', help='path to errorRecord')
    return parser.parse_args()

def main():
    args = arg_parse()
    return args.annotationResult, args.errorRecord

def extract_single_exon_gene_sequence_from_annotation_file(annotationResult, errorRecord):
    with open(errorRecord, 'w') as er:
        for sample in os.listdir(annotationResult):
            if os.path.isdir(os.path.join(annotationResult, sample, 'annotation')):
                annotationDir = os.path.join(annotationResult, sample, 'annotation')


            if 'single_exon_annotation.txt' in os.listdir(annotationDir):
                singleAnnotation = os.path.join(annotationDir, 'single_exon_annotation.txt')
                single_sequence = os.path.join(annotationDir, 'single_exon_gene.fasta')
            else:
                er.write(f'{sample}: singleAnnotation file not found\n')
                print(f'{sample}: singleAnnotation file not found')

            ### 单外显子序列提取
            with open(singleAnnotation, 'r') as s, open(single_sequence, 'w') as ss:
                for line in s:
                    if not 'gene	readsName	assemble_sequence_path	len_max	s_start	s_end' in line:
                        items = line.strip().split('\t')
                        if len(items) >= 6:
                            geneName = items[0]
                            ss.write(f'>{geneName}\n')
                            seqName = items[1]
                            seqPath = items[2]
                            s_site = int(items[4])
                            e_site = int(items[5])
                            for record in SeqIO.parse(seqPath, 'fasta'):
                                if seqName == record.id:
                                    if e_site > s_site:
                                        ss.write(str(record.seq[int(s_site) - 1:int(e_site)]) + '\n')
                                    else:
                                        ss.write(str(record.seq[int(e_site)-1:int(s_site)].reverse_complement()) + '\n')


                            if len(items) > 6:
                                er.write(f'{sample}: {geneName} {items[6]}\n')
                                print(f'{sample}: {geneName} {items[6]}')
                        else:
                            print(f'{sample}: {line.strip()}')
                            er.write(f'{sample}: {line.strip()}\n')
                            print(f'{sample}: {line.strip()}')

def extract_multiple_exon_gene_sequence_from_annotation_file(annotationResult, errorRecord):
    with open(errorRecord, 'a') as er:
        for sample in os.listdir(annotationResult):
            if os.path.isdir(os.path.join(annotationResult, sample, 'annotation')):
                annotationDir = os.path.join(annotationResult, sample, 'annotation')

                if 'multiple_exon_annotation.txt' in os.listdir(annotationDir):
                    multiAnnotation = os.path.join(annotationDir, 'multiple_exon_annotation.txt')
                    multi_sequence = os.path.join(annotationDir, 'multiple_exon_gene.fasta')
                else:
                    er.write(f'{sample}: multiAnnotation file not found\n')
                    print(f'{sample}: multiAnnotation file not found')

            ### 多外显子序列提取
            with open(multiAnnotation, 'r') as m, open(multi_sequence, 'w') as ms:
                geneName = ''
                siteList = []
                seqName = ''
                seqPath = ''
                contigsAss = ''
                geneSequence = ''
                for line in m:
                    if not 'gene	readsName	assemble_sequence_path	exon_length	exon_start	exon_end' in line:
                        items = line.strip().split('\t')
                        if len(items) >= 6:
                            if geneName != items[0].lower():
                                geneSequence = ''
                                contigsAss = ''
                                if siteList:
                                    for i in siteList:
                                        seqName = i[0]
                                        s_site = i[1]
                                        e_site = i[2]
                                        contigsAss = ''
                                        if not os.path.exists(seqPath):
                                            print(f'{sample}: {seqPath} file not found')
                                            print(line.split('\t')[2])
                                            er.write(f'{sample}: {seqPath} file not found\n')
                                            break
                                        with open(seqPath, 'r') as f:
                                            maker = False
                                            for line_s in f:
                                                if '>' in line_s:
                                                    if seqName in line_s:
                                                        maker = True
                                                else:
                                                    if maker:
                                                        contigsAss += line_s.strip()
                                            if contigsAss:
                                                if int(s_site) < int(e_site):
                                                    exon = contigsAss[int(s_site) - 1:int(e_site)]
                                                    geneSequence += exon
                                                else:
                                                    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                                                    contigsAss = ''.join(complement[n] for n in contigsAss)
                                                    if int(e_site) != 1:
                                                        exon = contigsAss[int(s_site) - 1:int(e_site) - 2:-1]
                                                        geneSequence += exon
                                                    else:
                                                        exon = contigsAss[int(s_site) - 1::-1]
                                                        geneSequence += exon
                                    ms.write(geneSequence + '\n')
                                    siteList = []

                                geneName = items[0].lower()
                                ms.write(f'>{geneName}\n')
                            if '---' not in items[2]:
                                seqPath = items[2]
                            s_site = items[4]
                            e_site = items[5]
                            seqName = items[1]
                            if seqPath:
                                siteList.append((seqName, s_site, e_site))
                        else:
                            print(f'{sample}: {line.strip()}')
                            er.write(f'{sample}: {line.strip()}\n')
                        if len(items) > 6:
                            print(f'{sample}: {geneName} {items[6]}')
                            er.write(f'{sample}: {geneName} {items[6]}\n')

                geneSequence = ''
                if siteList:
                    for i in siteList:
                        seqName = i[0]
                        s_site = i[1]
                        e_site = i[2]
                        contigsAss = ''
                        if not os.path.exists(seqPath):
                            print(f'{sample}: {seqPath} file not found')
                            print(line.split('\t')[2])
                            er.write(f'{sample}: {seqPath} file not found\n')
                            break
                        with open(seqPath, 'r') as f:
                            maker = False
                            for line_s in f:
                                if '>' in line_s:
                                    if seqName in line_s:
                                        maker = True
                                else:
                                    if maker:
                                        contigsAss += line_s.strip()
                            if contigsAss:
                                if int(s_site) < int(e_site):
                                    exon = contigsAss[int(s_site) - 1:int(e_site)]
                                    geneSequence += exon
                                else:
                                    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                                    contigsAss = ''.join(complement[i] for i in contigsAss)
                                    if int(e_site) != 1:
                                        exon = contigsAss[int(s_site) - 1:int(e_site) - 2:-1]
                                        geneSequence += exon
                                    else:
                                        exon = contigsAss[int(s_site) - 1::-1]
                                        geneSequence += exon
                    ms.write(geneSequence + '\n')




if __name__ == '__main__':
    annotationResult, errorRecord = main()
    extract_single_exon_gene_sequence_from_annotation_file(annotationResult, errorRecord)
    extract_multiple_exon_gene_sequence_from_annotation_file(annotationResult, errorRecord)


# annotationResult = sys.argv[1] #annotationResult表示样本结果所在目录（批量），annotationDir表示单个样本中注释文件所在的目录
#
# errorRecord = os.path.join(annotationResult, '..', 'errorRecord.txt')
#
# with open(errorRecord, 'w') as er:
#     for sample in os.listdir(annotationResult):
#         if os.path.isdir(os.path.join(annotationResult, sample, 'annotation')):
#             annotationDir = os.path.join(annotationResult, sample, 'annotation')
#
#             for file in os.listdir(annotationDir):
#                 if 'multiple_exon_annotation.txt' in file:
#                     multiAnnotation = os.path.join(annotationDir, file)
#                     multi_sequence = os.path.join(annotationDir, 'multiple_exon_gene.fasta')
#                 if 'single_exon_annotation.txt' in file:
#                     singleAnnotation = os.path.join(annotationDir, file)
#                     single_sequence = os.path.join(annotationDir, 'single_exon_gene.fasta')
#             if not os.path.exists(multiAnnotation):
#                 print(f'{sample}: multiAnnotation file not found')
#                 er.write(f'{sample}: multiAnnotation file not found\n')
#                 ptint(f'{sample}: multiAnnotation file not found')
#             if not os.path.exists(singleAnnotation):
#                 print(f'{sample}: singleAnnotation file not found')
#                 er.write(f'{sample}: singleAnnotation file not found\n')
#                 print(f'{sample}: singleAnnotation file not found')
#
#
#             ### 单外显子序列提取
#             with open(singleAnnotation, 'r') as s, open(single_sequence, 'w') as ss:
#                 for line in s:
#                     if not 'gene	readsName	assemble_sequence_path	len_max	s_start	s_end' in line:
#                         items = line.strip().split('\t')
#                         if len(items) >= 6:
#                             geneName = items[0]
#                             ss.write(f'>{geneName}\n')
#                             seqName = items[1]
#                             seqPath = items[2]
#                             s_site = items[4]
#                             e_site = items[5]
#                             with open(seqPath, 'r') as f:
#                                 maker = False
#                                 contigsAss = ''
#                                 for line in f:
#                                     if '>' in line:
#                                         maker = False
#                                         if seqName in  line:
#                                             maker = True
#                                     else:
#                                         if maker:
#                                             contigsAss += line.strip()
#
#                             if contigsAss:
#                                 if int(s_site) < int(e_site):
#                                     ss.write(contigsAss[int(s_site)-1:int(e_site)] + '\n')
#                                 else:
#                                     complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
#                                     contigsAss = ''.join(complement[i] for i in contigsAss)
#                                     if int(e_site) != 1:
#                                         ss.write(contigsAss[int(s_site)-1:int(e_site)-2:-1] + '\n')
#                                     else:
#                                         ss.write(contigsAss[int(s_site)-1::-1] + '\n')
#                             if len(items) > 6:
#                                 er.write(f'{sample}: {geneName} {items[6]}\n')
#                                 print(f'{sample}: {geneName} {items[6]}')
#                         else:
#                             print(f'{sample}: {line.strip()}')
#                             er.write(f'{sample}: {line.strip()}\n')
#                             print(f'{sample}: {line.strip()}')
#
#             ### 多外显子序列提取
#             with open(multiAnnotation, 'r') as m, open(multi_sequence, 'w') as ms:
#                 geneName = ''
#                 siteList = []
#                 seqName = ''
#                 seqPath = ''
#                 contigsAss = ''
#                 geneSequence = ''
#                 for line in m:
#                     if not 'gene	readsName	assemble_sequence_path	exon_length	exon_start	exon_end' in line:
#                         items = line.strip().split('\t')
#                         if len(items) >= 6:
#                             if geneName != items[0].lower():
#                                 geneSequence = ''
#                                 contigsAss = ''
#                                 if siteList:
#                                     for i in siteList:
#                                         seqName = i[0]
#                                         s_site = i[1]
#                                         e_site = i[2]
#                                         contigsAss = ''
#                                         if not os.path.exists(seqPath):
#                                             print(f'{sample}: {seqPath} file not found')
#                                             print(line.split('\t')[2])
#                                             er.write(f'{sample}: {seqPath} file not found\n')
#                                             break
#                                         with open(seqPath, 'r') as f:
#                                             maker = False
#                                             for line_s in f:
#                                                 if '>' in line_s:
#                                                     if seqName in line_s:
#                                                         maker = True
#                                                 else:
#                                                     if maker:
#                                                         contigsAss += line_s.strip()
#                                             if contigsAss:
#                                                 if int(s_site) < int(e_site):
#                                                     exon = contigsAss[int(s_site)-1:int(e_site)]
#                                                     geneSequence += exon
#                                                 else:
#                                                     complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
#                                                     contigsAss = ''.join(complement[n] for n in contigsAss)
#                                                     if int(e_site) != 1:
#                                                         exon = contigsAss[int(s_site)-1:int(e_site)-2:-1]
#                                                         geneSequence += exon
#                                                     else:
#                                                         exon = contigsAss[int(s_site)-1::-1]
#                                                         geneSequence += exon
#                                     ms.write(geneSequence + '\n')
#                                     siteList = []
#
#
#                                 geneName = items[0].lower()
#                                 ms.write(f'>{geneName}\n')
#                             if '---' not in items[2]:
#                                 seqPath = items[2]
#                             s_site = items[4]
#                             e_site = items[5]
#                             seqName = items[1]
#                             if seqPath:
#                                 siteList.append((seqName, s_site, e_site))
#                         else:
#                             print(f'{sample}: {line.strip()}')
#                             er.write(f'{sample}: {line.strip()}\n')
#                         if len(items) > 6:
#                             print(f'{sample}: {geneName} {items[6]}')
#                             er.write(f'{sample}: {geneName} {items[6]}\n')
#
#
#                 geneSequence = ''
#                 if siteList:
#                     for i in siteList:
#                         seqName = i[0]
#                         s_site = i[1]
#                         e_site = i[2]
#                         contigsAss = ''
#                         if not os.path.exists(seqPath):
#                             print(f'{sample}: {seqPath} file not found')
#                             print(line.split('\t')[2])
#                             er.write(f'{sample}: {seqPath} file not found\n')
#                             break
#                         with open(seqPath, 'r') as f:
#                             maker = False
#                             for line_s in f:
#                                 if '>' in line_s:
#                                     if seqName in line_s:
#                                         maker = True
#                                 else:
#                                     if maker:
#                                         contigsAss += line_s.strip()
#                             if contigsAss:
#                                 if int(s_site) < int(e_site):
#                                     exon = contigsAss[int(s_site) - 1:int(e_site)]
#                                     geneSequence += exon
#                                 else:
#                                     complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
#                                     contigsAss = ''.join(complement[i] for i in contigsAss)
#                                     if int(e_site) != 1:
#                                         exon = contigsAss[int(s_site) - 1:int(e_site) - 2:-1]
#                                         geneSequence += exon
#                                     else:
#                                         exon = contigsAss[int(s_site) - 1::-1]
#                                         geneSequence += exon
#                     ms.write(geneSequence + '\n')
# #
# #
# #
