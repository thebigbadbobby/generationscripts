import os
import re
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

def line_prepender(filename, line):
    with open(path + filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)    

    # convert genome to single string
referenceSequence = "ccaattagtcctattgagactgtaccagtaaaattaaagccaggaatggatggcccaaaggttaaacaatggccattgacagaagaaaaaataaaagcattaacagaaatttgtacagagatggaaaaggaaggaaaaatttcaaaaattgggcctgaaaatccatacaatactccaatatttgcgataaagaaaaaagatagtactaaatggaggaaattagtagatttcagagagctcaataaaagaacacaagacttttgggaagttcaattaggaataccgcatccagcgggcctaaaaaagaaaaaatcagtaacagtactagatgtgggggacgcatatttttcagttcctttacatgaaagctttagaaagtataccgcattcaccatacctagtacaaacaatgagacaccaggaatcaggtatcagtacaatgtgcttccacagggatggaaaggatcaccggcaatattccagagtagcatgacaaaaatcttagagccctttagatcaaaaaatccagaaataattatctatcaatacatggatgacttgtatgtaggatctgatttagaaatagggcagcatagaacaaaaatagaagagttaagagctcatctattgagctggggatttactacaccagacaaaaagcatcagaaagaacctccattcctttggatgggatatgagctccatcctgacaagtggacgtccagcctataatgctgccagaaaaagaaagctggactgtcaatgatatacagaaattagtggggaaactaaattgggcaagtcaaatttatgcagggattaaagtaaagcaattgtgtaaactcctcaggggagccaaagcactaacagatatagtaacattgactgaggaagcagaattagaattggcagagaacagggagattctaaaagaccctgtgcatggggtatattatgacccatcaaaggacttaatagcagaaatacagaaacaagggcaagaccaatggacatatcaaatttatcaagagccatttaaaaatctaaaaacagggaagtatgcaagaaaaaggtctgctcacactaatgatgtaaaacaattagcagaagtggtgcaaaaggtggtcatggaaagcatagtaatatggggaaagactcctaaatttaaactacccatacaaaaagaaacatgggaaacatggtgggtggactattggcaggctacctggattcctgaatgggagtttgtcaatacccctcctctagtaaaattgtggtaccaattagagaaagaccccatagcaggagcagagactttctatgtagatggggcagccaatagggagactaagctaggaaaagcagggtatgtcactgacaggggaagacaaaaggttgtttccctaactgagacaacaaatcaaaagactgaactacatgcaatccatctagccttgcaggattcaggatcagaagtaaatatagtaacagactcacagtatgcattaggaatcattcaggcacaaccagacagaagtgaatcagagttagtcaatcaaataatagagaagctaataggaaaggacaaagtctacctgtcatgggtaccagcacacaagggaattggaggaaatgaacaagtagataaattagttagctctggaatcaggaagatacta"
path="../nucleotide/"
counter=1
for (dirpath, dirnames, filenames) in os.walk(path):
    for filename in filenames:
        counter+=1
        with open(path + filename, 'r') as test_file:
            full_file=test_file.read()
            genomes=re.split(r'>.*\n', full_file)
            alignments = pairwise2.align.globalms(referenceSequence, genomes[1].replace("\n", ""), 2, -1, -30, -.1)
            res = alignments[0][1].find(next(filter(str.isalpha, alignments[0][1])))-genomes[1].replace("\n", "").find(next(filter(str.isalpha, genomes[1].replace("\n", ""))))
            line_prepender(filename, str(res))
            print(alignments[0][1])

