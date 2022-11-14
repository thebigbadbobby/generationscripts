import os
import re
import random
import csv

# print("Hi")
# Aligned data: all_time0_linsi.fasta

def randpicker(letter):
    randnum = random.randint(0, 12)
    if letter == "u" or letter == "U":
        return "t"
    if letter == "t" or letter == "T":
        return "t"
    if letter == "a" or letter == "A":
        return "a"
    if letter == "g" or letter == "G":
        return "g"
    if letter == "c" or letter == "C":
        return "c"
    if letter == "n" or letter == "N":
        return "n"
    if letter == "r" or letter == "R":
        if randnum > 6:
            return "a"
        else:
            return "g"
    if letter == "y" or letter == "Y":
        if randnum > 6:
            return "c"
        else:
            return "t"
    if letter == "s" or letter == "S":
        if randnum > 6:
            return "g"
        else:
            return "c"
    if letter == "w" or letter == "W":
        if randnum > 6:
            return "a"
        else:
            return "t"
    if letter == "k" or letter == "K":
        if randnum > 6:
            return "g"
        else:
            return "t"
    if letter == "m" or letter == "M":
        if randnum > 6:
            return "a"
        else:
            return "c"
    if letter == "b" or letter == "B":
        if randnum > 4:
            return "c"
        if randnum > 8:
            return "g"
        else:
            return "t"
    if letter == "d" or letter == "D":
        if randnum > 4:
            return "a"
        if randnum > 8:
            return "g"
        else:
            return "t"
    if letter == "h" or letter == "H":
        if randnum > 4:
            return "a"
        if randnum > 8:
            return "c"
        else:
            return "t"
    if letter == "v" or letter == "V":
        if randnum > 4:
            return "a"
        if randnum > 8:
            return "c"
        else:
            return "g"

    else:
        return "c"

def populateContexts(mutations):
    left = "A"
    right = "A"
    for key in mutations.keys():
        for i in range(4):
            if i % 4 == 0:
                left = "A"
            if i % 4 == 1:
                left = "C"
            if i % 4 == 2:
                left = "G"
            if i % 4 == 3:
                left = "T"
            for j in range(4):
                if j % 4 == 0:
                    right = "A"
                if j % 4 == 1:
                    right = "C"
                if j % 4 == 2:
                    right = "G"
                if j % 4 == 3:
                    right = "T"
                context = left + key[0] + right
                mutations[key][context] = 0

def categorize(omega, pairprob):
    # pairbiggest=0.95
    # pairbig=0.65
    # pairmid=0.28
    # pairsmall=0.08
    # omegamid=0.1
    tier=int(float(pairprob)/.1)
    return tier

# Split by time

baseline=[]
with open('mutationbaseline.csv', mode ='r') as file: 
    
# reading the CSV file 
    baseline = list(csv.reader(file))
# were already in directory: HIVSigs.py		all_time0_linsi.fasta

patientmax=3000
patientnum=0
outputfile = open('out.txt', 'w')
# for file in os.listdir("/Users/macbook/Desktop/Proj6/HIVMutationSignatures/hiv_longitudinal/nucleotide"):
#     with open('/Users/macbook/Desktop/Proj6/HIVMutationSignatures/hiv_longitudinal/nucleotide/' + file, 'r') as content_file:
#         content = content_file.read()
#         # outputfile.write(content + "\n\n\n")
referenceSeqs = "ccaattagtcctattgagactgtaccagtaaaattaaagccaggaatggatggcccaaaggttaaacaatggccattgacagaagaaaaaataaaagcattaacagaaatttgtacagagatggaaaaggaaggaaaaatttcaaaaattgggcctgaaaatccatacaatactccaatatttgcgataaagaaaaaagatagtactaaatggaggaaattagtagatttcagagagctcaataaaagaacacaagacttttgggaagttcaattaggaataccgcatccagcgggcctaaaaaagaaaaaatcagtaacagtactagatgtgggggacgcatatttttcagttcctttacatgaaagctttagaaagtataccgcattcaccatacctagtacaaacaatgagacaccaggaatcaggtatcagtacaatgtgcttccacagggatggaaaggatcaccggcaatattccagagtagcatgacaaaaatcttagagccctttagatcaaaaaatccagaaataattatctatcaatacatggatgacttgtatgtaggatctgatttagaaatagggcagcatagaacaaaaatagaagagttaagagctcatctattgagctggggatttactacaccagacaaaaagcatcagaaagaacctccattcctttggatgggatatgagctccatcctgacaagtggacgtccagcctataatgctgccagaaaaagaaagctggactgtcaatgatatacagaaattagtggggaaactaaattgggcaagtcaaatttatgcagggattaaagtaaagcaattgtgtaaactcctcaggggagccaaagcactaacagatatagtaacattgactgaggaagcagaattagaattggcagagaacagggagattctaaaagaccctgtgcatggggtatattatgacccatcaaaggacttaatagcagaaatacagaaacaagggcaagaccaatggacatatcaaatttatcaagagccatttaaaaatctaaaaacagggaagtatgcaagaaaaaggtctgctcacactaatgatgtaaaacaattagcagaagtggtgcaaaaggtggtcatggaaagcatagtaatatggggaaagactcctaaatttaaactacccatacaaaaagaaacatgggaaacatggtgggtggactattggcaggctacctggattcctgaatgggagtttgtcaatacccctcctctagtaaaattgtggtaccaattagagaaagaccccatagcaggagcagagactttctatgtagatggggcagccaatagggagactaagctaggaaaagcagggtatgtcactgacaggggaagacaaaaggttgtttccctaactgagacaacaaatcaaaagactgaactacatgcaatccatctagccttgcaggattcaggatcagaagtaaatatagtaacagactcacagtatgcattaggaatcattcaggcacaaccagacagaagtgaatcagagttagtcaatcaaataatagagaagctaataggaaaggacaaagtctacctgtcatgggtaccagcacacaagggaattggaggaaatgaacaagtagataaattagttagctctggaatcaggaagatacta"
path="../nucleotide/"
for (dirpath, dirnames, filenames) in os.walk(path):
    patientnum = 0
    mutations = []
    for filename in filenames:
        with open(path + filename, 'r') as test_file:
            full_file=test_file.read()
            genomes=re.split(r'>.*\n', full_file)
            #print("ekans", len(genomes))
            igloo=int(genomes[0].replace("\n",""))
            if igloo>10:
                continue
            for index1 in range(1,len(genomes)):
                referenceSequence=""
                A=0
                C=0
                G=0
                T=0
                for line in genomes[index1]:
                    referenceSequence += line.split("\n")[0]
                print(len(referenceSequence))
                if len(referenceSequence)<600:
                        continue
                referenceSequence=referenceSequence[:600]
                for base in referenceSequence:
                    print("ekans")
                    if base=="a":
                        A+=1
                    if base=="c":
                        C+=1
                    if base=="g":
                        G+=1
                    if base=="t":
                        T+=1
                mutations.append({"A":A,"C":C,"G":G,"T":T})
# opening the csv file in 'a+' mode 
csv_file = "Names.csv"
try:
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["A","C","G","T"])
        writer.writeheader()
        for data in mutations:
            writer.writerow(data)
except IOError:
    print("I/O error")
            
            
print(A,C,G,T)

