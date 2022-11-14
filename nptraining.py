import os
import re
import random
import csv
import json
# print("Hi")
# Aligned data: all_time0_linsi.fasta



# Split by time

# baseline=[]
# with open('mutationbaseline.csv', mode ='r') as file: 
    
# # reading the CSV file 
#     baseline = list(csv.reader(file))
# # were already in directory: HIVSigs.py		all_time0_linsi.fasta
def saveµtoµfile(filename, text):
    outFile = open(filename,'w+')
    outFile.write(text)

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
    mutations = {"x_training": [], "y_training": []}
    for filename in filenames:
        with open(path + filename, 'r') as test_file:
            full_file=test_file.read()
            identifiers=re.split(r'\n[a-z|\n|-]*', full_file)
            genomes=re.split(r'>.*\n', full_file)
            #print("ekans", len(genomes))
            igloo=int(genomes[0].replace("\n",""))
            if igloo>10:
                continue
            for index1 in range(1,len(genomes)):
                # print(genomes)
                referenceSequence=""
                A=0
                C=0
                G=0
                T=0
                for line in genomes[index1]:
                    referenceSequence += line.split("\n")[0]
                if "-" in referenceSequence: ##if "-" in substr and not include:
                        continue
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
                mutations["x_training"].append(referenceSequence[10-igloo:10-igloo+540])
                mutations["y_training"].append(identifiers[index1][1:])

textfile="training.txt"
saveµtoµfile(textfile, json.dumps(mutations))

# opening the csv file in 'w' write mode 
csv_file = "Names.csv"
try:
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["x_training","y_training"])
        writer.writeheader()
        for i in range(0,len(mutations["x_training"])):
            writer.writerow({"x_training": mutations["x_training"][i], "y_training": mutations["y_training"][i]})
except IOError:
    print("I/O error")
            
            
print(A,C,G,T)
