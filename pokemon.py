import os
import random
import csv
import itertools
import shutil


# Used to generate frequency CSVs based on reference-based comparisons!!
# Used same reference for all studies -- may want to split by study for git file cap

def generate_frequencies(divisions, patientmax, absmin, absmax, drug):
    try:
        os.mkdir(drug+"."+str(divisions))
    except:
        pass
    # -- = insert study reference
    study = "2031"
    referenceGenome = open('../nucleotide2/reference_subregion_rt.txt', 'r')

    path = '../nucleotide2/'
    

    # convert genome to single string
    referenceSequence = ""
    with referenceGenome as reference_file:
        for line in reference_file.readlines():
            referenceSequence += line.split("\n")[0]
    referenceSequence = referenceSequence[absmin:absmax]
    
    for iterations in range(0,divisions):
        print(divisions)
        print(iterations)
        
        outputFile = open(drug+"."+str(divisions)+"/"+"MutationDivision"+str(iterations+1)+"."+str(divisions)+drug+".txt", "w")
        foundContexts = []

        firstPatient = True
        for (dirpath, dirnames, filenames) in os.walk(path):
            patientnum = 0
            pokemon = {}
            for filename in filenames:
                if str(filename).endswith(".txt") and (drug in str(filename)): #and str(filename).startswith(study):
                    patientnum += 1
                    with open(path + filename, 'r') as test_file:
                        if patientnum>=patientmax:
                            break
                        print(str(patientnum) + "/" + str(patientmax))
                            #if patientnum >= 40:
                            #return This is a limitation to the number of

                        
                        
                        
                        

                        testSequence = ""
                        for line in test_file.readlines():
                            testSequence += line.split("\n")[0]
                        testSequence=testSequence[absmin:absmax]
                        
                        startIndex =round( len(testSequence)*iterations/divisions)
                        endIndex =round (len(testSequence)*(iterations+1)/divisions)
                        #print(startIndex,endIndex)
                        #print(testSequence)
                        arbok=""
                        for i in testSequence:
                            arbok+=randpicker(i)
                        testSequence=arbok
                        #print(testSequence)
                        if "-" in testSequence[startIndex:endIndex]: ##if "-" in substr and not include
                            
                            print("failed")
                            continue
                        if startIndex<800:
                            string = translate(testSequence[startIndex-startIndex%3+4:endIndex-endIndex%3-2])
                        else:
                            string=string = translate(testSequence[startIndex-startIndex%3+2:endIndex-endIndex%3-1])
                        print(string)
                        try:
                            pokemon[string]+=1
                        except:
                            pokemon[string]=1
                        
        print("finishing")
        ash=sorted(pokemon, key=pokemon.__getitem__, reverse=True)
        pikachu={}
        total=0
        i=0
        for name in ash:
            i+=1
            total+=pokemon[name]
            if i>10:
                break
        i=0
        for name in ash:
            i+=1
            pikachu[name]=float(pokemon[name])/float(total)
            if i>10:
                break
        outputFile.write(str(pikachu))

def randpicker(letter):
    if letter == "u" or letter == "U":
        return "A"
    if letter == "t" or letter == "T":
        return "T"
    if letter == "a" or letter == "A":
        return "A"
    if letter == "g" or letter == "G":
        return "G"
    if letter == "c" or letter == "C":
        return "C"
    return "-"

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
                
def translate(seq):
       
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            #print(seq[i:i + 3])
            protein+=table[codon]
    return protein

generate_frequencies(3,72000,2,1680, "AZT")
# To do
