import os
import re
import random
import csv

# print("Hi")
# Aligned data: all_time0_linsi.fasta
def getquartiles(genome):
    A=0
    C=0
    G=0
    T=0
    for base in genome:
        if base=="a":
            A+=1
        if base=="c":
            C+=1
        if base=="g":
            G+=1
        if base=="t":
            T+=1
    Aportion=A/(A+C+G+T)
    Cportion=C/(A+C+G+T)
    Gportion=G/(A+C+G+T)
    Tportion=T/(A+C+G+T)
    portions=0
    # if Aportion<0.3929146538:
    #     portions+=0
    # elif Aportion<0.3976190476:
    #     portions+=1
    # elif Aportion<0.4030121575:
    #     portions+=2
    # else:
    #     portions+=3

    # if Cportion<0.1644776049:
    #     portions+=0
    # elif Cportion<0.1672131148:
    #     portions+=1
    # elif Cportion<0.1702432046:
    #     portions+=2
    # else:
    #     portions+=3

    # if Gportion<0.1989730001:
    #     portions+=0
    # elif Gportion<0.203030303:
    #     portions+=1
    # elif Gportion<0.2095808383:
    #     portions+=2
    # else:
    #     portions+=3

    if Tportion<0.2217453505:
        portions+=0
    elif Tportion<0.2324840764:
        portions+=1
    elif Tportion<0.2378640777:
        portions+=2
    else:
        portions+=3
    print(portions)
    return portions
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
def synonymous(codon1, codon2):
    cdn1=""
    cdn2=""
    for i in range(0,3):
        cdn1+=(randpicker(codon1[i]).upper())
        cdn2+=(randpicker(codon2[i]).upper())
    # if "N" in cdn1 or "N" in cdn2:
    #     return False
    codontoprotein = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
                      'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                      'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
                      'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
                      'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                      'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                      'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
                      'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                      'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
                      'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                      'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
                      'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                      'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                      'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                      'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
                      'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    if codontoprotein[cdn1]==codontoprotein[cdn2]:
        return True
    else:
        return False

def getsynonymous(codon1):
    print(codon1)
    cdn1=""
    for i in range(0,3):
        cdn1+=(randpicker(codon1[i]).upper())
    # if "N" in cdn1 or "N" in cdn2:
    #     return False
    codontoprotein = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
                      'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                      'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
                      'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
                      'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                      'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                      'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
                      'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                      'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
                      'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                      'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
                      'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                      'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                      'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                      'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
                      'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    return codontoprotein[cdn1]

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
    for i in range(0,256):
        mutations.append({})
    A=0
    C=0
    G=0
    T=0
    for i in range(0,256):
        mutations[i]["A>C"] = 0
        mutations[i]["T>G"] = 0
        mutations[i]["A>G"] = 0
        mutations[i]["T>C"] = 0
        mutations[i]["A>T"] = 0
        mutations[i]["T>A"] = 0
        mutations[i]["C>A"] = 0
        mutations[i]["G>T"] = 0
        mutations[i]["C>G"] = 0
        mutations[i]["G>C"] = 0
        mutations[i]["C>T"] = 0
        mutations[i]["G>A"] = 0
        mutations[i]["A"] = 0
        mutations[i]["C"] = 0
        mutations[i]["G"] = 0
        mutations[i]["T"] = 0
        mutations[i]["Tot"] = 0
        mutations[i]["q*"] = 0
        mutations[i]["qA"] = 0
        mutations[i]["qB"] = 0
        mutations[i]["qC"] = 0
        mutations[i]["qD"] = 0
        mutations[i]["qE"] = 0
        mutations[i]["qF"] = 0
        mutations[i]["qG"] = 0
        mutations[i]["qH"] = 0
        mutations[i]["qI"] = 0
        mutations[i]["qJ"] = 0
        mutations[i]["qK"] = 0
        mutations[i]["qL"] = 0
        mutations[i]["qM"] = 0
        mutations[i]["qN"] = 0
        mutations[i]["qO"] = 0
        mutations[i]["qP"] = 0
        mutations[i]["qQ"] = 0
        mutations[i]["qR"] = 0
        mutations[i]["qS"] = 0
        mutations[i]["qT"] = 0
        mutations[i]["qU"] = 0
        mutations[i]["qV"] = 0
        mutations[i]["qW"] = 0
        mutations[i]["qX"] = 0
        mutations[i]["qY"] = 0
        mutations[i]["qZ"] = 0
        mutations[i]["/*"] = 0
        mutations[i]["/A"] = 0
        mutations[i]["/B"] = 0
        mutations[i]["/C"] = 0
        mutations[i]["/D"] = 0
        mutations[i]["/E"] = 0
        mutations[i]["/F"] = 0
        mutations[i]["/G"] = 0
        mutations[i]["/H"] = 0
        mutations[i]["/I"] = 0
        mutations[i]["/J"] = 0
        mutations[i]["/K"] = 0
        mutations[i]["/L"] = 0
        mutations[i]["/M"] = 0
        mutations[i]["/N"] = 0
        mutations[i]["/O"] = 0
        mutations[i]["/P"] = 0
        mutations[i]["/Q"] = 0
        mutations[i]["/R"] = 0
        mutations[i]["/S"] = 0
        mutations[i]["/T"] = 0
        mutations[i]["/U"] = 0
        mutations[i]["/V"] = 0
        mutations[i]["/W"] = 0
        mutations[i]["/X"] = 0
        mutations[i]["/Y"] = 0
        mutations[i]["/Z"] = 0
        mutations[i]["/*"] = 0
        mutations[i]["AAA"] = 0
        mutations[i]["AAC"] = 0
        mutations[i]["AAG"] = 0
        mutations[i]["AAT"] = 0
        mutations[i]["ACA"] = 0
        mutations[i]["ACC"] = 0
        mutations[i]["ACG"] = 0
        mutations[i]["ACT"] = 0
        mutations[i]["AGA"] = 0
        mutations[i]["AGC"] = 0
        mutations[i]["AGG"] = 0
        mutations[i]["AGT"] = 0
        mutations[i]["ATA"] = 0
        mutations[i]["ATC"] = 0
        mutations[i]["ATG"] = 0
        mutations[i]["ATT"] = 0
        mutations[i]["CAA"] = 0
        mutations[i]["CAC"] = 0
        mutations[i]["CAG"] = 0
        mutations[i]["CAT"] = 0
        mutations[i]["CCA"] = 0
        mutations[i]["CCC"] = 0
        mutations[i]["CCG"] = 0
        mutations[i]["CCT"] = 0
        mutations[i]["CGA"] = 0
        mutations[i]["CGC"] = 0
        mutations[i]["CGG"] = 0
        mutations[i]["CGT"] = 0
        mutations[i]["CTA"] = 0
        mutations[i]["CTC"] = 0
        mutations[i]["CTG"] = 0
        mutations[i]["CTT"] = 0
        mutations[i]["GAA"] = 0
        mutations[i]["GAC"] = 0
        mutations[i]["GAG"] = 0
        mutations[i]["GAT"] = 0
        mutations[i]["GCA"] = 0
        mutations[i]["GCC"] = 0
        mutations[i]["GCG"] = 0
        mutations[i]["GCT"] = 0
        mutations[i]["GGA"] = 0
        mutations[i]["GGC"] = 0
        mutations[i]["GGG"] = 0
        mutations[i]["GGT"] = 0
        mutations[i]["GTA"] = 0
        mutations[i]["GTC"] = 0
        mutations[i]["GTG"] = 0
        mutations[i]["GTT"] = 0
        mutations[i]["TAA"] = 0
        mutations[i]["TAC"] = 0
        mutations[i]["TAG"] = 0
        mutations[i]["TAT"] = 0
        mutations[i]["TCA"] = 0
        mutations[i]["TCC"] = 0
        mutations[i]["TCG"] = 0
        mutations[i]["TCT"] = 0
        mutations[i]["TGA"] = 0
        mutations[i]["TGC"] = 0
        mutations[i]["TGG"] = 0
        mutations[i]["TGT"] = 0
        mutations[i]["TTA"] = 0
        mutations[i]["TTC"] = 0
        mutations[i]["TTG"] = 0
        mutations[i]["TTT"] = 0

    for filename in filenames:
        with open(path + filename, 'r') as test_file:
            full_file=test_file.read()
            genomes=re.split(r'>.*\n', full_file)
            #print("ekans", len(genomes))
            igloo=int(genomes[0].replace("\n",""))
            for index1 in range(1,len(genomes)-1):
                referenceSequence=""
                for line in genomes[index1]:
                    referenceSequence += line.split("\n")[0]
                for index2 in range(index1+1,index1+2):
                    testSequence = ""
                    for line in genomes[index2]:
                        testSequence += line.split("\n")[0]
                    if "-" in testSequence: ##if "-" in substr and not include:
                        continue
                    patientnum += 1
                    #print(filename)
                    if patientnum>=patientmax:
                        #print("f",count)
                        break
                    # print(testSequence[0:5],referenceSeqs[igloo:igloo+5])
                    if len(testSequence)!=len(referenceSequence):
                        continue
                    for i in range(0, len(testSequence)):
                        # if :
                        if testSequence[i]=="a":
                            A+=1
                        if testSequence[i]=="c":
                            C+=1
                        if testSequence[i]=="g":
                            G+=1
                        if testSequence[i]=="t":
                            T+=1
                        if "n" not in referenceSequence[i-(i+igloo)%3:i-(i+igloo)%3+3] and "-" not in referenceSequence[i-(i+igloo)%3:i-(i+igloo)%3+3] :
                            if (i+igloo)%3==0:
                                mutations[categorize(baseline[i][2],baseline[i][3])]["q"+getsynonymous(referenceSequence[i-(i+igloo)%3:i-(i+igloo)%3+3])]+=1
                        
                        if testSequence[i] != referenceSequence[i] and testSequence[i] != "-":
                            test = randpicker(testSequence[i])
                            ref = randpicker(referenceSequence[i])
                            testCodon= testSequence[i-(i+igloo)%3:i-(i+igloo)%3+3]
                            referenceCodon = referenceSequence[i-(i+igloo)%3:i-(i+igloo)%3+3]
                            if "n" in testCodon or "n" in referenceCodon or test == ref:#test == "n" or test == "N" or ref == "n" or ref == "N" or test == ref: #  or not mutations.keys().__contains__(ref.upper() + ">" + test.upper()):
                                continue
                            if len(testCodon)<3 or len(referenceCodon)<3:
                                continue
                            tCodonCopy=testCodon
                            
                            if synonymous(testCodon, referenceCodon):
                                continue
                            else:
                                burp=getsynonymous(referenceCodon)
                                slurp=getsynonymous(testCodon)
                            if i>1676:
                                break
                            # if float(baseline[i][2])>.1:
                            #     continue
                            composition=categorize(baseline[i][2],baseline[i][3])
                            context=composition#categorize(baseline[i+igloo][2],baseline[i+igloo][3])
                            mut = ref + ">" + test
                            mut = mut.upper()
                            ekans=""
                            for base in tCodonCopy:
                                ekans+=randpicker(base).upper()
                            print(ekans)

                            temp=mut
                            mut=context
                            context=temp
                            mutations[mut][context[0]] = mutations[mut][context[0]] + 1
                            
                            if context.__contains__("-"):
                                continue
                                #print(mut)
                            if False:#if not mutations[mut].__contains__(context):
                                # print("opposite search")
                                #print(opposite(mut[0]))

                                # mut = opposite(mut[0]) + ">" + opposite(mut[2])
                                # context = opposite(context[0]) + opposite(context[1]) + opposite(context[2])
                                # mutations[mut] = mutations[mut] + 1
                                #print(mut)
                                pass

                            else:
                                # print("mut", mut)
                                # print("context", context)
                                # print("arbok")
                                mutations[mut][context] = mutations[mut][context] + 1
                                mutations[mut]["Tot"] = mutations[mut]["Tot"] + 1
                                if test=="t":
                                    mutations[mut]["/"+ slurp] = mutations[mut]["/"+ slurp]+1
                                    print(i)
                                mutations[mut][ekans] = mutations[mut][ekans]+1
                                
# opening the csv file in 'a+' mode 
csv_file = "Names.csv"
try:
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=mutations[0].keys())
        writer.writeheader()
        for data in mutations:
            writer.writerow(data)
except IOError:
    print("I/O error")
            
            
print(A,C,G,T)



