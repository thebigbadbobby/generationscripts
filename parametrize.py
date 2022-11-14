import dnds
import os
import csv 
def DnaCheck(sequence):
    return all(base in ('A', 'C', 'T', 'G') for base in sequence)

referenceGenome = open('../nucleotide2/reference_subregion_rt.txt', 'r')
referenceSequence = ""
with referenceGenome as reference_file:
    for line in reference_file.readlines():
        referenceSequence += line.split("\n")[0]
sequence_1 = "ACGACT"
sequence_2 = "ATCACT"
print(dnds.pnps(sequence_1,sequence_2))
path = '../nucleotide2/'
sequences=[]
omegas=[]
lenlimit=30
absmin=0
for seqnumber in range(0,559):
    sequences.append("")
for (dirpath, dirnames, filenames) in os.walk(path):
    patientnum = 0
    for filename in filenames:
        patientnum+=1
        if patientnum==72000:
            break
        print(patientnum)
        if str(filename).endswith(".txt"):
            with open(path + filename, 'r') as test_file:
                testSequence = ""
                for line in test_file.readlines():
                    testSequence += line.split("\n")[0]
                for baseidx in range(len(testSequence)-1,0,-1):
                    if referenceSequence[baseidx]=="-":
                        testSequence=testSequence[:baseidx]+testSequence[baseidx+1:]

                for codonindex in range(0,559):
                    codon=testSequence[absmin+3*codonindex:absmin+3*codonindex+3]
                    if DnaCheck(codon):
                        sequences[codonindex]+=codon
counter=-1
#print(referenceSequence)
for baseidx in range(len(referenceSequence)-1,0,-1):
    if referenceSequence[baseidx]=="-":
        referenceSequence=referenceSequence[:baseidx]+referenceSequence[baseidx+1:]
        #print(referenceSequence)


for sequence in sequences:
    counter+=1
    #print(counter)
    # if len(sequence_0)%2==1:
    #     sequence_0+=sequence_0[0:3]
    #print(int(len(sequence_0)/2))
    ACCURACY=2700
    sequence_1=sequence[0:min(len(sequence),3*ACCURACY)]
    refcodon=referenceSequence[absmin+3*counter:absmin+3*(counter+1)]
    # print(refcodon)
    # print(sequence_1[:3])
    sequence_2=""
    form=min(int(len(sequence)/3),ACCURACY)
    sequence_2=refcodon*form
    #if counter<3:
        #print("\n"+sequence_1, sequence_2)
    try:
        ratio=dnds.pnps(sequence_1,sequence_2)
        if ratio==0:
            omegas.append([int(ratio)])
            omegas.append([int(ratio)])
            omegas.append([int(ratio)])
            print(counter, sequence_1[:30], sequence_2[:30])
        else:
            omegas.append([float(ratio)])
            omegas.append([float(ratio)])
            omegas.append([float(ratio)])

    except:
        omegas.append(["Infinity"])
        omegas.append(["Infinity"])
        omegas.append(["Infinity"])
        #, sequence_2[30])
    if counter==558:
        print(counter, sequence_1[:30], sequence_2[:30])
print(referenceSequence)
print(omegas)

  
# data to be written row-wise in csv fil 
  
# opening the csv file in 'w+' mode 
file = open('dndsallreference.csv', 'w+', newline ='') 
  
# writing the data into the file 
with file:     
    write = csv.writer(file) 
    write.writerows(omegas) 

    

                
                    
