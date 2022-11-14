import os
import random
import csv
import itertools
from tempfile import NamedTemporaryFile
import shutil

LAPLACE_K = 0.3
# Used to generate frequency CSVs based on reference-based comparisons!!
# Used same reference for all studies -- may want to split by study for git file cap

def generate_frequencies(divisions, mod, patientmax, absmin, absmax, include, laPlace, hist):
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
        outputFile = open("MutationDivision"+str(iterations+1)+"."+str(divisions)+".csv", "w")
        count=0
        firstPatient = True
        while count<(absmax-absmin)/divisions:
            print(count)
            foundContexts = []
            
            for (dirpath, dirnames, filenames) in os.walk(path):
                patientnum = 0
                mutations = {}
                mutations["T>A"] = {}
                mutations["T>C"] = {}
                mutations["T>G"] = {}
                mutations["C>A"] = {}
                mutations["C>T"] = {}
                mutations["C>G"] = {}
                mutations["A>C"] = {}
                mutations["A>G"] = {}
                mutations["A>T"] = {}
                mutations["C>A"] = {}
                mutations["C>T"] = {}
                mutations["C>G"] = {}

                for filename in filenames:
                    if str(filename).endswith(".txt"): #and str(filename).startswith(study):
                        with open(path + filename, 'r') as test_file:
                            #print("opening file", patientnum)
                                #if patientnum >= 40:
                                #return This is a limitation to the number of

                            
                            
                            
                            
                            if patientnum%mod==0:
                                populateContexts(mutations)

                            testSequence = ""
                            for line in test_file.readlines():
                                testSequence += line.split("\n")[0]
                            testSequence=testSequence[absmin:absmax]
                            startIndex = max(1-absmin,len(testSequence)*iterations/divisions + count)
                            endIndex = min(len(testSequence)*iterations/divisions + count+hist, len(testSequence)*(iterations+1)/divisions)
                            if "-----" in testSequence[startIndex:endIndex] and not include: ##if "-" in substr and not include:
                                continue
                            patientnum += 1
                            if patientnum>=patientmax:
                                #print("f",count)
                                break
                            for i in range(startIndex, endIndex):
                                    # find mutations in line
                                    if testSequence[i] != referenceSequence[i] and testSequence[i] != "-":
                                        test = randpicker(testSequence[i])
                                        ref = randpicker(referenceSequence[i])
                                        #print(test)
                                        #print(test.upper())
                                        if test == "n" or ref == "n" or test == ref: #  or not mutations.keys().__contains__(ref.upper() + ">" + test.upper()):
                                            continue
                                        
                                        # get context of test
                                        # print("Mutation at nucleotide " + str(i) + " between: " + ref + " and: " + test)
                                        context = ""
                                        if i - 1 < 0:
                                            context = "-" + ref + referenceSequence[i+1]
                                            # print("Found a starting mutation!")
                                        if i + 1 >= len(referenceSequence):
                                            context = referenceSequence[i-1] + ref + "-"
                                            # print("Found an ending mutation!")
                                        else:
                                            k=1
                                            while referenceSequence[i-k] == "-" and k<4 and i>=k:
                                                k+=1
                                            q=1
                                            while referenceSequence[i+q] == "-" and q<4 and i+q<len(referenceSequence):
                                                q+=1
                                            context = referenceSequence[i-k] + ref + referenceSequence[i+q]
                                        #print(context)
                                        # print("Context " + context.upper() + "of mutation at i = " + str(i))
                                        mut = ref + ">" + test
                                        mut = mut.upper()
                                        
                                        context = context.upper()
                                        
                                        
                                        if context.__contains__("-"):
                                            continue
                                         #print(mut)
                                        if not mut in mutations:#if not mutations[mut].__contains__(context):
                                           # print("opposite search")
                                            #print(opposite(mut[0]))
                                            mut = opposite(mut[0]) + ">" + opposite(mut[2])
                                            context = opposite(context[0]) + opposite(context[1]) + opposite(context[2])
                                            mutations[mut][context] = mutations[mut][context] + 1

                                        else:
                                            mutations[mut][context] = mutations[mut][context] + 1

                        #print("Mutations for file" + filename + " : " + str(mutations))
                        
                        if (patientnum+1)%mod==0:
                            if laPlace:
                                mutations = laplace_smoothing(mutations, LAPLACE_K)
                                #print(patientnum)
                            if firstPatient is True:
                                #print("a")
                                with open("MutationDivision"+str(iterations+1)+"."+str(divisions)+".csv", "w") as f:
                                    w = csv.writer(f)
                                    w.writerow(["Mutation Type"] + ["Trinucleotide"] + ["Sample" + str(1)])
                                    for key in sorted(mutations.keys()):
                                        for context in mutations[key].keys():
                                            w.writerow([key] + [context] + [mutations[key][context]])
                                            foundContexts.append([key, context])
                                firstPatient = False
                                print("b")
                            else:  # move on to new column for next patient
                                
                                filename = "MutationDivision"+str(iterations+1)+"."+str(divisions)+".csv"
                                #print("Looking at patient number " + str(patientnum))
                                # now populate data
                                with open(filename) as csvFile:
                                    with open("dummyfile.csv", "w") as tempfile:
                                        # the following executes for every extracted mutation-context frequency
                                        # found from this patient
                                        csv_reader = csv.reader(csvFile, delimiter=',')
                                        csv_writer = csv.writer(tempfile)
                                        line_count = 0
                                        # Go over already found contexts
                                        
                                        for row in csv_reader:
                                            line_count += 1
                                            
                                            if line_count == 1:
                                                #csv_writer.writerow(row + ["Sample" + str((patientnum+1)/mod)])
                                                print("e",count)
                                                csv_writer.writerow(row + ["Sample" + str(patientnum/mod) + str(count)])
                                                #print("d", count)
                                                continue
                                            if mutations[row[0]].__contains__(row[1]):
                                                csv_writer.writerow(row + [mutations[row[0]][row[1]]])
                                                foundContexts.append([row[0], row[1]])
                                            else:
                                                csv_writer.writerow(row + [0])

                                        for key in sorted(mutations.keys()):
                                            for context in mutations[key].keys():
                                                if not foundContexts.__contains__([key, context]):
                                                    row = [key] + [context]
                                                    for i in range(patientnum - 1):
                                                        row += [0]
                                                    csv_writer.writerow(row + [mutations[key][context]])
                                shutil.move(tempfile.name, filename)
                            if hist != absmin-absmax and patientnum+1==patientmax:
                                #print("c")
                                count+=hist
                                #print(count)

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

def opposite(char):
    char = char.upper()
    if char == "A":
        return "T"
    if char == "T":
        return "A"
    if char == "C":
        return "G"
    if char == "G":
        return "C"
    if char == ">":
        return ">"

def laplace_smoothing(mutations, k):
    # find # total mutations
    num_total_mutations = 0
    for key in mutations.keys():
        for context in mutations[key].keys():
            num_total_mutations += mutations[key][context]
    #print("got here")
    # generate laplace values
    laplaceTable = {}
    for key in mutations.keys():
        laplaceTable[key] = {}
        #print("Assigning key : " + key)

    for key in mutations.keys():
        for context in mutations[key].keys():
            laplaceTable[key][context] = int(round(laplace(mutations[key][context], 96, num_total_mutations, k) * 1000))
    return laplaceTable


def laplace(count_x, magnitude_X, N, k):
    print(count_x,(count_x + k) / (N + k * magnitude_X))
    return (count_x + k) / (N + k * magnitude_X)

generate_frequencies(1, 10000, 10000, 0, 1710, True, False, 190)
# To do
