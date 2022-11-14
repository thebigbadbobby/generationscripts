import csv 
# opening the CSV file 
with open('mutationbaseline.csv', mode ='r') as file: 
    
    # reading the CSV file 
    csvFile = csv.reader(file) 
  
    # displaying the contents of the CSV file 
    for lines in csvFile: 
        print(lines)

    study = "2031"
    referenceGenome = open('../nucleotide2/reference_subregion_rt.txt', 'r')

    path = '../nucleotide2/'
