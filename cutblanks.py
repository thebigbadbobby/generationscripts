import os
import re

def line_prepender(filename, line):
    with open(path + filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)    

    # convert genome to single string
path="../nucleotide2/"
with open(path + "reference_subregion_rt.txt", 'r') as reference_file:
    full_file=reference_file.read()
    indexes = [i for i in range(len(full_file)) if full_file.startswith("-", i)] 
counter=1
for (dirpath, dirnames, filenames) in os.walk(path):
    for filename in filenames:
        counter+=1
        if filename[0]==".":
            continue
        with open(path + filename, 'r+') as test_file:
            print(filename)
            full_file=test_file.read()
            for index in indexes[::-1]:
                full_file=full_file[:index]+full_file[index+1:]
            test_file.seek(0,0)    
            test_file.write(full_file)
            test_file.truncate()