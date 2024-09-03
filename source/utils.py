# Functions and classes for ampAnalysis.py
import os
import subprocess
import time
from collections import defaultdict
import pandas as pd
from Bio import SeqIO

#----------------------------------------------------------------------------------------------------------------------------------------
# Create and fill progress file
#----------------------------------------------------------------------------------------------------------------------------------------

def createProgressFile(project_name):
    progressFile = open('./' + project_name + '/progress_file.txt', 'x')
    progressFile.write('USEARCH: Not done.' + '\n')
    progressFile.write('REMOVE NOT ALPHA: Not done.' + '\n')
    progressFile.write('SEARCH dsOLIGO: Not done.' + '\n')
    progressFile.write('SEARCH INDELS: Not done.' + '\n')


#----------------------------------------------------------------------------------------------------------------------------------------
# Check sentence in file
#----------------------------------------------------------------------------------------------------------------------------------------

def checkProgressFile(project_name, sentence):
    progressFile = open('./' + project_name + '/progress_file.txt', 'r')
    lines = progressFile.readlines()
    if sentence in lines:
        return(True)


#----------------------------------------------------------------------------------------------------------------------------------------
# Change status in progress file
#----------------------------------------------------------------------------------------------------------------------------------------

def changeProgressFile(project_name, sentence, new_sentence):
    file = open('./' + project_name + '/progress_file.txt', 'r')
    lines = file.read()
    lines = lines.replace(sentence, new_sentence)
    new_file = open('./' + project_name + '/progress_file.txt', 'w')
    new_file.write(lines)


#----------------------------------------------------------------------------------------------------------------------------------------
# Color code
#----------------------------------------------------------------------------------------------------------------------------------------

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


#----------------------------------------------------------------------------------------------------------------------------------------
# Append part to all the lists indicated
#----------------------------------------------------------------------------------------------------------------------------------------

def appendList(lists, part, keys):
    for element in lists:
        if type(element) is defaultdict: # if the element in the list is a defaultdict (sgRNAefficiency), append the part per hit of each sgRNA separately (keeped in keys)
            for k in keys: element[k].append(part)
        else: element.append(part)


#----------------------------------------------------------------------------------------------------------------------------------------
# Run command and print
#----------------------------------------------------------------------------------------------------------------------------------------

def running(command):
    print(f"{bcolors.OKGREEN}{command}{bcolors.ENDC}")
    os.system(command)


#----------------------------------------------------------------------------------------------------------------------------------------
# Write in temp jobs file
#----------------------------------------------------------------------------------------------------------------------------------------

def temp_jobs(project_name, command):
    file = open('./' + project_name + '/results_usearch/temp_jobs.txt', 'a')
    file.write(command + '\n')


#----------------------------------------------------------------------------------------------------------------------------------------
# Concat dataframes of stats of map and indels generated in indels.py
#----------------------------------------------------------------------------------------------------------------------------------------

def mergedfs(tmp_folder, suffix, column_index, output):
    files_list = [tmp_folder + file for file in os.listdir(tmp_folder) if file.endswith(suffix)]
    df_list = []
    for file in files_list:
        df_list.append(pd.read_csv(file, sep = "\t", header = 0))
    dfs = [df.set_index(column_index) for df in df_list]
    df_complete = pd.concat(dfs, axis=1)
    df_complete.to_csv(output, sep="\t")

#----------------------------------------------------------------------------------------------------------------------------------------
# Generate reads and freq per group of otus
#----------------------------------------------------------------------------------------------------------------------------------------

def calculateStats(otutable, nameFilter, otus_list, totalreads, sample):
    otutable[nameFilter] = ['yes' if otu in otus_list else 'no' for otu in otutable['#OTU ID']]
    reads = otutable.loc[otutable[nameFilter] == 'yes', sample].sum()
    if totalreads == 0.0:
        freq = 0.0
    else:
        freq = float((reads/totalreads)*100)
    return(reads, freq)


#----------------------------------------------------------------------------------------------------------------------------------------
# Generate fasta files per group of otus: NOT USED
#----------------------------------------------------------------------------------------------------------------------------------------

def generateFasta(denoised, nameOutput, otus_list):
    denoisedFasta = SeqIO.to_dict(SeqIO.parse(denoised, 'fasta')) # otus in denoised
    outputFasta = open(nameOutput, 'w+')
    for k, v in denoisedFasta.items():
        otu = k.split(";")[0]
        if otu in otus_list:
            outputFasta.write('>' + otu + '\n')
            outputFasta.write(str(v.seq) + '\n')


#----------------------------------------------------------------------------------------------------------------------------------------
# Generate a dict with items per number of lists in which this item appear
#----------------------------------------------------------------------------------------------------------------------------------------

def countIntersections(a):
    b = defaultdict(list)
    for item in nonRedundant([x for l in list(a.values()) for x in l]):
        counter = len([l for l in list(a.values()) if item in l])
        b[counter].append(item)
    return(b)


#----------------------------------------------------------------------------------------------------------------------------------------
# Generate non-redundant list
#----------------------------------------------------------------------------------------------------------------------------------------

def nonRedundant(listItem):
    return(list(set(listItem)))


#----------------------------------------------------------------------------------------------------------------------------------------
# Generate non-redundant list for more complex samples: item[0]
#----------------------------------------------------------------------------------------------------------------------------------------

def nonRedundantItem0(listItem):
    newlistItem = []
    for list_item in listItem:
        if list_item[0] not in [item[0] for item in newlistItem]:
            newlistItem.append(list_item)
    return(newlistItem)


#----------------------------------------------------------------------------------------------------------------------------------------
# Check if slurm process finished
#----------------------------------------------------------------------------------------------------------------------------------------

def checkrunning(job_name):
    print(job_name) # borrar
    count = 0
    check = 'Not completed'
    while True:
        if count > 172800: # 48 h maximum, denoising for some samples (protoplasts minampsize 1) takes many hours
            break
        try:
            time.sleep(2) # It takes few seconds in start the running
            status = int(subprocess.run('squeue -n ' + job_name + ' -o "%t" | grep -c "R"', shell=True, capture_output=True, text=True).stdout.strip())
        except:
            print('ERROR')
            print(int(subprocess.run('squeue -n ' + job_name + ' -o "%t" | grep -c "R"', shell=True, capture_output=True, text=True).stdout.strip()))
        print('Running processes: ' + str(status))
        if status == 0:
            check = 'Completed'
            break
        else:
            count += 20
            time.sleep(20)
    return(check)
