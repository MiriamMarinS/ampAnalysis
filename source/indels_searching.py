# Searching InDels

import argparse
import os
import sys
from utils import changeProgressFile, bcolors, running, checkrunning, mergedfs
import pandas as pd

parser = argparse.ArgumentParser()

# Arguments
parser.add_argument("-p", "--project_name")
parser.add_argument("-a", "--alpha")
parser.add_argument("-m", "--metafile")
parser.add_argument("-g", "--dsoligo")
parser.add_argument("-c", "--num_diff")

args = parser.parse_args()

PATHSCRIPTS="/home/tonin/mimi/protoplasts_alpha/Scripts/ampAnalysisv2/source"

print('\n'
        '#####################################################\n'
        '##                                                 ##\n'
        '##                       INDELS                    ##\n'
        '##                                                 ##\n'
        '#####################################################\n'
        '\n')

path = './' + args.project_name + '/indels/'
statsFile = path + 'stats_indels.txt'
statsmapped = path + 'stats_mapped.txt'
statsamptype_path = path + 'stats_per_ampType/'
statsamptype_files = os.listdir(statsamptype_path)
if os.path.exists(statsFile): os.remove(statsFile)
if os.path.exists(statsmapped): os.remove(statsmapped)
if len(statsamptype_files) != 0: # remove stats files of amp types if they exist
    for file in statsamptype_files:
        if file.endswith('_stats.txt'):
            os.remove(statsamptype_path + file)
os.system('touch ' + statsFile)
os.system('touch ' + statsmapped)

if args.alpha == 'True':
    pathdenoised = './' + args.project_name + '/alpha_removing/'
else:
    pathdenoised = './' + args.project_name + '/results_usearch/'

denoisedFiles = [f for f in os.listdir(pathdenoised) if os.path.isfile(os.path.join(pathdenoised, f)) if f.endswith('_denoised.fa')]

bwa_path = '/home/tonin/mimi/database/BW208_Amp/bwa/BW208_Amp' # ORIGINAL
bbmap_path = '/home/tonin/mimi/database/BW208_Amp/bbmap/BW208_Amp.fasta' # ORIGINAL
#bwa_path = '/home/tonin/mimi/database/Back_T590+BW208/bwa/Back_T590+BW208' # DH
#bbmap_path = '/home/tonin/mimi/database/Back_T590+BW208/bbmap/Back_T590+BW208.fasta' # DH
#bwa_path = '/home/tonin/mimi/database/Back_AN+BW208/bwa/Back_AN+BW208' # DH
#bbmap_path = '/home/tonin/mimi/database/Back_AN+BW208/bbmap/Back_AN+BW208.fasta' # DH


try:
    # Align amplicons to references with BWA-mem and BBmap alignment tools
    tempjobsPath = './' + args.project_name + '/indels/array_script_align.slurm'
    if os.path.exists(tempjobsPath): os.remove(tempjobsPath)
    os.system(PATHSCRIPTS + '/generate_slurm_align.sh ' + args.project_name + ' ' + pathdenoised + ' ' + bwa_path + ' ' + bbmap_path)
    os.system('sbatch ./' + args.project_name + '/indels/array_script_align.slurm') 
    if checkrunning(args.project_name) != 'Completed': sys.exit() # if all the processes finished, continue
    
    # Search indels
    tempjobsPath = './' + args.project_name + '/indels/array_script_indels.slurm'
    if os.path.exists(tempjobsPath): os.remove(tempjobsPath)
    os.system(PATHSCRIPTS + '/generate_slurm_indels.sh ' + args.project_name + ' ' + path + ' ' + pathdenoised + ' ' + args.metafile + ' ' + args.dsoligo + ' ' + args.num_diff)
    os.system('sbatch ./' + args.project_name + '/indels/array_script_indels.slurm') 
    if checkrunning(args.project_name) != 'Completed': sys.exit() # if all the processes finished, continue
    # Merge stats files
    mergedfs('./' + args.project_name + '/indels/temp_stats_map/', '_statsmapped.txt', 'Params', path + 'stats_mapped.txt') # merge stats map
    mergedfs('./' + args.project_name + '/indels/temp_stats_indels/', '_statsindels.txt', 'Params', path + 'stats_indels.txt') # merge stats indels
    mergedfs('./' + args.project_name + '/indels/temp_stats_indels/', '_missingEpitopes.txt', 'Epitopes', path + 'stats_missingEpitopes.txt') # merge stats missing epitopes
    mergedfs('./' + args.project_name + '/indels/temp_stats_indels/', '_epitopes.txt', 'Epitopes', path + 'stats_epitopes.txt') # merge stats epitopes
    mergedfs('./' + args.project_name + '/indels/temp_stats_indels/', '_lengthindels.txt', 'Events', path + 'stats_lengthindels.txt') # merge stats length indel events
    mergedfs('./' + args.project_name + '/indels/temp_stats_indels/', '_pamefficiency.txt', 'PAM+guide', path + 'stats_pamefficiency.txt') # merge stats pam efficiency
    if args.dsoligo == 'True':
        mergedfs('./' + args.project_name + '/indels/temp_stats_indels/', '_statsnumbersgRNAs.txt', 'Params', path + 'stats_dsoligo_NgRNAs_rname.txt') # merge stats dsoligo per number of sgRNAs in rname, calculate the probability of insertion by the number of sgRNAs hits in the reference
        mergedfs('./' + args.project_name + '/indels/temp_stats_indels/', '_statsODNinNumbersgRNAs.txt', 'Params', path + 'stats_ODNInNumbersgRNAs.txt') # merge stats dsoligo per number of sgRNAs with the insertion, calculate the frequency of insertion in x number of sgRNAs (between 2 sgRNAs, in 1, 3...)
    type_tmpFolder = './' + args.project_name + '/indels/temp_stats_indels/'
    amp_types = [file.split('_')[1] for file in os.listdir(type_tmpFolder) if file.endswith('_typestats.txt')]
    for amp_type in amp_types:
        mergedfs(type_tmpFolder, amp_type + '_typestats.txt', 'Params', path + 'stats_per_ampType/' + amp_type + '_stats.txt') # merge stats per amp type
    amp_types_number = [file.split('_')[1] for file in os.listdir(type_tmpFolder) if file.endswith('_typenumberstats.txt')]
    for amp_type_number in amp_types_number:
        mergedfs(type_tmpFolder, amp_type_number + '_typenumberstats.txt', 'Params', path + 'stats_per_ampType/' + amp_type_number + '_stats.txt') # merge stats per number of sgRNAs-amp type
    changeProgressFile(args.project_name, 'SEARCH INDELS: Not done.\n', 'SEARCH INDELS: Done.\n')
    if args.dsoligo == 'True': changeProgressFile(args.project_name, 'SEARCH dsOLIGO: Not done.\n', 'SEARCH dsOLIGO: Done.\n')
except:
    print('\n'
            f"{bcolors.FAIL}ERROR: INDELS SEARCHING can not be completed.{bcolors.ENDC}"
            '\n')
    sys.exit()
