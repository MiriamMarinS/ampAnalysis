# Searching dsOligo

import argparse
import os
import sys
from utils import changeProgressFile, bcolors, running
import pandas as pd

parser = argparse.ArgumentParser()

# Arguments
parser.add_argument("-p", "--project_name")
parser.add_argument("-a", "--alpha")
parser.add_argument("-t", "--metafile")
parser.add_argument("-c", "--com_target_diff")

args = parser.parse_args()

PATHSCRIPTS="/home/tonin/mimi/protoplasts_alpha/Scripts/ampAnalysisv2/source"

print('\n'
        '#####################################################\n'
        '##                                                 ##\n'
        '##                      dsOLIGO                    ##\n'
        '##                                                 ##\n'
        '#####################################################\n'
        '\n')

path = './' + args.project_name + '/dsOligo/'
statsFile = path + 'stats_introgression.txt'
if os.path.exists(statsFile): os.remove(statsFile)
os.system('touch ' + statsFile)

if args.com_target_diff == '0':
    com_target_diff = 0
else:
    com_target_diff = int(args.com_target_diff)

if args.alpha == 'True':
    pathdenoised = './' + args.project_name + '/alpha_removing/'
else:
    pathdenoised = './' + args.project_name + '/results_usearch/'

denoisedFiles = [f for f in os.listdir(pathdenoised) if os.path.isfile(os.path.join(pathdenoised, f)) if f.endswith('_denoised.fa')]

try:
    for filename in denoisedFiles:
        file = pathdenoised + filename
        print('Searching dsOligo in: ' + file)
        print('Complete dsOligo will be searched with ' + str(com_target_diff) + '% of differences.')
        project = filename.replace('_denoised.fa', '').strip()

        # Samples
        otutable = pd.read_csv('./' + args.project_name + '/results_usearch/otutable_' + project + '.txt', sep='\t', header=0)
        samples = otutable.columns.to_list()[1:]
        for sample in samples:
            print(sample)
            running('python ' + PATHSCRIPTS + '/introgression.py -d ' + file + ' -o ./' + args.project_name + '/results_usearch/otutable_' + project + '.txt -i ' + args.project_name + '/dsOligo/stats_introgression.txt -c ' + str(com_target_diff) + ' -s ' + sample + ' -m ' + args.metafile)
    changeProgressFile(args.project_name, 'SEARCH dsOLIGO: Not done.\n', 'SEARCH dsOLIGO: Done.\n')
except:
    print('\n'
            f"{bcolors.FAIL}ERROR: dsOLIGO SEARCHING can not be completed.{bcolors.ENDC}"
            '\n')
    sys.exit()