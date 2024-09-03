# Usearch pipeline

import os
import sys
from utils import changeProgressFile, bcolors, checkrunning
import argparse

parser = argparse.ArgumentParser(
                    prog = 'usearch.py',
                    description = 'Usearch methods')

parser.add_argument('-i', '--input',
                    help='<string> Path/to/input prefix, to select paired-end files. Alternative to -ip option.')
parser.add_argument('-f', '--multi_inputs',
                    help='If there are more than one sample to process, do not put prefix in -i option instead.')
parser.add_argument('-s', '--per_sample',
                    help='Run usearch per sample (parallelize).')
parser.add_argument('-p', '--project_name',
                    help='<string> Project name.')
parser.add_argument('-m', '--minampsize',
                    help='<int> Value of minampsize Usearch option of denoising.')

args = parser.parse_args()

PATHSCRIPTS="/home/tonin/mimi/protoplasts_alpha/Scripts/ampAnalysisv2/source"

print('\n'
        '#####################################################\n'
        '##                                                 ##\n'
        '##               USEARCH PIPELINE                  ##\n'
        '##                                                 ##\n'
        '#####################################################\n'
        '\n')
os.system('mkdir -p ./' + args.project_name + '/results_usearch')
if args.multi_inputs == 'False':
    print('\n'
            '******                INPUT                   ******'
            '\n')
    # For one sample: usearch_onesample.sh
    file_name = os.path.basename(args.input).split('.')[0]
    file=args.input + '*_R1*.fastq'
    try:
        os.system(PATHSCRIPTS + '/usearch_onesample.sh ' + file + ' ' + file_name + ' ' + args.minampsize + ' ' + args.project_name)
        changeProgressFile(args.project_name, 'USEARCH: Not done.\n', 'USEARCH: Done.\n')
    except:
        print('\n'
                f"{bcolors.FAIL}ERROR: USEARCH can not be completed.{bcolors.ENDC}"
                '\n')
        sys.exit()
else:
    print('\n'
        '******             INPUT FILES               ******'
        '\n')
    
    if args.per_sample == 'False':
        print('\n'
            f"{bcolors.OKBLUE}ADVISE: all the samples will be analyzed together by Usearch.{bcolors.ENDC}"
            '\n')
        try:
            os.system(PATHSCRIPTS + '/usearch_allsamples.sh ' + args.input + ' ' + args.minampsize + ' ' + args.project_name)
            changeProgressFile(args.project_name, 'USEARCH: Not done.\n', 'USEARCH: Done.\n')
        except:
            print('\n'
                f"{bcolors.FAIL}ERROR: USEARCH can not be completed.{bcolors.ENDC}"
                '\n')
            sys.exit()
    else:
        print('\n'
            f"{bcolors.OKBLUE}ADVISE: all the samples will be analyzed separately by Usearch.{bcolors.ENDC}"
            '\n')
        try:
            tempjobsPath = './' + args.project_name + '/results_usearch/array_script_usearch.slurm'
            if os.path.exists(tempjobsPath): os.remove(tempjobsPath)
            os.system(PATHSCRIPTS + '/generate_slurm_usearch.sh ' + args.project_name + ' ' + args.input + ' ' + args.minampsize)
            os.system('sbatch ./' + args.project_name + '/results_usearch/array_script_usearch.slurm')
            if checkrunning(args.project_name) != 'Completed': sys.exit() # if all the processes finished, continue
            changeProgressFile(args.project_name, 'USEARCH: Not done.\n', 'USEARCH: Done.\n')
        except:
            print('\n'
                f"{bcolors.FAIL}ERROR: USEARCH can not be completed.{bcolors.ENDC}"
                '\n')
            sys.exit()
