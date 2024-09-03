# Amplicon analysis: complete pipeline
# Updates: no truncate primers, search exact, local+global alignments, DBPC+DBAm + sgRNA+PAM
# Miriam Marin Sanz, 2023

############################################################
# GLOBAL VARIABLES                                         #
############################################################
PATHSCRIPTS="/home/tonin/mimi/protoplasts_alpha/Scripts/ampAnalysisv2/source"

############################################################
# Libraries                                                #
############################################################
import os
import sys
import argparse
from source.utils import createProgressFile, checkProgressFile, bcolors, running

############################################################
# Options                                                  #
############################################################
parser = argparse.ArgumentParser(
                    prog = 'ampAnalysis.py',
                    description = 'Analysis of InDels from NGS sequence')

parser.add_argument('-i', '--input',
                    help='<string> Path/to/input prefix, to select paired-end files. Alternative to -ip option.',
                    required=True)
parser.add_argument('-f', '--multi_inputs',
                    help='If there are more than one sample to process, do not put prefix in -i option instead.',
                    action='store_true',
                    required=False)
parser.add_argument('-s', '--per_sample',
                    help='Run usearch per sample (parallelize).',
                    action='store_true',
                    required=False)
parser.add_argument('-p', '--project_name',
                    help='<string> Project name.',
                    required=True)
parser.add_argument('-t', '--metafile',
                    help='<string> Metafile with the dsOligo information per sample. Required if --dsOligo option is used.',
                    required=False)
parser.add_argument('-m', '--minampsize',
                    help='<int> Value of minampsize Usearch option of denoising.',
                    required=True)
parser.add_argument('-a', '--alpha',
                    action='store_true',
                    help='Remove denoised amplicons not annotated as alpha-gliadins.',
                    required=False)
parser.add_argument('-d', '--dsOligo',
                    action='store_true',
                    help='Search the dsOligo (both) in denoised amplicons.',
                    required=False)
parser.add_argument('-c', '--num_diff',
                    help='<int> Value of % diffs in complete target searching process. Default: 0 (no diffs).',
                    required=False)
parser.add_argument('-n', '--indels',
                    action='store_true',
                    help='Search InDels in denoised amplicons.',
                    required=False)

args = parser.parse_args()


def main():
    #####################################################
    # PROGRESS FILE                                     #
    #####################################################
    os.system('mkdir -p ./' + args.project_name) # create folder of project
    if not os.path.isfile('./' + args.project_name + '/progress_file.txt'): # check if progress file exists
        createProgressFile(args.project_name)
    
    #####################################################
    # USEARCH PIPELINE                                  #
    #####################################################
    if checkProgressFile(args.project_name, 'USEARCH: Not done.\n'):
        if args.multi_inputs == True:
            multisample = 'True'
        else:
            multisample = 'False'
        if args.per_sample == True:
            per_sample = 'True'
        else:
            per_sample = 'False'
        print(args.multi_inputs)
        running('python ' + PATHSCRIPTS + '/usearch.py -i ' + args.input + ' -f ' + multisample + ' -s ' + per_sample + ' -p ' + args.project_name + ' -m ' + args.minampsize)
        running('python ' + PATHSCRIPTS + '/stats_usearch.py -p ' + args.project_name + ' -m ' + args.minampsize)
    
    #####################################################
    # REMOVING NOT ALPHA OTUS                           #
    #####################################################
    if args.alpha == True:
        if checkProgressFile(args.project_name, 'REMOVE NOT ALPHA: Not done.\n'):
            os.system('mkdir -p ./' + args.project_name + '/mmseqs2')
            os.system('mkdir -p ./' + args.project_name + '/alpha_removing')
            os.system('python ' + PATHSCRIPTS + '/alpha_searching.py -p ' + args.project_name)
            os.system('gzip -f ./' + args.project_name + '/mmseqs2/*_mmseqs2.txt')
    
    #####################################################
    # SEARCHING DSOLIGO                                 #
    #####################################################
    dsoligo_status = 'False'
    num_diff = 0 # default number of differences in complete dsOligo (not short target)
    if args.dsOligo == True:
        if args.metafile is None:
            print('\n'
                f"{bcolors.FAIL}ERROR: for --dsOligo option, --metafile is required.{bcolors.ENDC}"
                '\n')
            sys.exit()
        if checkProgressFile(args.project_name, 'SEARCH dsOLIGO: Not done.\n'):
            #os.system('mkdir -p ./' + args.project_name + '/dsOligo/')
            dsoligo_status = 'True'
            if args.num_diff is not None:
                num_diff = args.num_diff       

    #####################################################
    # SEARCHING INDELS                                  #
    #####################################################
    if args.indels == True:
        if args.metafile is None:
            print('\n'
                f"{bcolors.FAIL}ERROR: for --indels option, --metafile is required.{bcolors.ENDC}"
                '\n')
            sys.exit()
        if checkProgressFile(args.project_name, 'SEARCH INDELS: Not done.\n'):
            os.system('mkdir -p ./' + args.project_name + '/indels/')
            os.system('mkdir -p ./' + args.project_name + '/indels/stats_per_ampType/')
            if args.alpha == True:
                print('Searching InDels in alpha-gliadin filtered otus...')
                running('python ' + PATHSCRIPTS + '/indels_searching.py -p ' + args.project_name + ' -a ' + 'True -m ' + args.metafile + ' -g ' + dsoligo_status + ' -c ' + str(num_diff))
            else:
                print('Searching InDels in non-filtered otus...')
                running('python ' + PATHSCRIPTS + '/indels_searching.py -p ' + args.project_name + ' -a ' + 'False -m ' + args.metafile + ' -g ' + dsoligo_status + ' -c ' + str(num_diff))
        


if __name__ == "__main__":
    main()