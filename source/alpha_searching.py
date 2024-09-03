import os
import argparse
import sys
from utils import changeProgressFile, bcolors, checkrunning


parser = argparse.ArgumentParser(
                    prog = 'stats_usearch.py',
                    description = 'Usearch stats methods')

parser.add_argument('-p', '--project_name',
                    help='<string> Project name.')

args = parser.parse_args()

PATHSCRIPTS="/home/tonin/mimi/protoplasts_alpha/Scripts/ampAnalysisv2/source"

print('\n'
        '#####################################################\n'
        '##                                                 ##\n'
        '##          REMOVING NOT ALPHA AMPLICONS           ##\n'
        '##                                                 ##\n'
        '#####################################################\n'
        '\n')

path='./' + args.project_name + '/results_usearch/'
pathmmseqs2 = './' + args.project_name + '/mmseqs2'
denoisedFiles = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) if f.endswith('_denoised.fa')]
pathalphaDB = '/home/tonin/mimi/database/alpha_gliadinmmseqs2/alphaDB'
pathalphatemp = pathmmseqs2 + '/tmp'

statsmmseqs2otusFile = pathmmseqs2 + '/mmseqs2_stats.txt'
statsmmseqs2readsFile = pathmmseqs2 + '/reads_stats.txt'
if os.path.exists(statsmmseqs2otusFile): os.remove(statsmmseqs2otusFile)
if os.path.exists(statsmmseqs2readsFile): os.remove(statsmmseqs2readsFile)

try:
    print('Amplicons against alpha DB...')
    tempjobsPath = './' + args.project_name + '/mmseqs2/array_script_mmseqs2.slurm'
    if os.path.exists(tempjobsPath): os.remove(tempjobsPath)
    os.system(PATHSCRIPTS + '/generate_slurm_mmseqs2.sh ' + args.project_name + ' ' + path + ' ' + pathalphaDB + ' ' + pathalphatemp)
    os.system('sbatch ./' + args.project_name + '/mmseqs2/array_script_mmseqs2.slurm') 
    if checkrunning(args.project_name) != 'Completed': sys.exit() # if all the processes finished, continue
    for filename in denoisedFiles:
        file = path + filename
        sample = filename.replace('_denoised.fa', '')
        print('Removing aplha from: ' + filename)
        print('DEBUG: ' + 'python ' + PATHSCRIPTS + '/alpha_remove.py -f ./' + args.project_name + '/mmseqs2/' + sample + '_mmseqs2.txt -d ./' + args.project_name + '/results_usearch/' + sample + '_denoised.fa -o ./' + args.project_name + '/alpha_removing/' + sample + '_denoised.fa -s ./' + args.project_name + '/mmseqs2/mmseqs2_stats.txt')
        os.system('python ' + PATHSCRIPTS + '/alpha_remove.py -f ./' + args.project_name + '/mmseqs2/' + sample + '_mmseqs2.txt -d ./' + args.project_name + '/results_usearch/' + sample + '_denoised.fa -o ./' + args.project_name + '/alpha_removing/' + sample + '_denoised.fa -s ./' + args.project_name + '/mmseqs2/mmseqs2_stats.txt')
        print('DEBUG: ' + 'python ' + PATHSCRIPTS + '/reads_alpha.py -d ./' + args.project_name + '/alpha_removing/' + sample + '_denoised.fa -o ./' + args.project_name + '/results_usearch/otutable_' + sample + '.txt -s ./' + args.project_name + '/mmseqs2/reads_stats.txt')
        os.system('python ' + PATHSCRIPTS + '/reads_alpha.py -d ./' + args.project_name + '/alpha_removing/' + sample + '_denoised.fa -o ./' + args.project_name + '/results_usearch/otutable_' + sample + '.txt -s ./' + args.project_name + '/mmseqs2/reads_stats.txt')
    changeProgressFile(args.project_name, 'REMOVE NOT ALPHA: Not done.\n', 'REMOVE NOT ALPHA: Done.\n')
except:
    print('\n'
            f"{bcolors.FAIL}ERROR: NOT ALPHA REMOVING process can not be completed.{bcolors.ENDC}"
            '\n')
    sys.exit()

