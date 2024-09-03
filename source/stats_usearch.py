import os
import argparse
import pandas as pd
from utils import bcolors


parser = argparse.ArgumentParser(
                    prog = 'stats_usearch.py',
                    description = 'Usearch stats methods')

parser.add_argument('-p', '--project_name',
                    help='<string> Project name.')
parser.add_argument('-m', '--minampsize',
                    help='<int> Value of minampsize Usearch option of denoising.')

args = parser.parse_args()

print('\n'
        '#####################################################\n'
        '##                                                 ##\n'
        '##                   STATS USEARCH                 ##\n'
        '##                                                 ##\n'
        '#####################################################\n'
        '\n')

# log files
path='./' + args.project_name + '/results_usearch/'
logFiles = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) if f.endswith('_log.txt')]

statsFile = path + 'stats_usearch.txt'
if os.path.exists(statsFile): os.remove(statsFile)

stats_samples = {}
for filename in logFiles:
    print(filename)
    try:
        file = open(path + filename, 'r')
        lines = file.readlines()
        sample = filename.replace('_log.txt', '')
        denoised_list = []
        for line in lines:
            if 'Merged (' in line:
                merged = float(line.split(',')[1].replace('%)', '').strip())
            if 'Filtered' in line:
                filtered = float(line.split('(')[1].split(',')[1].replace('%)', ''))
            if 'singletons' in line:
                uniques = int(line.split(',')[1].replace(' uniques', '').strip())
                singletons = float(line.split(',')[2].split('(')[1].replace('%)', '').strip())
            if 'good' in line and '100.0%' in line:
                denoised_list.append(int(line.split('%')[1].split(',')[0].replace(' good', '').strip()))
            if 'mapped to OTUs' in line:
                searched = float(line.split('(')[1].replace('%)', '').strip())
        denoised = denoised_list[-1]
        stats_samples[sample + 'minampsize = ' + str(args.minampsize)] = [merged, filtered, uniques, singletons, denoised, searched]
    except:
        print('\n'
                f"{bcolors.FAIL}ERROR: STATS OF USEARCH can not be completed.{bcolors.ENDC}"
                '\n')
        #sys.exit()
        continue

stats = pd.DataFrame.from_dict(stats_samples)
stats.index = ['Merge %', 'Filtered %', 'N uniques', 'Singletons %', 'N denoised', 'Searched %']

stats.to_csv(path + 'stats_usearch.txt', sep="\t")
print(sample + ': DONE.')
