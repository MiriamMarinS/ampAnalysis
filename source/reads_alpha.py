from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from progressbar import *
import pandas as pd
import os

parser = argparse.ArgumentParser()

# Arguments
parser.add_argument("-d", "--denoised_alpha")
parser.add_argument("-o", "--otu_table")
parser.add_argument("-s", "--stats")

args = parser.parse_args()

widgets = ['Test: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),
           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

# Read the otutable
print("Reading the otutable information.")
otutable = pd.read_csv(args.otu_table, sep="\t", header=0)
project = os.path.basename(args.denoised_alpha).split('/')[-1].split("_")[0].strip()

# Read denoised file with only alpha-gliadin otus
otus_alpha = [otu.split(";")[0].strip() for otu in list(SeqIO.to_dict(SeqIO.parse(args.denoised_alpha, "fasta")).keys())]

# Total reads and alpha-gliadin reads per sample
# Write stats
statsfile = open(args.stats, "a")
samples = otutable.columns.to_list()[1:]
for sample in samples:
    totalreads = otutable[sample].sum() # sum of reads of all the otus
    temp = otutable
    temp['Status'] = ['alpha' if otu in otus_alpha else 'other' for otu in temp['#OTU ID']]
    alphareads = otutable.loc[temp['Status'] =='alpha', sample].sum()
    freq = float(alphareads/totalreads)*100
    statsfile.write(sample + ": " + str(freq) + "% of reads" + "\n")
    print(sample + ': DONE.')
