# Remove not alpha-gliadin denoised amplicons from the fasta file based on blast results

from Bio import SeqIO
import argparse
import pandas as pd
from progressbar import *
import os

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', help='/path/to/mmseqs2 descript output')
parser.add_argument('-d', '--denoised', help='/path/to/denoised input')
parser.add_argument('-o', '--output', help='/path/to/output denoised new file')
parser.add_argument('-s', '--stats', help='/path/to/stats file')
args = parser.parse_args()

widgets = ['Test: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),
           ' ', ETA(), ' ', FileTransferSpeed()]

mmseqs2Output = pd.read_csv(args.file, sep = "\t", header=None)

# Number of denoised amplicons
ndenoised = len([record for record in SeqIO.parse(args.denoised, "fasta")])

# Number of alpha-gliadin otus
nalphadenoised = len(list(set(mmseqs2Output[0].values.tolist())))

# Percentage of alpha-gliadins
palpha = (float(nalphadenoised/ndenoised))*100

# Otus of alpha-gliadins
conservedOtus = list(otuid for otuid in list(set(mmseqs2Output[0].values.tolist())) if "Otu" in otuid)

# Write new denoised file with alpha-gliadins otus
newdenoised = open(args.output, "w+")

count = 0
pbar = ProgressBar(widgets=widgets, maxval=len([record.id for record in SeqIO.parse(args.denoised, "fasta")]))
pbar.start()
for record in SeqIO.parse(args.denoised, "fasta"):
   if str(record.id) in conservedOtus:
      SeqIO.write(record, newdenoised, "fasta")
   count += 1
   pbar.update(count)

pbar.finish()

# Write stats
statsfile = open(args.stats, "a")
statsfile.write(os.path.basename(args.denoised).split('/')[-1].split("_")[0] + ": " + str(palpha) + "% of otus" + "\n")
