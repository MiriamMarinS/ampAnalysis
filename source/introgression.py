from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from analysisSeq import Otu
from progressbar import *
import pandas as pd
import sys
import os

parser = argparse.ArgumentParser()

# Arguments
parser.add_argument("-d", "--denoised")
parser.add_argument("-o", "--otu_table")
parser.add_argument("-i", "--input")
parser.add_argument("-c", "--complete_diffs")
parser.add_argument("-s", "--sample")
parser.add_argument("-m", "--metafile")

args = parser.parse_args()

widgets = ['Test: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),
           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

sample = args.sample

# Reading the otutable
print('Reading otutable information')
otutable = pd.read_csv(args.otu_table, sep="\t", header=0)
otutable = otutable[otutable[sample] > 0] # Filter if otu is present in this sample (in case of usearch_allsamples.sh)
otus_sample = otutable['#OTU ID'].to_list() # List of otus present in sample: > 0 reads mapped

# Reading metafile information
print("Reading the metafile information.")
metafile = pd.read_csv(args.metafile, sep="\t", header=None)
try:
   targets = [metafile[metafile[0].str.contains(sample)][1].values[0], metafile[metafile[0].str.contains(sample)][2].values[0]]
   print(sample + ' targets:')
   print(targets[0])
   print(targets[1])
except:
   print("ERROR: there are no sample information in the metafile, check the sample name in the file or maybe this sample does not have dsOligos.\n")
   sys.exit()

# Read fasta of denoised: filter otus by alpha-gliadins (if apply) and by sample
print("Reading the denoised file.")
otus_denoised = SeqIO.to_dict(SeqIO.parse(args.denoised, "fasta"))
otus_list = []
otus_names = []
for k, v in otus_denoised.items(): # only otus present in denoised file
    otu = k.split(";")[0].strip()
    if otu in otus_sample: # only otus present in sample
        otus_list.append(Otu(otu, v, sample))
        otus_names.append(otu)



# Searching the dsOligo in the denoised file
print("Searching the dsOligo in the denoised file.")
otus_shorttarget = []
otus_largetarget = []

# Percentage of diffs
if float(args.complete_diffs) <= 100.0:
   p_diffs = float(args.complete_diffs)
else:
   print("ERROR: the percentage of differences must be equal or lower than 100%.")
   sys.exit()

count = 0
pbar = ProgressBar(widgets=widgets, maxval=len(list(otus_denoised.keys())))
pbar.start()
for otu in otus_list:
   if otu.targetOtu(targets[1], 0) != "NA":
      otus_shorttarget.append(otu.otu_id) # 0, becouse the short target has to be exact.
      if otu.targetOtu(targets[0], p_diffs) != "NA": otus_largetarget.append(otu.otu_id)
   count += 1
   pbar.update(count)

pbar.finish()

# Statistics of dsOligo
print("Calculating statistics for the dsOligo...")

# Number of otus with targets, it is not depend of the alpha option
notusshort = len(otus_shorttarget)
notuslarge = len(otus_largetarget)

# Reads with targets, it depends of the alpha option
otutable['Status_denoised'] = ['denoised' if otu in otus_names else 'no' for otu in otutable['#OTU ID']] # Filter if otu is alpha-gliadin
otutable = otutable[otutable['Status_denoised'] == 'denoised']

totalreads = otutable[sample].sum() # sum of reads of all the otus
otutable['Status_short'] = ['short' if otu in otus_shorttarget else 'no' for otu in otutable['#OTU ID']]
otutable['Status_large'] = ['large' if otu in otus_largetarget else 'no' for otu in otutable['#OTU ID']]
readsotusshort = otutable.loc[otutable['Status_short'] =='short', sample].sum()
readsotuslarge = otutable.loc[otutable['Status_large'] =='large', sample].sum()

# Percentage of reads with targets, it depends of the alpha option
potusshort = float((readsotusshort/totalreads))*100
potuslarge = float((readsotuslarge/totalreads))*100

print([notusshort, notuslarge, readsotusshort, readsotuslarge, potusshort, potuslarge])

# Append to existing stats
if os.stat(args.input).st_size == 0:
   stats = pd.DataFrame({'Params': ['Notus_short', 'Notus_large', 'Reads_short', 'Reads_large', 'Freq_short', 'Freq_large']})
else:
   stats = pd.read_csv(args.input, sep = "\t", header=0)
stats_temp = pd.DataFrame({sample: [notusshort, notuslarge, readsotusshort, readsotuslarge, potusshort, potuslarge]})
stats = pd.concat([stats.reset_index(drop=True), stats_temp], axis=1)

stats.to_csv(args.input, sep="\t", index=False)
print(sample + ': DONE.')
