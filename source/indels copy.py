from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from progressbar import *
import pandas as pd
from analysisSeq import analysisSamv2, epitopeFind, sgRNA
from utils import calculateStats
import os
from functools import reduce
from collections import defaultdict
import itertools

parser = argparse.ArgumentParser()

# Arguments
parser.add_argument("-s", "--sam_bwa")
parser.add_argument("-b", "--sam_bbmap")
parser.add_argument("-o", "--otu_table")
parser.add_argument("-t", "--stats")
parser.add_argument("-m", "--stats_mapped")
parser.add_argument("-a", "--sample")
parser.add_argument("-d", "--denoised")
parser.add_argument("-f", "--metafile")
parser.add_argument("-g", "--dsoligo")
parser.add_argument("-c", "--num_diff")

args = parser.parse_args()

widgets = ['Test: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),
           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options


#------------------------------------------------------------------------------------------------------
# Params
#------------------------------------------------------------------------------------------------------
# Get the sample name
sample =  args.sample

# Get the status of dsOligo
dsOligo_status = args.dsoligo

# Fasta with references
refsFasta = '/home/tonin/mimi/database/BW208_Amp/BW208_Amp.fasta'


#------------------------------------------------------------------------------------------------------
# Search the epitopes in the Amp database
#------------------------------------------------------------------------------------------------------
amp_types, hits_epitopes = epitopeFind('/home/tonin/mimi/database/epitopes/epitopes.fasta', refsFasta)


#------------------------------------------------------------------------------------------------------
# Search Indels oand/or dsoligo
#------------------------------------------------------------------------------------------------------
# Read SAM file from bwa-mem and BBmap, filter by otus in sample
otutable = pd.read_csv(args.otu_table, sep="\t", header=0)
otutable = otutable[otutable[sample] > 0] # Filter otus in sample
otus_sample = otutable['#OTU ID'].to_list() # otus in sample

# Search indels/dsoligo
bwotus_wt, bwotus_indels, bwotus_extraction, bwotus_all, bwundetermined_otus, bwotus_sgRNAs, bwotus_per_ref, bwotus_dsoligo = analysisSamv2(args.sam_bwa, otus_sample, args.metafile, refsFasta, sample, hits_epitopes, dsOligo_status, int(args.num_diff))
bbotus_wt, bbotus_indels, bbotus_extraction, bbotus_all, bbundetermined_otus, bbotus_sgRNAs, bbotus_per_ref, bbotus_dsoligo = analysisSamv2(args.sam_bbmap, otus_sample, args.metafile, refsFasta, sample, hits_epitopes, dsOligo_status, int(args.num_diff))


#------------------------------------------------------------------------------------------------------
# Process otus per category
#------------------------------------------------------------------------------------------------------
# Combined indels/dsoligo
bwotus_indels_values = list(set(reduce(lambda i,j:i+j,list(bwotus_indels.values()),[])))
bbotus_indels_values = list(set(reduce(lambda i,j:i+j,list(bbotus_indels.values()),[])))
otus_indels = list(set(bwotus_indels_values + bbotus_indels_values)) # combined list (non-redundant) of otus with indels in BWA-mem or BBmap
otus_dels = list(set(bwotus_indels['Dels'] + bbotus_indels['Dels']))
otus_ins = list(set(bwotus_indels['Ins'] + bbotus_indels['Ins']))
otus_wt = [otu for otu in list(set(bwotus_wt + bbotus_wt)) if otu not in otus_indels] # combined list (non-redundant) of otus without indels in BW-mem or BBmap
otus_undeterminedtemp = list(set(bwundetermined_otus + bbundetermined_otus))
otus_undetermined = [otu for otu in otus_undeterminedtemp if otu not in otus_indels and otu not in otus_wt]
bwotus_dsoligo_values = list(set(reduce(lambda i,j:i+j,list(bwotus_dsoligo.values()),[])))
bbotus_dsoligo_values = list(set(reduce(lambda i,j:i+j,list(bbotus_dsoligo.values()),[])))
otus_dsoligo = list(set(bwotus_dsoligo_values + bbotus_dsoligo_values)) # combined list (non-redundant) of otus with dsoligo in BWA-mem or BBmap

# Combined indels per sgRNA
otus_sgRNAs = {}
for k, v in bwotus_sgRNAs.items():
   otus_sgRNAs[k] = list(set(v + bbotus_sgRNAs[k]))
otus_indels_allsgRNAs = list(set.intersection(*map(set,list(otus_sgRNAs.values())))) # otus with indels in all the sgRNAs
count = 0
while len(otus_sgRNAs) < 3:
   otus_sgRNAs[count] = 'NA'
   count += 1

# Combined extractions
otus_extraction = defaultdict(list)
for amptype in list(amp_types.keys()):
   for element in list(set(bwotus_extraction[amptype] + bbotus_extraction[amptype])):
      otus_extraction[amptype].append(element)

# Combined dsoligo per amp type/rname
if dsOligo_status == 'True':
   otus_dsoligo_rname = defaultdict(list)
   otus_dsoligo_amp = defaultdict(list)
   sgRNAs_number = {} # number of sgRNAs per Amp (rname)
   metaFile = pd.read_csv(args.metafile, sep="\t", header=None)
   plasmid = metaFile[metaFile[0] == sample][3].values[0]
   refs = SeqIO.to_dict(SeqIO.parse(refsFasta, 'fasta'))
   for rname in list(amp_types.values()):
      for element in list(set(bwotus_dsoligo[rname] + bbotus_dsoligo[rname])):
         otus_dsoligo_rname[rname].append(element)
         sgRNAs_number[rname] = len(sgRNA(str(refs[rname].seq), plasmid))
         otus_dsoligo_amp[[i for i in amp_types if rname in amp_types[i]][0]].append(element)


#------------------------------------------------------------------------------------------------------
# Stats of dsoligo per number of sgRNAs
#------------------------------------------------------------------------------------------------------
if dsOligo_status == 'True':
   print('Generating stats of dsoligo per number of sgRNAs...')
   stats_dsoligo_number_sgRNAs = {}
   for number in list(set(sgRNAs_number.values())):
      rnames_ = [i for i in sgRNAs_number if sgRNAs_number[i] == number] # all rnames with x number of sgRNAs
      
      otus_NsgRNAstotal = list(set(itertools.chain(*[bwotus_per_ref[rname] + bbotus_per_ref[rname] for rname in rnames_]))) # otus aligned to rname with x number of sgRNAs
      otus_dsoligoNsgRNAs = list(set(itertools.chain(*[v for k, v in otus_dsoligo_rname.items() if k in rnames_]))) # otus aligned to rname with x number of sgRNAs and dsoligo

      otutable_stats = otutable.copy() # only contain otus in sample (and denoised if apply (alpha option))

      # Calculate total reads
      otutable_stats['Status_NsgRNA'] = ['NsgRNA' if otu in otus_NsgRNAstotal else 'no' for otu in otutable_stats['#OTU ID']]
      totalreads = otutable_stats[otutable_stats['Status_NsgRNA'] == 'NsgRNA'].sum() # sum of reads of all otus by rname

      reads_NsgRNA, freq_NsgRNA = calculateStats(otutable_stats, 'Status_dsoligo_NsgRNA', otus_dsoligoNsgRNAs, totalreads)
      #otutable_stats['Status_dsoligo_NsgRNA'] = ['dsoligo' if otu in otus_dsoligoNsgRNAs else 'no' for otu in otutable_stats['#OTU ID']]
      #reads_NsgRNA = otutable_stats[otutable_stats['Status_dsoligo_NsgRNA'] == 'dsoligo'].sum()
      #freq_NsgRNA = float((reads_NsgRNA/totalreads)*100)

      stats_dsoligo_number_sgRNAs['otus_ODN_NsgRNAs' + str(number)] = len(otus_dsoligoNsgRNAs)
      stats_dsoligo_number_sgRNAs['reads_ODN_NsgRNAs' + str(number)] = reads_NsgRNA
      stats_dsoligo_number_sgRNAs['freq_ODN_NsgRNAs' + str(number)] = freq_NsgRNA

   stats_dsoligo_pernumber = pd.DataFrame({'Params': ['otus_ODN_NsgRNAs' + str(item) for item in list(set(sgRNAs_number.values()))] +
                                          ['reads_ODN_NsgRNAs' + str(item) for item in list(set(sgRNAs_number.values()))] +
                                          ['freq_ODN_NsgRNAs' + str(item) for item in list(set(sgRNAs_number.values()))], 
                                          sample: [stats_dsoligo_number_sgRNAs['otus_ODN_NsgRNAs' + str(item)] for item in list(set(sgRNAs_number.values()))] +
                                                   [stats_dsoligo_number_sgRNAs['reads_ODN_NsgRNAs' + str(item)] for item in list(set(sgRNAs_number.values()))] +
                                                   [stats_dsoligo_number_sgRNAs['freq_ODN_NsgRNAs' + str(item)] for item in list(set(sgRNAs_number.values()))]})
   stats_dsoligo_pernumber.to_csv(args.stats.split("_")[0] + '_statsnumbersgRNAs.txt', sep="\t", index=False)


#------------------------------------------------------------------------------------------------------
# Stats mapped, filter by sample, by alpha if apply
#------------------------------------------------------------------------------------------------------
print('Generating stats of mapping...')
otus_denoised = [k.split(";")[0].strip() for k in SeqIO.to_dict(SeqIO.parse(args.denoised, 'fasta')).keys()] # otus in denoised
all_otus = [otu for otu in otus_sample if otu in otus_denoised] # otus in sample and denoise (filter by alpha or not)
freq_mappedbwa = float((len(bwotus_all)/len(all_otus))*100)
freq_mappedbbmap = float((len(bbotus_all)/len(all_otus))*100)
# Generate stats mapped
stats_mapped = pd.DataFrame({'Params': ['Freq otus mapped bwa-mem', 'Freq otus mapped bbmap'], sample: [freq_mappedbwa, freq_mappedbbmap]})
stats_mapped.to_csv(args.stats_mapped, sep="\t", index=False)


#------------------------------------------------------------------------------------------------------
# Generate principal stats: % of indels, % of indels per sgRNA, % of dsoligo if apply
#------------------------------------------------------------------------------------------------------
print('Generating stats of indels/dsoligo...')
# BWA-mem and BBmap: filter by sample, by alpha if apply, and by mapped
otutable_stats = otutable.copy() # only contain otus in sample (and denoised if apply (alpha option))
otutable_stats['Status_mapped'] = ['mapped' if otu in list(set(otus_indels + otus_wt + otus_undetermined)) else 'no' for otu in otutable_stats['#OTU ID']]
otutable_stats = otutable_stats[otutable_stats['Status_mapped'] == 'mapped'] # filter otutable by otus mapped with bwa-mem or BBmap

totalreads = otutable_stats[sample].sum() # sum of reads of all the otus in sample mapped

reads_wt, freq_wt = calculateStats(otutable_stats, 'Status_wt', otus_wt, totalreads)
reads_indels, freq_indels = calculateStats(otutable_stats, 'Status_indels', otus_indels, totalreads)
reads_dels, freq_dels = calculateStats(otutable_stats, 'Status_dels', otus_dels, totalreads)
reads_ins, freq_ins = calculateStats(otutable_stats, 'Status_ins', otus_ins, totalreads)
reads_allsgRNAs, freq_allsgRNAs = calculateStats(otutable_stats, 'Status_allsgRNAs', otus_indels_allsgRNAs, totalreads)
reads_undetermined, freq_undetermined = calculateStats(otutable_stats, 'Status_undetermined', otus_undetermined, totalreads)

#otutable_stats['Status_wt'] = ['wt' if otu in otus_wt else 'no' for otu in otutable_stats['#OTU ID']]
#otutable_stats['Status_indels'] = ['indel' if otu in otus_indels else 'no' for otu in otutable_stats['#OTU ID']]
#otutable_stats['Status_dels'] = ['del' if otu in otus_dels else 'no' for otu in otutable_stats['#OTU ID']]
#otutable_stats['Status_ins'] = ['in' if otu in otus_ins else 'no' for otu in otutable_stats['#OTU ID']]
for k, v in otus_sgRNAs.items():
   if v == 'NA': continue
   else:
      otutable_stats['Status_' + k] = ['indel' if otu in v else 'no' for otu in otutable_stats['#OTU ID']]
otutable_stats['Status_allsgRNAs'] = ['indel' if otu in otus_indels_allsgRNAs else 'no' for otu in otutable_stats['#OTU ID']]
otutable_stats['Status_undetermined'] = ['undetermined' if otu in otus_undetermined else 'no' for otu in otutable_stats['#OTU ID']]
#reads_wt = otutable_stats.loc[otutable_stats['Status_wt'] =='wt', sample].sum()
#reads_indels = otutable_stats.loc[otutable_stats['Status_indels'] =='indel', sample].sum()
#reads_dels = otutable_stats.loc[otutable_stats['Status_dels'] =='del', sample].sum()
#reads_ins = otutable_stats.loc[otutable_stats['Status_ins'] =='in', sample].sum()
reads_sgRNAs = {}
for k, v in otus_sgRNAs.items():
   if v == 'NA': reads_sgRNAs[k] = 'NA'
   else:
      reads_sgRNA = otutable_stats.loc[otutable_stats['Status_' + k] =='indel', sample].sum()
      reads_sgRNAs[k] = reads_sgRNA
#reads_allsgRNAs = otutable_stats.loc[otutable_stats['Status_allsgRNAs'] =='indel', sample].sum()
#reads_undetermined = otutable_stats.loc[otutable_stats['Status_undetermined'] =='undetermined', sample].sum()
# Percentage of reads with targets, it depends of the alpha option
#freq_wt = float((reads_wt/totalreads))*100
#freq_indels = float((reads_indels/totalreads))*100
#freq_dels = float((reads_dels/totalreads))*100
#freq_ins = float((reads_ins/totalreads))*100
freq_sgRNAs = {}
for k, v in otus_sgRNAs.items():
   if v == 'NA':
      freq_sgRNAs[k] = 'NA'
   else:
      freq_sgRNA = float((reads_sgRNAs[k]/totalreads))*100
      freq_sgRNAs[k] = freq_sgRNA
#freq_allsgRNAs = float((reads_allsgRNAs/totalreads))*100
#freq_undetermined = float((reads_undetermined/totalreads))*100

# Generate stats of indels
if dsOligo_status == 'True':
   otutable_stats['Status_dsoligo'] = ['dsoligo' if otu in otus_dsoligo else 'no' for otu in otutable_stats['#OTU ID']]
   reads_dsoligo = otutable_stats.loc[otutable_stats['Status_dsoligo'] =='dsoligo', sample].sum()
   freq_dsoligo = float((reads_dsoligo/totalreads))*100
   stats = pd.DataFrame({'Params': ['wt_otus', 'indels_otus', 'dels_otus', 'ins_otus', 'un_otus'] + ['sgRNA1_otus', 'sgRNA2_otus', 'sgRNA3_otus'] + ['allsgRNAs_otus'] +
                         ['Reads_wt', 'Reads_indels', 'Reads_dels', 'Reads_ins', 'Reads_un'] + ['Reads_sgRNA1', 'Reads_sgRNA2', 'Reads_sgRNA3'] + ['Reads_allsgRNAs'] +
                         ['Freq_wt', 'Freq_indels', 'Freq_dels', 'Freq_ins', 'Freq_un'] + ['Freq_sgRNA1', 'Freq_sgRNA2', 'Freq_sgRNA3'] + ['Freq_allsgRNAs'],
                         sample: [len(set(otus_wt)), len(set(otus_indels)), len(set(otus_dels)), len(set(otus_ins)), len(set(otus_undetermined))] + [len(set(item)) if not item == 'NA' else 'NA' for item in list(otus_sgRNAs.values())] +[len(otus_indels_allsgRNAs), len(otus_dsoligo)] +
                         [reads_wt, reads_indels, reads_dels, reads_ins, reads_undetermined] + [item if not item == 'NA' else 'NA' for item in list(reads_sgRNAs.values())] + [reads_allsgRNAs, reads_dsoligo] +
                         [freq_wt, freq_indels, freq_dels, freq_ins, freq_undetermined] + [item if not item == 'NA' else 'NA' for item in list(freq_sgRNAs.values())] + [freq_allsgRNAs, freq_dsoligo]})
else:
   stats = pd.DataFrame({'Params': ['wt_otus', 'indels_otus', 'dels_otus', 'ins_otus', 'un_otus'] + ['sgRNA1_otus', 'sgRNA2_otus', 'sgRNA3_otus'] + ['allsgRNAs_otus'] +
                         ['Reads_wt', 'Reads_indels', 'Reads_dels', 'Reads_ins', 'Reads_un'] + ['Reads_sgRNA1', 'Reads_sgRNA2', 'Reads_sgRNA3'] + ['Reads_allsgRNAs'] +
                         ['Freq_wt', 'Freq_indels', 'Freq_dels', 'Freq_ins', 'Freq_un'] + ['Freq_sgRNA1', 'Freq_sgRNA2', 'Freq_sgRNA3'] + ['Freq_allsgRNAs'],
                         sample: [len(set(otus_wt)), len(set(otus_indels)), len(set(otus_dels)), len(set(otus_ins)), len(set(otus_undetermined))] + [len(set(item)) if not item == 'NA' else 'NA' for item in list(otus_sgRNAs.values())] + [len(otus_indels_allsgRNAs)] +
                         [reads_wt, reads_indels, reads_dels, reads_ins, reads_undetermined] + [item if not item == 'NA' else 'NA' for item in list(reads_sgRNAs.values())] + [reads_allsgRNAs] +
                         [freq_wt, freq_indels, freq_dels, freq_ins, freq_undetermined] + [item if not item == 'NA' else 'NA' for item in list(freq_sgRNAs.values())] + [freq_allsgRNAs]})
stats.to_csv(args.stats, sep="\t", index=False)


#------------------------------------------------------------------------------------------------------
# Generate fasta files
#------------------------------------------------------------------------------------------------------
print('Generating fasta files...')
# Generate fasta with undetermined otus
denoisedFasta = SeqIO.to_dict(SeqIO.parse(args.denoised, 'fasta')) # otus in denoised
denoisedFasta_undetermined = open(os.path.dirname(args.sam_bwa) + '/' + sample + 'undetermined.fa', 'w+')
for k, v in denoisedFasta.items():
   otu = k.split(";")[0]
   if otu in otus_undetermined:
      denoisedFasta_undetermined.write('>' + otu + '\n')
      denoisedFasta_undetermined.write(str(v.seq) + '\n')

# Generate fasta with indels otus
denoisedFasta = SeqIO.to_dict(SeqIO.parse(args.denoised, 'fasta')) # otus in denoised
denoisedFasta_indels = open(os.path.dirname(args.sam_bwa) + '/' + sample + 'indels.fa', 'w+')
for k, v in denoisedFasta.items():
   otu = k.split(";")[0]
   if otu in otus_indels:
      denoisedFasta_indels.write('>' + otu + '\n')
      denoisedFasta_indels.write(str(v.seq) + '\n')


#------------------------------------------------------------------------------------------------------
# Stats of % indels, % of extraction and % of dsoligo if apply per Amp type
#------------------------------------------------------------------------------------------------------
print('Generating stats per Amp type...')
otus_amp_types = defaultdict(list)
for k, v in amp_types.items(): # generate a dict with otus per Amps types
   otus_in_amp = []
   for amp_name in v:
      otus_in_amp += list(set(bwotus_per_ref[amp_name] + bbotus_per_ref[amp_name]))
   otus_amp_types[k] = list(set(otus_in_amp))
stats_per_amp = defaultdict(list)
for k, v in otus_amp_types.items():
   # Number of otus
   wt_in_amp = [otu for otu in otus_wt if otu in v]
   dels_in_amp = [otu for otu in otus_dels if otu in v and not otu in wt_in_amp] # Take into account: if a otu assigned to this Amp type is a del here, but a wt in other Amp type, thi otu will be considered as wt for the current Amp type, ignoring its status for the other Amps types
   # In the case of extraction is different from dels, the extractions are found mainly in BBmap alignments, not BWA-mem, but otus aligned with BWA are assigned also to Amp types, so there are otus with extractions in Alpha7 with BBmap that are also aligned with BWA to Alpha0 amps... Corrected in extractionFragment()
   extraction_in_amp = [otu for otu in otus_extraction[k] if otu in v and not otu in wt_in_amp] # if otu in v is redundant...
   stats_per_amp['otus_wt'].append(len(wt_in_amp))
   stats_per_amp['otus_del'].append(len(dels_in_amp))
   stats_per_amp['otus_ext'].append(len(extraction_in_amp))

   # Reads
   # Filter by otus mapped
   amps_stats = otutable.copy() # only contain otus in sample (and denoised if apply (alpha option))
   amps_stats['Status_mapped'] = ['mapped' if otu in list(set(otus_indels + otus_wt + otus_undetermined)) else 'no' for otu in amps_stats['#OTU ID']]
   amps_stats = amps_stats[amps_stats['Status_mapped'] == 'mapped'] # filter otutable by otus mapped with bwa-mem or BBmap
   # Filter by otus present in the Amp type
   amps_stats['Otus_in_amp'] = ['amp' if otu in v else 'no' for otu in amps_stats['#OTU ID']]
   amps_stats = amps_stats[amps_stats['Otus_in_amp'] == 'amp']
   totalreads_in_amp = amps_stats[sample].sum() # sum of reads of all the otus in sample mapped
   # Reads of wt and dels
   amps_stats['Status_wt'] = ['wt' if otu in wt_in_amp else 'no' for otu in amps_stats['#OTU ID']]
   amps_stats['Status_dels'] = ['del' if otu in dels_in_amp else 'no' for otu in amps_stats['#OTU ID']]
   amps_stats['Status_ext'] = ['ext' if otu in extraction_in_amp else 'no' for otu in amps_stats['#OTU ID']]
   reads_wt_in_amp = amps_stats.loc[amps_stats['Status_wt'] =='wt', sample].sum()
   reads_dels_in_amp = amps_stats.loc[amps_stats['Status_dels'] =='del', sample].sum()
   reads_ext_in_amp = amps_stats.loc[amps_stats['Status_ext'] =='ext', sample].sum()
   stats_per_amp['reads_wt'].append(reads_wt_in_amp)
   stats_per_amp['reads_dels'].append(reads_dels_in_amp)
   stats_per_amp['reads_ext'].append(reads_ext_in_amp)

   # Frequencies
   freq_wt_in_amp = float((reads_wt_in_amp/totalreads_in_amp))*100
   freq_dels_in_amp = float((reads_dels_in_amp/totalreads_in_amp))*100
   freq_ext_in_amp = float((reads_ext_in_amp/totalreads_in_amp))*100
   stats_per_amp['freq_wt'].append(freq_wt_in_amp)
   stats_per_amp['freq_dels'].append(freq_dels_in_amp)
   stats_per_amp['freq_ext'].append(freq_ext_in_amp)

   # Stats for dsoligo
   if dsOligo_status == 'True': 
      dsoligo_in_amp = [otu for otu in otus_dsoligo_amp[k] if otu in v] # if otu in v is redundant... a dsoligo can be in wt otu
      stats_per_amp['otus_dsoligo'].append(len(dsoligo_in_amp))
      amps_stats['Status_dsoligo'] = ['dsoligo' if otu in dsoligo_in_amp else 'no' for otu in amps_stats['#OTU ID']]
      reads_dsoligo_in_amp = amps_stats.loc[amps_stats['Status_dsoligo'] =='dsoligo', sample].sum()
      stats_per_amp['reads_dsoligo'].append(reads_dsoligo_in_amp)
      freq_dsoligo_in_amp = float((reads_dsoligo_in_amp/totalreads_in_amp))*100
      stats_per_amp['freq_dsoligo'].append(freq_dsoligo_in_amp)
      

   # Generate csv per sample
   if dsOligo_status == 'False':
      stats_amp = pd.DataFrame({'Params': ['wt_otus', 'dels_otus', 'ext_otus', 'Reads_wt', 'Reads_dels', 'Reads_ext', 'Freq_wt', 'Freq_dels', 'Freq_ext'],
                                sample: [len(set(wt_in_amp)), len(set(dels_in_amp)), len(set(extraction_in_amp)),
                                         reads_wt_in_amp, reads_dels_in_amp, reads_ext_in_amp,
                                         freq_wt_in_amp, freq_dels_in_amp, freq_ext_in_amp]})
      stats_amp.to_csv(os.path.dirname(args.sam_bwa) + '/temp_stats_indels/' + sample + '_' + k + '_typestats.txt', sep="\t", index=False)
   else:
      stats_amp = pd.DataFrame({'Params': ['wt_otus', 'dels_otus', 'ext_otus', 'dsoligo_otus', 'Reads_wt', 'Reads_dels', 'Reads_ext', 'Reads_dsoligo', 'Freq_wt', 'Freq_dels', 'Freq_ext', 'Freq_dsoligo'],
                                sample: [len(set(wt_in_amp)), len(set(dels_in_amp)), len(set(extraction_in_amp)), len(set(dsoligo_in_amp)),
                                         reads_wt_in_amp, reads_dels_in_amp, reads_ext_in_amp, reads_dsoligo_in_amp,
                                         freq_wt_in_amp, freq_dels_in_amp, freq_ext_in_amp, freq_dsoligo_in_amp]})
      stats_amp.to_csv(os.path.dirname(args.sam_bwa) + '/temp_stats_indels/' + sample + '_' + k + '_typestats.txt', sep="\t", index=False)



print(sample + ': DONE.')
