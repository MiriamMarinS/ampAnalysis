from Bio import SeqIO
import argparse
from progressbar import *
import pandas as pd
from analysisSeq import analysisSamv2, epitopeFind, sgRNA, missingEpitopes, totalEpitopes, length_indel_events, PAMefficiency
from utils import calculateStats, generateFasta, nonRedundant, nonRedundantItem0, countIntersections
from workingwithdfs import creatingDFfromOtus
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
refsFasta = '/home/tonin/mimi/database/BW208_Amp/BW208_Amp.fasta' # ORIGINAL
#refsFasta = '/home/tonin/mimi/database/Back_T590+BW208/Back_T590+BW208.fasta' # DH
#refsFasta = '/home/tonin/mimi/database/Back_AN+BW208/Back_AN+BW208.fasta' # DH

#------------------------------------------------------------------------------------------------------
# Search the epitopes in the Amp database
#------------------------------------------------------------------------------------------------------
epitopesFile = '/home/tonin/mimi/database/epitopes/epitopes.fasta'
amp_types, hits_epitopes = epitopeFind(epitopesFile, refsFasta)


#------------------------------------------------------------------------------------------------------
# Search Indels oand/or dsoligo
#------------------------------------------------------------------------------------------------------
# Read SAM file from bwa-mem and BBmap, filter by otus in sample
otutable = pd.read_csv(args.otu_table, sep="\t", header=0)
otutable = otutable[otutable[sample] > 0] # Filter otus in sample
otus_sample = otutable['#OTU ID'].to_list() # otus in sample

# Search indels/dsoligo
bwotus_wt, bwotus_indels, bwotus_extraction, bwotus_all, bwundetermined_otus, bwotus_sgRNAs, bwotus_per_ref, bwotus_dsoligo, bwlength_indels, bwpam_efficiency = analysisSamv2(args.sam_bwa, otus_sample, args.metafile, refsFasta, sample, hits_epitopes, dsOligo_status, int(args.num_diff))
bbotus_wt, bbotus_indels, bbotus_extraction, bbotus_all, bbundetermined_otus, bbotus_sgRNAs, bbotus_per_ref, bbotus_dsoligo, bblength_indels, bbpam_efficiency = analysisSamv2(args.sam_bbmap, otus_sample, args.metafile, refsFasta, sample, hits_epitopes, dsOligo_status, int(args.num_diff))


#------------------------------------------------------------------------------------------------------
# Process otus per category
#------------------------------------------------------------------------------------------------------

# Combined indels/dsoligo
bwotus_indels_values = nonRedundant(reduce(lambda i,j:i+j,list(bwotus_indels.values()),[]))
bbotus_indels_values = nonRedundant(reduce(lambda i,j:i+j,list(bbotus_indels.values()),[]))
otus_indels = nonRedundant(bwotus_indels_values + bbotus_indels_values) # combined list (non-redundant) of otus with indels in BWA-mem or BBmap
otus_dels = nonRedundant(bwotus_indels['Dels'] + bbotus_indels['Dels'])
otus_ins = nonRedundant(bwotus_indels['Ins'] + bbotus_indels['Ins'])
otus_wt = [otu for otu in nonRedundant(bwotus_wt + bbotus_wt) if otu not in otus_indels] # combined list (non-redundant) of otus without indels in BW-mem or BBmap
otus_undeterminedtemp = nonRedundant(bwundetermined_otus + bbundetermined_otus)
otus_undetermined = [otu for otu in otus_undeterminedtemp if otu not in otus_indels and otu not in otus_wt]
bwotus_dsoligo_values = reduce(lambda i,j:i+j,list(bwotus_dsoligo.values()),[])
bbotus_dsoligo_values = reduce(lambda i,j:i+j,list(bbotus_dsoligo.values()),[])
otus_dsoligo = nonRedundant([item[0] for item in bwotus_dsoligo_values] + [item[0] for item in bbotus_dsoligo_values]) # combined list (non-redundant) of otus with dsoligo in BWA-mem or BBmap
otus_dsoligo_fwd = nonRedundant([item[0] for item in bwotus_dsoligo_values if item[1] == 'fwd'] + [item[0] for item in bbotus_dsoligo_values if item[1] == 'fwd'])
otus_dsoligo_rev = nonRedundant([item[0] for item in bwotus_dsoligo_values if item[1] == 'rev'] + [item[0] for item in bbotus_dsoligo_values if item[1] == 'rev'])

# Include otus with perfect substitution of part of Amp with the dsoligo in InDels statistics
if dsOligo_status == 'True':
   for otu in otus_dsoligo:
      # if otu with dsoligo not classified as InDel, then remove from the other lists and include in indels list as putative deletion followed with an insertion of dsoligo
      if any(otu in list_otus for list_otus in [otus_wt, otus_undetermined]):
         if otu in otus_wt: otus_wt.remove(otu)
         if otu in otus_undetermined: otus_undetermined.remove(otu)
         # include in otus indels
         if otu not in otus_indels:
            otus_indels.append(otu)
            otus_dels.append(otu)
            otus_ins.append(otu)

# Combined indels per sgRNA
otus_sgRNAs = {}
otus_per_guide = {}
for k in nonRedundant(list(bwotus_sgRNAs.keys()) + list(bbotus_sgRNAs.keys())):
   # otus with perfect substitution, add to otus_sgRNAs
   otus_per_guide[k] = nonRedundant([item[0] for item in bwotus_dsoligo_values if k in item[2]] + [item[0] for item in bbotus_dsoligo_values if k in item[2]])
   # otus with indels per sgRNAs
   otus_sgRNAs[k] = nonRedundant(bwotus_sgRNAs[k] + bbotus_sgRNAs[k])
# Combine otus from dsoligo and perfect substitution with otus per sgRNAs
if dsOligo_status == 'True':
   for k, v in otus_per_guide.items():
      for otu in v:
         if otu not in list(otus_sgRNAs[k]): otus_sgRNAs[k].append(otu)
# if all the sgRNAs of otu has indels
otus_indels_allsgRNAs = list(set.intersection(*map(set,list(otus_sgRNAs.values())))) # otus with indels in all the sgRNAs

# Combined extractions
otus_extraction = defaultdict(list)
for amptype in list(amp_types.keys()):
   for element in nonRedundant(bwotus_extraction[amptype] + bbotus_extraction[amptype]):
      otus_extraction[amptype].append(element)

# Combined dsoligo per amp type/rname
if dsOligo_status == 'True':
   otus_dsoligo_rname, otus_dsoligo_amp = defaultdict(list), defaultdict(list)
   for rname in nonRedundant(itertools.chain(*[ampName for _, ampName in amp_types.items()])):
      for element in nonRedundantItem0(bwotus_dsoligo[rname] + bbotus_dsoligo[rname]):
         otus_dsoligo_rname[rname].append(element)
         otus_dsoligo_amp[[i for i in amp_types if rname in amp_types[i]][0]].append(element)


#------------------------------------------------------------------------------------------------------
# Get otutable filtered by otus mapped to the Amp ref database
#------------------------------------------------------------------------------------------------------

# BWA-mem and BBmap: filter by sample, by alpha if apply, and by mapped
otutable_stats = otutable.copy() # only contain otus in sample (and denoised if apply (alpha option))
otutable_stats['Status_mapped'] = ['mapped' if otu in nonRedundant(otus_indels + otus_wt + otus_undetermined) else 'no' for otu in otutable_stats['#OTU ID']]
otutable_stats = otutable_stats[otutable_stats['Status_mapped'] == 'mapped'] # filter otutable by otus mapped with bwa-mem or BBmap

#------------------------------------------------------------------------------------------------------
# Get otus per Amp type and number of sgRNAs
#------------------------------------------------------------------------------------------------------
sgRNAs_number = {} # number of sgRNAs hits per Amp name
for rname in nonRedundant(itertools.chain(*[ampName for _, ampName in amp_types.items()])):
   metaFile = pd.read_csv(args.metafile, sep="\t", header=None)
   plasmid = metaFile[metaFile[0] == sample][3].values[0]
   refs = SeqIO.to_dict(SeqIO.parse(refsFasta, 'fasta'))
   list_hits, _ = sgRNA(str(refs[rname].seq), plasmid)
   sgRNAs_number[rname] = len(list_hits)


#------------------------------------------------------------------------------------------------------
# Stats of dsoligo per number of sgRNAs
#------------------------------------------------------------------------------------------------------

# The stats will be calculated as (reads otus with dsoligo and x number of sgRNAs / reads total otus with x number of sgRNAs)*100
# Per number of hits: some hits of different sgRNAs are considered different hits
if dsOligo_status == 'True':
   print('Generating stats of dsoligo per number of sgRNAs...')
   stats_dsoligo_number_sgRNAs = {}
   
   for number in nonRedundant(sgRNAs_number.values()):
      # otus aligned in ref with x number of sgRNAs
      rnames_ = [i for i in sgRNAs_number if sgRNAs_number[i] == number] # all rnames with x number of sgRNAs
      otus_NsgRNAstotal = nonRedundant(itertools.chain(*[bwotus_per_ref[rname] + bbotus_per_ref[rname] for rname in rnames_])) # otus aligned to rname with x number of sgRNAs

      # otus with dsoligo inserted in equal number of 
      otus_dsoligoNsgRNAs_temp = nonRedundantItem0(itertools.chain(*[v for k, v in otus_dsoligo_rname.items() if k in rnames_])) # otus aligned to rname with x number of sgRNAs and dsoligo
      otus_dsoligoNsgRNAs = [item[0] for item in otus_dsoligoNsgRNAs_temp]

      otutable_stats_sgRNAs = otutable_stats.copy() # only contain otus in sample (and denoised if apply (alpha option))

      # Calculate total reads
      otutable_stats_sgRNAs['Status_NsgRNA'] = ['NsgRNA' if otu in otus_NsgRNAstotal else 'no' for otu in otutable_stats_sgRNAs['#OTU ID']]
      otutable_rname = otutable_stats_sgRNAs[otutable_stats_sgRNAs['Status_NsgRNA'] == 'NsgRNA']
      totalreads = otutable_rname[sample].sum() # sum of reads of all otus by rname in the sample

      # Calculate reads and frequencies of otus with dsOligo (aligned to ref Amps with x number od sgRNAs in each iterations)
      reads_NsgRNA, freq_NsgRNA = calculateStats(otutable_stats_sgRNAs, 'Status_dsoligo_NsgRNA', otus_dsoligoNsgRNAs, totalreads, sample)

      stats_dsoligo_number_sgRNAs['otus_ODN_NsgRNAs' + str(number)] = len(otus_dsoligoNsgRNAs)
      stats_dsoligo_number_sgRNAs['reads_ODN_NsgRNAs' + str(number)] = reads_NsgRNA
      stats_dsoligo_number_sgRNAs['freq_ODN_NsgRNAs' + str(number)] = freq_NsgRNA

   stats_dsoligo_pernumber = pd.DataFrame({'Params': ['otus_ODN_NsgRNAs' + str(item) for item in nonRedundant(sgRNAs_number.values())] +
                                          ['reads_ODN_NsgRNAs' + str(item) for item in nonRedundant(sgRNAs_number.values())] +
                                          ['freq_ODN_NsgRNAs' + str(item) for item in nonRedundant(sgRNAs_number.values())], 
                                          sample: [stats_dsoligo_number_sgRNAs['otus_ODN_NsgRNAs' + str(item)] for item in nonRedundant(sgRNAs_number.values())] +
                                                   [stats_dsoligo_number_sgRNAs['reads_ODN_NsgRNAs' + str(item)] for item in nonRedundant(sgRNAs_number.values())] +
                                                   [stats_dsoligo_number_sgRNAs['freq_ODN_NsgRNAs' + str(item)] for item in nonRedundant(sgRNAs_number.values())]})
   stats_dsoligo_pernumber.to_csv(os.path.dirname(args.stats) + '/' + sample + '_statsnumbersgRNAs.txt', sep="\t", index=False)


#------------------------------------------------------------------------------------------------------
# Stats of dsoligo per number of different sgRNAs covered (if it is inserted in 1 sgRNA, or between 2 sgRNAs....)
#------------------------------------------------------------------------------------------------------

if dsOligo_status == 'True':
   print('Generating stats of sgRNAs covered by the dsoligo insertion...')
   stats_dsoligo_sgRNAs_cover = {}
   otusPerNumber = countIntersections(otus_per_guide)
   
   otutable_stats_N = otutable_stats.copy()
   # total reads
   totalotus = nonRedundant(itertools.chain(*[otus for otus in list(otus_per_guide.values())])) # otus with insertion of dsoligo in any of their sgRNAs
   otutable_stats_N['Status_NsgRNAwithIn'] = ['yes' if otu in totalotus else 'no' for otu in otutable_stats_N['#OTU ID']]
   otutable_stats_N = otutable_stats_N[otutable_stats_N['Status_NsgRNAwithIn'] == 'yes']
   totalreads = otutable_stats_N[sample].sum() # sum of reads of all otus by dsoligo In in the sample
   
   # Calculate stats for ODN insertion in x number of sgRNAs
   for k, v in otusPerNumber.items():
      stats_dsoligo_sgRNAs_cover['otus_ODN_in' + str(k)] = len(v)
      
      reads_dsIn, freq_dsIn = calculateStats(otutable_stats_N, 'Status_dsoligo_In', v, totalreads, sample)

      stats_dsoligo_sgRNAs_cover['reads_ODN_in_' + str(k)] = reads_dsIn
      stats_dsoligo_sgRNAs_cover['freq_ODN_in_' + str(k)] = freq_dsIn
   
   stats_dsoligo_perN = pd.DataFrame({'Params': ['otus_ODN_in' + str(item) for item in list(otusPerNumber.keys())] +
                                          ['reads_ODN_in' + str(item) for item in list(otusPerNumber.keys())] +
                                          ['freq_ODN_in' + str(item) for item in list(otusPerNumber.keys())], 
                                          sample: [stats_dsoligo_sgRNAs_cover['otus_ODN_in' + str(item)] for item in list(otusPerNumber.keys())] +
                                                   [stats_dsoligo_sgRNAs_cover['reads_ODN_in_' + str(item)] for item in list(otusPerNumber.keys())] +
                                                   [stats_dsoligo_sgRNAs_cover['freq_ODN_in_' + str(item)] for item in list(otusPerNumber.keys())]})
   stats_dsoligo_perN.to_csv(os.path.dirname(args.stats) + '/' + sample + '_statsODNinNumbersgRNAs.txt', sep="\t", index=False)

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
# Calculate frequency of missing epitopes in the sample
#------------------------------------------------------------------------------------------------------
print('Generating stats of missing epitopes...')
denoised = {k.split(";")[0].strip(): v for k, v in SeqIO.to_dict(SeqIO.parse(args.denoised, 'fasta')).items()} # denoised
missingEp_dict = missingEpitopes(otus_indels, bwotus_per_ref, bbotus_per_ref, hits_epitopes, refs, amp_types, sample, otus_sample, otutable_stats, epitopesFile, denoised) # Dict with % of missing epitopes
missinfEp_df = pd.DataFrame({'Epitopes': list(missingEp_dict.keys()),
                              sample: list(missingEp_dict.values())})
missinfEp_df.to_csv(os.path.dirname(args.stats) + '/' + sample + '_missingEpitopes.txt', sep="\t", index=False)


#------------------------------------------------------------------------------------------------------
# Generate principal stats: % of indels, % of indels per sgRNA, % of dsoligo if apply
#------------------------------------------------------------------------------------------------------
print('Generating stats of indels/dsoligo...')

totalreads = otutable_stats[sample].sum() # sum of reads of all the otus in sample mapped

reads_wt, freq_wt = calculateStats(otutable_stats, 'Status_wt', otus_wt, totalreads, sample)
reads_indels, freq_indels = calculateStats(otutable_stats, 'Status_indels', otus_indels, totalreads, sample)
reads_dels, freq_dels = calculateStats(otutable_stats, 'Status_dels', otus_dels, totalreads, sample)
reads_ins, freq_ins = calculateStats(otutable_stats, 'Status_ins', otus_ins, totalreads, sample)
reads_undetermined, freq_undetermined = calculateStats(otutable_stats, 'Status_undetermined', otus_undetermined, totalreads, sample)

# Get otus per each sgRNA (based on ref Amps with each sgRNA hit) and otus with all the different sgRNAs of the plasmid
allotus_per_sgRNA = defaultdict(list)
sgRNAs = {'pSSL9': ['GTTGTGATGGAAATGGTTG', 'GTTGTGGTCGAAATGGTTG'],
              'pSSL10': ['GTTGTGATGGAAATGGTTG', 'GTCGAAATGGTTGCGGCTG'],
              'Cpf1': ['CATCACAACAACCATATCTG', 'CGCAGCCGCAACTACCATAT', 'GACCACAACAACCATATCCA'],
              'pAlpha2': ['GTTGTGATGGAAATGGTTG'],
              'pAlpha1+Alpha2': ['CCACAAGAGCAAGTTCCAT', 'GTTGTGATGGAAATGGTTG'],
              'pAlpha1': ['CCACAAGAGCAAGTTCCAT'],
              'primeDelAlpha': ['CCACAAGAGCAAGTTCCAT', 'GTTGTGGTCGAAATGGTTG'],
              'primeEditingAlpha': ['GGCTGCGGATATGGTAGTTG'],
              'pFB046': ['TACAACAACAACAATTTCTA', 'GTTGTGGTCGAAATGGTTG'],
              'pFB047': ['TACAACAACAACAATTTCTA', 'AATGGTTGTGGCTGCGAATA'],
              'pFB048+50': ['TACAACAACAACAATTTCTA','GTTGTGGTCGAAATGGTTG', 'AATGGTTGTGGCTGCGAATA']}
for rname, v in refs.items():
    hits_list, hits = sgRNA(str(v.seq), plasmid) # (z, w) positions of sgRNAs in ref sequence
    all_sgRNAs = sgRNAs[plasmid]
    for item in list(set(hits.keys())): # for each sgRNA found in Amp
        otus_r = nonRedundant(bwotus_per_ref[rname] + bbotus_per_ref[rname])
        for otu_i in otus_r: allotus_per_sgRNA[item].append(otu_i)
    if len(list(set(hits.keys()))) == len(all_sgRNAs):
        otus_r = nonRedundant(bwotus_per_ref[rname] + bbotus_per_ref[rname])
        for otu_i in otus_r: allotus_per_sgRNA['all_sgRNAs'].append(otu_i)

# Stats of all different sgRNAs efficiency
# Instead of totalreads, the reads of otus mapped to Amps with all the different sgRNAs (sgAlpha2 + sgAlpha9, etc) hits in they
otutable_stats['Status_' + 'all_sgRNAs'] = ['sgRNA' if otu in nonRedundant(allotus_per_sgRNA['all_sgRNAs']) else 'no' for otu in otutable_stats['#OTU ID']]
otutable_stats_allsgRNA = otutable_stats[otutable_stats['Status_' + 'all_sgRNAs'] == 'sgRNA'] # filter otutable by otus with all the different sgRNAs
totalreads_allsgRNAs = otutable_stats_allsgRNA[sample].sum() # sum of reads of all the otus in sample mapped to Amps with all the different sgRNAs of the plasmid
# to calculate the frequency of indels in all the different sgRNAs of each construct (sgAlpha2 + sgAlpha9, etc) based on all the reads with these sgRNAs hits
reads_allsgRNAs, freq_allsgRNAs = calculateStats(otutable_stats, 'Status_allsgRNAs', otus_indels_allsgRNAs, totalreads_allsgRNAs, sample) # frequency of reads with InDels in all the different types of sgRNAs of the plasmid (not all the hits). If Cpf1 has 3 different sgRNAs, only reads with InDels in those three will be counted

# Stats of sgRNA efficiency
reads_sgRNAs, freq_sgRNAs = {}, {}
for k, v in otus_sgRNAs.items():
   otutable_stats['Status_' + k] = ['sgRNA' if otu in nonRedundant(allotus_per_sgRNA[k]) else 'no' for otu in otutable_stats['#OTU ID']]
   otutable_stats_sgRNA = otutable_stats[otutable_stats['Status_' + k] == 'sgRNA'] # filter otutable by otus with sgRNA
   totalreads_sgRNA = otutable_stats_sgRNA[sample].sum() # sum of reads of all the otus in sample with this sgRNA
   # to calculate the frequency of indels in sgRNA based on all the reads with this sgRNA hit
   reads_sgRNA, freq_sgRNA = calculateStats(otutable_stats, 'Status_' + k, v, totalreads_sgRNA, sample) # Instead of totalreads, indicate totalreads with each sgRNA
   reads_sgRNAs[k], freq_sgRNAs[k] = reads_sgRNA, freq_sgRNA

# Generate stats of indels
if dsOligo_status == 'True':
   stats_dsoligo = {}
   for otus_list in [('Status_dsoligo', otus_dsoligo), ('Status_dsoligofwd', otus_dsoligo_fwd), ('Status_dsoligorev', otus_dsoligo_rev)]:
      typeEvent = otus_list[0].split('_')[1]
      reads_dsoligo, freq_dsoligo = calculateStats(otutable_stats, otus_list[0], otus_list[1], totalreads, sample)
      stats_dsoligo['otus_' + typeEvent] = len(otus_list[1])
      stats_dsoligo['Reads_' + typeEvent] = reads_dsoligo
      stats_dsoligo['Freq_' + typeEvent] = freq_dsoligo
      stats_dsoligo['FreqPerIns_' + typeEvent] = float((reads_dsoligo/reads_ins)*100)


   stats = pd.DataFrame({'Params': ['wt_otus', 'indels_otus', 'dels_otus', 'ins_otus', 'un_otus'] + [item + '_otus' for item in list(otus_sgRNAs.keys())] + ['allsgRNAs_otus', 'otus_dsoligo', 'otus_dsoligofwd', 'otus_dsoligorev'] +
                         ['Reads_wt', 'Reads_indels', 'Reads_dels', 'Reads_ins', 'Reads_un'] + ['Reads_' + item for item in list(reads_sgRNAs.keys())] + ['Reads_allsgRNAs', 'Reads_dsoligo', 'Reads_dsoligofwd', 'Reads_dsoligorev'] +
                         ['Freq_wt', 'Freq_indels', 'Freq_dels', 'Freq_ins', 'Freq_un'] + ['Freq_' + item for item in list(freq_sgRNAs.keys())] + ['Freq_allsgRNAs', 'Freq_dsoligo', 'Freq_dsoligofwd', 'Freq_dsoligorev'] +
                         ['FreqPerIns_dsoligo', 'FreqPerIns_dsoligofwd', 'FreqPerIns_dsoligorev'],
                         sample: [len(set(otus_wt)), len(set(otus_indels)), len(set(otus_dels)), len(set(otus_ins)), len(set(otus_undetermined))] + [len(set(item)) for item in list(otus_sgRNAs.values())] + [len(otus_indels_allsgRNAs), stats_dsoligo['otus_dsoligo'], stats_dsoligo['otus_dsoligofwd'], stats_dsoligo['otus_dsoligorev']] +
                         [reads_wt, reads_indels, reads_dels, reads_ins, reads_undetermined] + [item for item in list(reads_sgRNAs.values())] + [reads_allsgRNAs, stats_dsoligo['Reads_dsoligo'], stats_dsoligo['Reads_dsoligofwd'], stats_dsoligo['Reads_dsoligorev']] +
                         [freq_wt, freq_indels, freq_dels, freq_ins, freq_undetermined] + [item for item in list(freq_sgRNAs.values())] + [freq_allsgRNAs, stats_dsoligo['Freq_dsoligo'], stats_dsoligo['Freq_dsoligofwd'], stats_dsoligo['Freq_dsoligorev']] +
                         [stats_dsoligo['FreqPerIns_dsoligo'], stats_dsoligo['FreqPerIns_dsoligofwd'], stats_dsoligo['FreqPerIns_dsoligorev']]})
else:
   stats = pd.DataFrame({'Params': ['wt_otus', 'indels_otus', 'dels_otus', 'ins_otus', 'un_otus'] + [item + '_otus' for item in list(otus_sgRNAs.keys())] + ['allsgRNAs_otus'] +
                         ['Reads_wt', 'Reads_indels', 'Reads_dels', 'Reads_ins', 'Reads_un'] + ['Reads_' + item for item in list(reads_sgRNAs.keys())] + ['Reads_allsgRNAs'] +
                         ['Freq_wt', 'Freq_indels', 'Freq_dels', 'Freq_ins', 'Freq_un'] + ['Freq_' + item for item in list(freq_sgRNAs.keys())] + ['Freq_allsgRNAs'],
                         sample: [len(set(otus_wt)), len(set(otus_indels)), len(set(otus_dels)), len(set(otus_ins)), len(set(otus_undetermined))] + [len(set(item)) for item in list(otus_sgRNAs.values())] + [len(otus_indels_allsgRNAs)] +
                         [reads_wt, reads_indels, reads_dels, reads_ins, reads_undetermined] + [item for item in list(reads_sgRNAs.values())] + [reads_allsgRNAs] +
                         [freq_wt, freq_indels, freq_dels, freq_ins, freq_undetermined] + [item for item in list(freq_sgRNAs.values())] + [freq_allsgRNAs]})
stats.to_csv(args.stats, sep="\t", index=False)


#------------------------------------------------------------------------------------------------------
# Calculate frequency of epitopes in the sample
#------------------------------------------------------------------------------------------------------
print('Generating stats of epitopes...')
ep_dict = totalEpitopes(bwotus_per_ref, bbotus_per_ref, refs, amp_types, sample, otus_sample, otutable_stats, epitopesFile, denoised)
ep_df = pd.DataFrame({'Epitopes': list(ep_dict.keys()),
                              sample: list(ep_dict.values())})
ep_df.to_csv(os.path.dirname(args.stats) + '/' + sample + '_epitopes.txt', sep="\t", index=False)


#------------------------------------------------------------------------------------------------------
# Calculate the frequency of each indel event based on their length in the sample
#------------------------------------------------------------------------------------------------------
print('Generating stats for indel event length...')
indel_ev_len_dict = length_indel_events(bwlength_indels, bblength_indels, otus_indels, sample, otutable_stats)
ev_df = pd.DataFrame({'Events': list(indel_ev_len_dict.keys()),
                              sample: list(indel_ev_len_dict.values())})
ev_df.to_csv(os.path.dirname(args.stats) + '/' + sample + '_lengthindels.txt', sep="\t", index=False)


#------------------------------------------------------------------------------------------------------
# Calculate the efficiency of PAM per guide in the sample
#------------------------------------------------------------------------------------------------------
print('Generating stats of PAM efficiency...')
# Get otus per each PAM+sgRNA (based on ref Amps with each PAM+sgRNA hit)
allotus_perPAM_sgRNA = defaultdict(list)
# sgRNAs dict defined in "Generate principal stats" section
for rname, v in refs.items():
   hits_list, hits = sgRNA(str(v.seq), plasmid) # (z, w) positions of sgRNAs in ref sequence + PAM seq found
   all_sgRNAs = sgRNAs[plasmid]
    
   # Get all combinations PAM+sgRNA
   PAMplussgRNA = []
   for guideName, items in hits.items(): # for each sgRNA found in Amp
      for item in items: PAMplussgRNA.append(guideName + '+' + item[2]) # being item[2] = PAM seq
   # Get all otus per PAM+sgRNA combinations
   otus_r = nonRedundant(bwotus_per_ref[rname] + bbotus_per_ref[rname]) # otus assigned to this ref Amp
   for PAMsgRNA in PAMplussgRNA:
      for otu_i in otus_r: allotus_perPAM_sgRNA[PAMsgRNA].append(otu_i)


pam_efficiency_dict = PAMefficiency(bwpam_efficiency, bbpam_efficiency, otus_indels, sample, otutable_stats, allotus_perPAM_sgRNA) # Not considering sgRNAs in otus with ODN insertion and no D or I, and considering multiple events of I or D in one sgRNAs
pam_df = pd.DataFrame({'PAM+guide': list(pam_efficiency_dict.keys()),
                       sample: list(pam_efficiency_dict.values())})
pam_df.to_csv(os.path.dirname(args.stats) + '/' + sample + '_pamefficiency.txt', sep="\t", index=False)


#------------------------------------------------------------------------------------------------------
# Generate fasta files
#------------------------------------------------------------------------------------------------------
print('Generating fasta files...')
generateFasta(args.denoised, os.path.dirname(args.sam_bwa) + '/' + sample + 'undetermined.fa', otus_undetermined) # Generate fasta with undetermined otus
generateFasta(args.denoised, os.path.dirname(args.sam_bwa) + '/' + sample + 'indels.fa', otus_indels) # Generate fasta with indels otus
if dsOligo_status == 'True':
   generateFasta(args.denoised, os.path.dirname(args.sam_bwa) + '/' + sample + 'dsoligo.fa', otus_dsoligo) # Generate fasta with dsoligo otus


#------------------------------------------------------------------------------------------------------
# Stats of % indels, % of extraction and % of dsoligo if apply per Amp type
#------------------------------------------------------------------------------------------------------
print('Generating stats per Amp type...')

# Generate a dict with otus per Amp type
otus_amp_types = defaultdict(list)
for ampType, ampNames in amp_types.items():
   otus_in_amp = []
   for amp_name in ampNames:
      otus_in_amp += nonRedundant(bwotus_per_ref[amp_name] + bbotus_per_ref[amp_name])
   otus_amp_types[ampType] = nonRedundant(otus_in_amp)

# Generate stats of extraction, dels and dsoligo if apply per Amp type
for ampType, otusperType in otus_amp_types.items():
   
   # ----------------------- First part -----------------------
   # Stats per Amp type
   # ----------------------------------------------------------

   outputname = sample + '_' + ampType + '_typestats.txt'
   if dsOligo_status == 'True':
      creatingDFfromOtus(otusperType, ampType, otus_indels, otus_wt, otus_undetermined, otutable, sample, dsOligo_status, args.sam_bwa, outputname, otus_dels, otus_extraction, otus_dsoligo_amp)
   else:
      creatingDFfromOtus(otusperType, ampType, otus_indels, otus_wt, otus_undetermined, otutable, sample, dsOligo_status, args.sam_bwa, outputname, otus_dels, otus_extraction, '')
   

   # ----------------------- Second part -----------------------
   # Stats per Amp type and per number of sgRNAs
   # -----------------------------------------------------------

   # Get otu list per number of sgRNAs and per Amp type
   for number in nonRedundant(sgRNAs_number.values()):
      rnames_ = [i for i in sgRNAs_number if sgRNAs_number[i] == number] # all rnames with x number of sgRNAs
      rnamesNsgRNAs_per_amp = [item for item in rnames_ if item in amp_types[ampType]] # all rnames with x number of sgRNAs in this Amp type
      otus_NsgRNAstotal_per_amp = nonRedundant(itertools.chain(*[bwotus_per_ref[rname] + bbotus_per_ref[rname] for rname in rnamesNsgRNAs_per_amp])) # otus aligned to rname with x number of sgRNAs and from this Amp type
      
      outputname = sample + '_' + ampType + '-' + str(number) + '_typenumberstats.txt'
      if dsOligo_status == 'True':
         creatingDFfromOtus(otus_NsgRNAstotal_per_amp, ampType, otus_indels, otus_wt, otus_undetermined, otutable, sample, dsOligo_status, args.sam_bwa, outputname, otus_dels, otus_extraction, otus_dsoligo_amp)
      else:
         creatingDFfromOtus(otus_NsgRNAstotal_per_amp, ampType, otus_indels, otus_wt, otus_undetermined, otutable, sample, dsOligo_status, args.sam_bwa, outputname, otus_dels, otus_extraction, '')
   

   



print(sample + ': DONE.')
