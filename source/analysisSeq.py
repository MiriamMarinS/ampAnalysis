# Functions and classes for ampAnalysis.py
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re
from collections import defaultdict
import itertools
from utils import appendList, nonRedundant
from functools import reduce
from Bio import Align

#----------------------------------------------------------------------------------------------------------------------------------------
# Align a pattern in a text with x mismatches
#----------------------------------------------------------------------------------------------------------------------------------------

def naive_mm(p, t, n_diffs):
    occurrences = []
    for i in range(len(t) - len(p) + 1): # loop over alignments
        match = True
        error_num = 0
        for j in range(len(p)):
            if t[i+j] != p[j]:
                if error_num == n_diffs:
                    match = False
                    break
                else:
                    error_num += 1
                    continue
        if match:
            occurrences.append(i)
    return(occurrences)


#----------------------------------------------------------------------------------------------------------------------------------------
# Classify Amps from DB by their number of epitopes in Amp types
#----------------------------------------------------------------------------------------------------------------------------------------

def epitopeFind(epitope_db, alpha_db):
    # Import the epitope fasta
    epitopesFile = SeqIO.to_dict(SeqIO.parse(epitope_db, 'fasta'))
    # Import the Amp fasta
    refs = SeqIO.to_dict(SeqIO.parse(alpha_db, 'fasta'))
    
    # Search the epitopes in the peptide sequences of Amps
    eps_in_alpha = defaultdict(list)
    hits_epitopes = defaultdict(list)
    for k, v in refs.items():
        while len(v.seq)%3 != 0:
            v = v + Seq('N')
        for ep_name, ep_seq in epitopesFile.items():
            eps = naive_mm(ep_seq.seq, v.translate().seq, 0)
            eps_in_alpha[k].append(len(eps))
            hits_epitopes[k].append(eps) # hits of the epitopes per amp to use in extractionFragment()
    
    # Classify the Amps by their number of epitopes in alpha-gliadins
    alpha_types = defaultdict(list)
    for k, v in eps_in_alpha.items():
        number_epitopes = sum(v[0:4]) # alpha epitopes
        type_amp = 'Alpha' + str(number_epitopes)
        alpha_types[type_amp].append(k)
    return(alpha_types, hits_epitopes) # return a dict with the Amp types as keys and all the Amps classified as values


#----------------------------------------------------------------------------------------------------------------------------------------
# Check if the PAM alignment is perfect
#----------------------------------------------------------------------------------------------------------------------------------------

def PAMsearch(plasmid, hit, ref, guide, sense):
    PAM = {'pSSL9': ('NGG', '3\''),
           'pSSL10': ('NGG', '3\''),
           'Cpf1': ('TTTN', '5\''),
           'pAlpha2': ('NGG', '3\''),
           'pAlpha1+Alpha2': ('NGG', '3\''),
           'pAlpha1': ('NGG', '3\''),
           'primeDelAlpha': ('NGG', '3\''),
           'primeEditingAlpha': ('NGG', '3\''),
           'pFB046': ('NGG', '3\''),
           'pFB047': ('NGG', '3\''),
           'pFB048+50': ('NGG', '3\'')}
    pam = PAM[plasmid][0]
    site = PAM[plasmid][1]
    found = 'False'
    pam_seq = ''
    for nucleotide in ['A', 'C', 'T', 'G']:
        if site == '5\'':
            if sense == 'fwd':
                if ref[hit-4:hit] == pam.replace('N', nucleotide):
                    found = 'True'
                    pam_seq = ref[hit-4:hit]
                    break
            else:
                if ref[hit+len(guide):hit+len(guide)+4] == str(Seq(pam).reverse_complement()).replace('N', nucleotide):
                    found = 'True'
                    pam_seq = ref[hit+len(guide):hit+len(guide)+4]
                    break
        else:
            if sense == 'fwd':
                if ref[hit+len(guide):hit+len(guide)+3] == pam.replace('N', nucleotide):
                    found = 'True'
                    pam_seq = ref[hit+len(guide):hit+len(guide)+3]
                    break
            else:
                if ref[hit-3:hit] == str(Seq(pam).reverse_complement()).replace('N', nucleotide):
                    found = 'True'
                    pam_seq = ref[hit-3:hit]
                    break
    return((found, pam_seq))


#----------------------------------------------------------------------------------------------------------------------------------------
# Check if the seed sequence alignment is perfect: seed sequence 10 bp downstream of the PAM sequence
#----------------------------------------------------------------------------------------------------------------------------------------

def perfectSeed(plasmid, hit, p, t, sense):
    found = 'False'
    PAM = {'pSSL9': ('NGG', '3\''),
           'pSSL10': ('NGG', '3\''),
           'Cpf1': ('TTTN', '5\''),
           'pAlpha2': ('NGG', '3\''),
           'pAlpha1+Alpha2': ('NGG', '3\''),
           'pAlpha1': ('NGG', '3\''),
           'primeDelAlpha': ('NGG', '3\''),
           'primeEditingAlpha': ('NGG', '3\''),
           'pFB046': ('NGG', '3\''),
           'pFB047': ('NGG', '3\''),
           'pFB048+50': ('NGG', '3\'')}
    
    seedLength = 12
    # Check if the seed sequence has not mismatches
    if PAM[plasmid][1] == '5\'':
        if sense == 'fwd' and p[:seedLength] == t[hit:hit+seedLength]: # seed sequence
            found = 'True'
        if sense == 'rev' and p[len(p)-seedLength:len(p)] == t[hit+len(p)-seedLength:hit+len(p)]: # seed sequence
            found = 'True'
    else:
        if sense == 'fwd' and p[len(p)-seedLength:len(p)] == t[hit+len(p)-seedLength:hit+len(p)]: # seed sequence
            found = 'True'
        if sense == 'rev' and p[:seedLength] == t[hit:hit+seedLength]: # seed sequence
            found = 'True'
    return(found)


#----------------------------------------------------------------------------------------------------------------------------------------
# Search the sgRNAs hits: taking into account the perfect PAM and the perfect seed sequence
#----------------------------------------------------------------------------------------------------------------------------------------

def sgRNA(ref, plasmid):
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
    hits = defaultdict(list)
    for guide in sgRNAs[plasmid]:
        hits_tuple_list = []
        for hit in naive_mm(guide, ref, 3): # Naive algorithm is 0-based
            if PAMsearch(plasmid, hit, ref, guide, 'fwd')[0] == 'True' and perfectSeed(plasmid, hit, guide, ref, 'fwd') == 'True':
                hits_tuple = (hit+1, hit+1+len(guide), PAMsearch(plasmid, hit, ref, guide, 'fwd')[1])
                if hits_tuple not in hits_tuple_list:
                    hits[guide].append(hits_tuple) # append 1-based start position of guide and 1-based end position of guide
                    hits_tuple_list.append(hits_tuple)
        for hit in naive_mm(Seq(guide).reverse_complement(), ref, 3): # Search the reverse complement of sgRNA
            if PAMsearch(plasmid, hit, ref, guide, 'rev')[0] == 'True' and perfectSeed(plasmid, hit, Seq(guide).reverse_complement(), ref, 'rev') == 'True':
                hits_tuple = (hit+1, hit+1+len(guide), PAMsearch(plasmid, hit, ref, guide, 'rev')[1])
                if hits_tuple not in hits_tuple_list:
                    hits[guide].append(hits_tuple)
                    hits_tuple_list.append(hits_tuple)
    # List with hits with no sgRNAs identification
    hits_list = []
    for item in list(hits.values()):
        for element in item:
            hits_list.append(element)
    hits_list = list(set(hits_list))

    return(hits_list, hits) # hits (z, w) of sgRNA in ref in hits_list, hits is a dict with unique hits per sgRNA to calculate the efficiency per sgRNA

#----------------------------------------------------------------------------------------------------------------------------------------
# Checking if the deletion is a CRISPR InDel
#----------------------------------------------------------------------------------------------------------------------------------------

def workingwithDels(Pi, Piq, part, m, z, w, hit, hits, query, ref, indels, sgRNAefficiency, hits_del, pam_efficiency):
    PiD = Pi # Intial position of Del
    PfD = PiD + m # Final position of Del
    pam_seq = hit[2]
    if PiD < w and PfD > z:
        appendList([indels['Dels'], sgRNAefficiency], part, [k for k,v in hits.items() if hit in v])
        appendList([pam_efficiency], pam_seq, [k for k,v in hits.items() if hit in v])
        hits_del.append(hit)
    else: # check if there is an InDel moving the alignment
        num = 0
        while True:
            if PiD < z: # deletion in the left of sgRNA
                try:
                    if query[Piq+num-1] == ref[PiD+num-1]:
                        num += 1
                        if PiD+num < w and PfD+num > z:
                            appendList([indels['Dels'], sgRNAefficiency], part, [k for k,v in hits.items() if hit in v])
                            appendList([pam_efficiency], pam_seq, [k for k,v in hits.items() if hit in v])
                            hits_del.append(hit)
                            break
                    else:
                        break
                except:
                    break
            else: # deletion in the right of sgRNA
                try:
                    if query[Piq-num-2] == ref[PfD-num-2]:
                        num += 1
                        if PiD-num < w and PfD-num > z:
                            appendList([indels['Dels'], sgRNAefficiency], part, [k for k,v in hits.items() if hit in v])
                            appendList([pam_efficiency], pam_seq, [k for k,v in hits.items() if hit in v])
                            hits_del.append(hit)
                            break
                    else:
                        break
                except:
                    break
    return(PfD)


#----------------------------------------------------------------------------------------------------------------------------------------
# Checking if the insertion is a CRISPR InDel
#----------------------------------------------------------------------------------------------------------------------------------------

def workingwithIns(Pi, Piq, part, m, z, w, hit, hits, query, ref, indels, sgRNAefficiency, pam_efficiency):
    PiI = Pi # Initial position of Ins
    pam_seq = hit[2]
    if PiI > z and PiI < w:
        appendList([indels['Ins'], sgRNAefficiency], part, [k for k,v in hits.items() if hit in v])
        appendList([pam_efficiency], pam_seq, [k for k,v in hits.items() if hit in v])
    else: # check if there is an InDel moving the alignment
        num = 0
        while True:
            if PiI < z:
                try:
                    if query[Piq+num-1] == ref[PiI+num-1]:
                        num += 1
                        if PiI+num > z and PiI+num < w:
                            appendList([indels['Ins'], sgRNAefficiency], part, [k for k,v in hits.items() if hit in v])
                            appendList([pam_efficiency], pam_seq, [k for k,v in hits.items() if hit in v])
                            break
                    else:
                        break
                except:
                    break
            else: # insertion in the right of sgRNA
                try:
                    if query[Piq+m-num-2] == ref[PiI-num-2]:
                        num += 1
                        if PiI-num > z and PiI-num < w:
                            appendList([indels['Ins'], sgRNAefficiency], part, [k for k,v in hits.items() if hit in v])
                            appendList([pam_efficiency], pam_seq, [k for k,v in hits.items() if hit in v])
                            break
                    else:
                        break
                except:
                    break


#----------------------------------------------------------------------------------------------------------------------------------------
# Check if a big fragment covering two different sgRNAs was deleted
#----------------------------------------------------------------------------------------------------------------------------------------

def extractionFragment(hits_del, hits_epitopes, part, extraction):
    condition1 = 'False' # for the same Dels, there are more than one hit: affects more than one sgRNA hit (maybe the same sgRNA twice, or two different sgRNAs, ...)
    condition2 = 'False' # one sgRNA has to be further to the left than the leftmost hit of epitopes (not DQ2.5_glia_a3 (not included in the 33-mer))
    condition3 = 'False' # one sgRNA has to be further to the right than the rightmost hit of epitopes, not taking into account the DQ2.5_glia_a3 (not included in the 33-mer)
    
    # Check condition1
    if len(hits_del) > 1: condition1 = 'True'
    
    # Calculate the Alpha type
    alphatype = 'Alpha' + str(sum([len(item) for item in hits_epitopes]))

    # Processing hits_epitopes
    pro_hits_epitopes = []
    for hit in hits_epitopes:
        if len(hit) == 0:
            pro_hits_epitopes.append([0])
        else:
            pro_hits_epitopes.append(hit)

    # Check condition2 and condition3
    # The value was multiplied by 3 because the conversion from peptide to DNA
    eps_first_hit = min(list(itertools.chain(*pro_hits_epitopes[0:3])))*3 # leftmost hit of epitopes, DQ2.5_glia_a3 has not been take into account (not included in the 33-mer)
    eps_last_hit = max(list(itertools.chain(*pro_hits_epitopes[0:3])))*3 # rightmost hit of epitope (not DQ2.5_glia_a3)
    for hit in hits_del:
        z = hit[0] # initial position of the sgRNA
        w = hit[1] # last position of the sgRNA
        # Because sgRNAs are 1-based and epitopes hits are 0-based: z-1 and w-1
        if z-1 < eps_first_hit: condition2 = 'True' # if there are no epitopes (a1a, a1b or a2): eps_first_hit = 0
        if w-1 > eps_last_hit and eps_last_hit != 0: condition3 = 'True' # if there are no epitopes (a1a, a1b or a2): eps_last_hit = 0
    
    if condition1 == 'True' and condition2 == 'True' and condition3 == 'True':
        extraction.append((alphatype, part))


#----------------------------------------------------------------------------------------------------------------------------------------
# Search dsOligo
#----------------------------------------------------------------------------------------------------------------------------------------

def search_dsoligo(query, cigar, rname, Pi, long_target, short_target, hits, otu, dsoligo, num_diff, indels_dsoligo, sgRNAefficiency, pam_efficiency):    
    sgRNAs_list = [] # list of sgRNAs in which the insertion of dsoligo occurred

    # Condition for searching dsOligo
    strand = 'NA'
    # if short target is perfect in query
    # It sould be taken into account that dsODN is searched directly into the query present in the sam file: this query is represented always as the 'fwd' versioin of the otu compared to the Amp ref used for the alignment.
    if len(naive_mm(short_target, query, 0)) > 0: strand = 'fwd'
    if len(naive_mm(Seq(short_target).reverse_complement(), query, 0)) > 0: strand = 'rev'
    
    dsoligo_hits = []
    if strand == 'fwd': dsoligo_hits = naive_mm(long_target, query, num_diff)
    if strand == 'rev': dsoligo_hits = naive_mm(Seq(long_target).reverse_complement(), query, num_diff)
    
    dsoligo_status = 'None'
    if len(dsoligo_hits) > 0:
        for hit in dsoligo_hits: # search long target in query with 3 mismatches
            Pref_start_dsoligo = 0
            Pquery_dsoligo,  Pq  = hit+1, 1 # Pquery_dsoligo and Pq, as Pi0, z and w: 1-based
            parts_dict = {'M': [0,0], 'S': [0,0], 'I': [0,0], 'D': [0,0], 'H': [0,0]}
            for part in re.findall(r'\d+[A-Z]', cigar):  
                if Pq < Pquery_dsoligo:
                    parts_dict[re.findall(r'[a-zA-Z]+', part)[0]][0] += int(re.findall(r'\d+', part)[0])
                else:
                    parts_dict[re.findall(r'[a-zA-Z]+', part)[0]][1] += int(re.findall(r'\d+', part)[0])
                if 'D' not in part: Pq += int(re.findall(r'\d+', part)[0])
            Pref_start_dsoligo = (Pi-1) - parts_dict['S'][0] - parts_dict['I'][0] + parts_dict['D'][0] + Pquery_dsoligo # position of dsoligo in ref, to compare with z and w
            Pref_end_dsoligo = Pref_start_dsoligo - parts_dict['I'][1] + parts_dict['D'][1] + len(long_target)

            for guide, positions in hits.items():
                for pos in positions:
                    z, w, pam_seq = pos[0], pos[1], pos[2]
                    if Pref_start_dsoligo < w and z <= Pref_end_dsoligo:
                        if guide not in sgRNAs_list: sgRNAs_list.append(guide)
                        dsoligo_status = 'hit'
                        appendList([indels_dsoligo['Ins'], sgRNAefficiency], str(len(long_target))+'I', [guide]) # add both sgRNAs
                        appendList([pam_efficiency], pam_seq, [guide])
            
    if dsoligo_status == 'hit':
        dsoligo[rname].append((otu, strand, sgRNAs_list))

#----------------------------------------------------------------------------------------------------------------------------------------
# Analysis of sam files, looking for indels and dsOligo
#----------------------------------------------------------------------------------------------------------------------------------------

def analysisSamv2(samFile, otus_sample, metafile, refsFasta, sample, hits_epitopes, dsOligo_status, num_diff):
    sam = open(samFile, 'r')
    otus_wt, otus_indels, otus_extraction, otus_all, otus_undetermined, otus_sgRNAs, otus_dsoligo, length_indels, otus_pam = [], defaultdict(list), defaultdict(list), [], [], defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
    # metafile contains plasmid
    metaFile = pd.read_csv(metafile, sep="\t", header=None)
    plasmid = metaFile[metaFile[0] == sample][3].values[0]
    refs = SeqIO.to_dict(SeqIO.parse(refsFasta, 'fasta'))
    otus_per_ref = defaultdict(list) # Generate a dict with otus per ref to use in stats per amp type
    long_target = metaFile[metaFile[0] == sample][1].values[0] # to search dsOligo
    short_target = metaFile[metaFile[0] == sample][2].values[0] # to search dsOligo
    if plasmid != 'na':
        for line in sam.readlines():
            if not line.startswith('@'):
                otu = line.split("\t")[0].split(";")[0].strip()
                rname = line.split("\t")[2] # Name of the reference, if is equal to '*', the otu is not mapped
                cigar = line.split("\t")[5]
                query = line.split("\t")[9]
                Pi = int(line.split("\t")[3]) # initial position of M in ref
                Piq = 1 # initial position of query
                if otu in otus_sample: # if otu is present in the sample (after apply otus filter in searching step of Usearch)
                    if rname != '*': # if otu is mapped
                        ref = str(refs[rname].seq)
                        otus_per_ref[rname].append(otu) # append otus to the ref
                        if otu not in otus_all: otus_all.append(otu) # all the otus in sam file present in the sample
                        
                        # Position of sgRNAs in the ref, taking into account the plasmid used. Retrieve a list of hits (hits_list) and a list of hits per sgRNA (hits)
                        hits_list, hits = sgRNA(str(refs[rname].seq), plasmid) # (z, w) positions of sgRNAs in ref sequence
                        
                        indels, indels_dsoligo, extraction, sgRNAefficiency, pam_efficiency = defaultdict(list), defaultdict(list), [], defaultdict(list), defaultdict(list) # Lists with otus will be filled

                        # Search dsOligo
                        if dsOligo_status == 'True' and long_target != 'NA' and short_target != 'NA':
                            search_dsoligo(query, cigar, rname, Pi, long_target, short_target, hits, otu, otus_dsoligo, num_diff, indels_dsoligo, sgRNAefficiency, pam_efficiency)

                        # if the otu aligned does not present InDels. Add to otus wt, and remove otu from the other list (indels lists) if it would be there
                        if all(item not in cigar for item in ['D', 'I']):
                            if otu not in [part_ds[0] for item_ds in otus_dsoligo.values() for part_ds in item_ds] and len(indels_dsoligo['Ins']) == 0: # if this otu does not contain InDels and not contain dsoligo
                                otus_wt.append(otu)
                                if otu in otus_indels['Dels']: otus_indels['Dels'].remove(otu)
                                if otu in otus_indels['Ins']: otus_indels['Ins'].remove(otu)
                                for k, v in otus_extraction.items():
                                    if otu in v: otus_extraction[k].remove(otu)
                                for k, v in otus_sgRNAs.items():
                                    if otu in v: v.remove(otu)
                                for k, v in otus_pam.items():
                                    if otu in v: v.remove(otu) # if the otu does not present InDels, it will be removed from otus_pam
                                if otu in otus_undetermined: otus_undetermined.remove(otu)

                        # if the otu aligned presents InDels: Calculate (PiD/PiI, PfD) positions of InDels
                        if any(item in cigar for item in ['D', 'I']):
                            
                            # Looking for true InDels in otus
                            if len(hits_list) > 0: # there are sgRNAs in the reference
                                for part in re.findall(r'\d+[A-Z]', cigar): # loop for the parts in cigar
                                    if any(indel in part for indel in ['I', 'D']): # if the part is an InDel: start
                                        m = int(re.findall(r'\d+', part)[0]) # the number item in this part of the cigar, used to calculate PfD and Piq
                                        # we are going to iterate over all the hits to take into account all the sgRNAs and calculate their efficiency
                                        hits_del = [] # list of hits of sgRNAs with the deletion, take into account for extractionFragment()
                                        for hit in hits_list: # there will be redundancy of otus with InDels if there are InDels in more than one sgRNA, but in indels.py the redundancy is removed with set()
                                            z = hit[0] # initial position of the sgRNA
                                            w = hit[1] # final position of the sgRNA
                                            if 'D' in part: # deletion
                                                PfD = workingwithDels(Pi, Piq, part, m, z, w, hit, hits, query, ref, indels, sgRNAefficiency, hits_del, pam_efficiency)
                                            else: # insertion
                                                workingwithIns(Pi, Piq, part, m, z, w, hit, hits, query, ref, indels, sgRNAefficiency, pam_efficiency)
                                        if 'D' in part:
                                            Pi = PfD # in case of deletions: Pi is equal to y (Pi (=PiD) + digit of deletion)
                                            extractionFragment(hits_del, hits_epitopes[rname], part, extraction) # if the Del is an extraction of the 33-mer (completed or incompleted) region
                                        else: Piq += m # only if part is not a Del, Piq increase
                                    else: # if the part is not an InDel, Pi and Piq increase but any InDel is checked
                                        if 'M' in part:
                                            Pi += int(re.findall(r'\d+', part)[0])
                                        if 'H' not in part:
                                            Piq += int(re.findall(r'\d+', part)[0])
                            
                            # Check if the otu with InDels are already present in otus wt list: if they are actually a wt otu, they are remove from the InDel list (WT has preference)
                            if any(len(item) > 0 for item in list(indels.values())): # if there is any InDel in the otu
                                if otu not in otus_wt:
                                    if len(indels['Dels']) > 0:
                                        otus_indels['Dels'].append(otu)
                                        for indel_event in indels['Dels']: length_indels[otu].append(indel_event)
                                    if len(indels['Ins']) > 0:
                                        otus_indels['Ins'].append(otu)
                                        for indel_event in indels['Ins']: length_indels[otu].append(indel_event)
                                    if len(extraction) > 0:
                                        for element in extraction:
                                            alphatype = element[0] # the Amp types in which the fragment extraction were found for this otu
                                            otus_extraction[alphatype].append(otu)
                                    for k, v in sgRNAefficiency.items():
                                        if len(v) > 0: otus_sgRNAs[k].append(otu)
                                    for k, v in pam_efficiency.items():
                                        if len(v) > 0: otus_pam[k + '_' + ';'.join(v)].append(otu)
                            else: # if the otu has not InDels in the sgRNAs and it is not a wt otu: this is an undetermined otu, an otu with InDel but not by CRISPR
                                if otu not in otus_wt:
                                    otus_undetermined.append(otu)
                        
                        if any(len(item) > 0 for item in list(indels_dsoligo.values())): # if there is dsoligo in the otu, depite it does not contain any InDel
                            if otu in otus_wt: otus_wt.remove(otu)
                            if otu in otus_undetermined: otus_undetermined.remove(otu)
                            otus_indels['Ins'].append(otu)
                            if not any(item in cigar for item in ['D', 'I']):
                                otus_indels['Dels'].append(otu) # In the case this is a perfect substitution, we consider there was a previous excision before the dsoligo insertion
                                for indel_event in indels_dsoligo['Ins']: length_indels[otu].append(indel_event) # as there are no insertion length information about this type of otu without 'InDels', we include the length of dsoligo as length insertion value
                            for k, v in sgRNAefficiency.items():
                                if len(v) > 0: otus_sgRNAs[k].append(otu)
                            for k, v in pam_efficiency.items():
                                if len(v) > 0: otus_pam[k + '_' + ';'.join(v)].append(otu)
    return(otus_wt, otus_indels, otus_extraction, otus_all, otus_undetermined, otus_sgRNAs, otus_per_ref, otus_dsoligo, length_indels, otus_pam)


#----------------------------------------------------------------------------------------------------------------------------------------
# Calculate the frequency of missing epitopes per sample
#----------------------------------------------------------------------------------------------------------------------------------------
def missingEpitopes(otus_indels, bwotus_per_ref, bbotus_per_ref, hits_epitopes, refs, amp_types, sample, otus_sample, otutable, epitopesFile, otus_denoised):
    aligner = Align.PairwiseAligner()
    epitopesdict = SeqIO.to_dict(SeqIO.parse(epitopesFile, 'fasta'))

    # Get the otus per ref
    otus_per_ref = {}
    for rname in nonRedundant(reduce(lambda i,j:i+j,list(amp_types.values()),[])):
        otus_per_ref[rname] = nonRedundant(bwotus_per_ref[rname] + bbotus_per_ref[rname])
    
    # Calculate the frequency of missing epitopes in the sample: only the otus with CRISPR InDels are considered as otus with putative missing epitopes
    otus_missing, total_otus = defaultdict(list), defaultdict(list)
    for rname, otus in otus_per_ref.items():
        rname_seq = refs[rname].seq
        while len(rname_seq)%3 != 0: rname_seq = rname_seq + Seq('N')
        rname_pep = str(rname_seq.translate()) # translate the reference to check below the translation of the otu
        
        # iterate otus to get the reads of missing epitopes
        for otu in [item for item in otus if item in otus_sample]:
            ab_otu = int(otutable[otutable['#OTU ID'] == otu][sample].iloc[0])

            # iterate epitopes
            count = 0
            for ep_name, ep_seq in epitopesdict.items():
                Nep_ref = len(hits_epitopes[rname][count]) # number of hits of this epitope in the reference
                total_otus[ep_name].append(ab_otu*Nep_ref) # abundance of this epitope based on abundance of otu and number of epitope hits

                # if otu has CRISPR indels: check if the epitope hits are missing
                if otu in otus_indels:
                    # Check number of missing epitopes
                    otu_seq = otus_denoised[otu].seq
                    best_alignment = ('', 0.0)
                    for frame in range(3): # try the translation in different frames: select the one that match better with the Amp reference peptide
                        otu_seq_frame = otu_seq[frame:]
                        while len(otu_seq_frame)%3 != 0: otu_seq_frame = otu_seq_frame + Seq('N')
                        otu_pep = str(otu_seq_frame.translate())
                        score = aligner.score(rname_pep, otu_pep)
                        if score > best_alignment[1]: best_alignment = (otu_pep, score)
                    query_pep = best_alignment[0]
                    Nep_amps = len(naive_mm(ep_seq.seq, query_pep, 0))
                    diff = Nep_ref - Nep_amps # diff value can be negative if the amp has more epitopes than the ref
                    if diff < 0: diff = 0 # only the reads with less epitopes in amp than in the ref are considered reads with missing epitopes, so diff can't be less than zero
                    otus_missing[ep_name].append(ab_otu*diff)
                count += 1

    # Calculate the frequency of missing epitopes
    missing_freq = {}
    for ep_name, ep_seq in epitopesdict.items():
        totalreads = sum(total_otus[ep_name])
        missing_reads = sum(otus_missing[ep_name])
        if totalreads != 0:
            freq_missing = float(missing_reads/totalreads)*100
        else:
            freq_missing = 0.0
        missing_freq[ep_name] = freq_missing
    return(missing_freq)

#----------------------------------------------------------------------------------------------------------------------------------------
# Calculate the frequency of epitopes in sample
#----------------------------------------------------------------------------------------------------------------------------------------
def totalEpitopes(bwotus_per_ref, bbotus_per_ref, refs, amp_types, sample, otus_sample, otutable, epitopesFile, otus_denoised):
    aligner = Align.PairwiseAligner()
    epitopesdict = SeqIO.to_dict(SeqIO.parse(epitopesFile, 'fasta'))

    # Get the otus per ref
    otus_per_ref = {}
    for rname in nonRedundant(reduce(lambda i,j:i+j,list(amp_types.values()),[])):
        otus_per_ref[rname] = nonRedundant(bwotus_per_ref[rname] + bbotus_per_ref[rname])
    
    # Found epitopes in Amps
    ep_freqs = defaultdict(list)
    total_reads = otutable[sample].sum() # otutable filtered by otus aligned to BW208 Amps
    for otu in otus_sample:
        # Get the peptide sequence of the Amp ref
        rnames = [k for k, v in otus_per_ref.items() if otu in v]
        if len(rnames) > 0:
            rname = rnames[0] # one ref as example
            rname_seq = refs[rname].seq
            while len(rname_seq)%3 != 0: rname_seq = rname_seq + Seq('N')
            rname_pep = str(rname_seq.translate()) # translate the reference to check below the translation of the otu

            # Abundance of otu
            ab_otu = int(otutable[otutable['#OTU ID'] == otu][sample].iloc[0])

            # Get the otu peptide
            otu_seq = otus_denoised[otu].seq
            best_alignment = ('', 0.0)
            for frame in range(3): # try the translation in different frames: select the one that match better with the Amp reference peptide
                otu_seq_frame = otu_seq[frame:]
                while len(otu_seq_frame)%3 != 0: otu_seq_frame = otu_seq_frame + Seq('N')
                otu_pep = str(otu_seq_frame.translate())
                score = aligner.score(rname_pep, otu_pep)
                if score > best_alignment[1]: best_alignment = (otu_pep, score)
            query_pep = best_alignment[0]

            # Get the epitopes in otu
            for ep_name, ep_seq in epitopesdict.items():
                Nep_amps = len(naive_mm(ep_seq.seq, query_pep, 0))
                ep_freqs[ep_name].append(float(ab_otu/total_reads)*100*Nep_amps)
        
        # Get the epitope frequencies
        ep_freq = {}
        for ep_name, ep_seq in epitopesdict.items():
            ep_freq[ep_name] = sum(ep_freqs[ep_name])
    return(ep_freq)


#----------------------------------------------------------------------------------------------------------------------------------------
# Calculate the frequency of each indel event based on their length
#----------------------------------------------------------------------------------------------------------------------------------------
def length_indel_events(bwlength_indels, bblength_indels, otus_indels, sample, otutable):
    ab_indel_events = {}
    totalreads = otutable[sample].sum()
    for otu in otus_indels: # preference: long length indels, so if the otu is in bblength_indels, it has preference
        ab_otu = int(otutable[otutable['#OTU ID'] == otu][sample].iloc[0])
        if otu in list(bblength_indels.keys()):
            for item in nonRedundant(bblength_indels[otu]):
                ab_item = float(ab_otu/totalreads)*100
                ab_indel_events[item] = ab_indel_events.get(item, 0) + ab_item
        else:
            for item in nonRedundant(bwlength_indels[otu]):
                ab_item = float(ab_otu/totalreads)*100
                ab_indel_events[item] = ab_indel_events.get(item, 0) + ab_item
    return(ab_indel_events)


#----------------------------------------------------------------------------------------------------------------------------------------
# Calculate the efficiency of PAM per guide COMPROBAR 12/07/2024
#----------------------------------------------------------------------------------------------------------------------------------------
def PAMefficiency(bwpam_efficiency, bbpam_efficiency, otus_indels, sample, otutable, allotus_perPAM_sgRNA):
    ab_pam = {}
    
    totalreads_sgRNAPAM_dict = {}
    for k, v in allotus_perPAM_sgRNA.items():
        otutable['Status' + k] = ['sgRNA+PAM' if otu in v else 'no' for otu in otutable['#OTU ID']]
        otutable_stats_sgRNAPAM = otutable[otutable['Status' + k] == 'sgRNA+PAM']
        totalreads_sgRNAPAM = otutable_stats_sgRNAPAM[sample].sum()
        totalreads_sgRNAPAM_dict[k] = totalreads_sgRNAPAM
    
    # Calculate edition efficiency of each sgRNA+PAM combination
    for otu in otus_indels: # preference: long length indels, so if the otu is in bblength_indels, it has preference
        ab_otu = int(otutable[otutable['#OTU ID'] == otu][sample].iloc[0])
        if len([k for k, v in bbpam_efficiency.items() if otu in v]) > 0:
            guidePAM_list = []
            for item in [k for k, v in bbpam_efficiency.items() if otu in v]: # the otu is in bbmap, this has preference
                guide = item.split('_')[0]
                pams = item.split('_')[1].split(';')
                for pam in list(set(pams)):
                    guidePAM_list.append(guide + '+' + pam)
            for guidePAM in list(set(guidePAM_list)):
                totalreads_pam = totalreads_sgRNAPAM_dict[guidePAM]# totalreads of all otus mapped containing this PAM+sgRNA
                ab_item = float(ab_otu/totalreads_pam)*100
                ab_pam[guidePAM] = ab_pam.get(guidePAM, 0) + ab_item
        else:
            guidePAM_list = []
            for item in [k for k, v in bwpam_efficiency.items() if otu in v]: # the otu is not in bbmap, in bwa instead
                guide = item.split('_')[0]
                pams = item.split('_')[1].split(';')
                for pam in list(set(pams)):
                    guidePAM_list.append(guide + '+' + pam)
            for guidePAM in list(set(guidePAM_list)):        
                totalreads_pam = totalreads_sgRNAPAM_dict[guidePAM]# totalreads of all otus mapped containing this PAM+sgRNA
                ab_item = float(ab_otu/totalreads_pam)*100
                ab_pam[guidePAM] = ab_pam.get(guidePAM, 0) + ab_item
    return(ab_pam)