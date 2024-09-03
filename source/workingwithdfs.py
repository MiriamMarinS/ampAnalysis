import pandas as pd
from utils import nonRedundant, nonRedundantItem0, calculateStats
import os
import itertools


def creatingDFfromOtus(otusperType, ampType, otus_indels, otus_wt, otus_undetermined, otutable, sample, dsOligo_status, sam_bwa, outputname, otus_dels, otus_extraction, otus_dsoligo_amp):
    stats_per_amp = {}
    # Number of otus
    wt_in_amp = [otu for otu in otus_wt if otu in otusperType]
    dels_in_amp = [otu for otu in otus_dels if otu in otusperType and not otu in wt_in_amp] # Take into account: if a otu assigned to this Amp type is a del here, but a wt in other Amp type, thi otu will be considered as wt for the current Amp type, ignoring its status for the other Amps types 
    # In the case of extraction is different from dels, the extractions are found mainly in BBmap alignments, not BWA-mem, but otus aligned with BWA are assigned also to Amp types, so there are otus with extractions in Alpha7 with BBmap that are also aligned with BWA to Alpha0 amps... Corrected in extractionFragment()
    extraction_in_amp = [otu for otu in otus_extraction[ampType] if otu in otusperType and not otu in wt_in_amp] # if otu in v is redundant...

    # Total reads
    # Filter by otus mapped
    amps_stats = otutable.copy() # only contain otus in sample (and denoised if apply (alpha option))
    amps_stats['Status_mapped'] = ['mapped' if otu in nonRedundant(otus_indels + otus_wt + otus_undetermined) else 'no' for otu in amps_stats['#OTU ID']]
    amps_stats = amps_stats[amps_stats['Status_mapped'] == 'mapped'] # filter otutable by otus mapped with bwa-mem or BBmap
    # Filter by otus present in the Amp type
    amps_stats['Otus_in_amp'] = ['amp' if otu in otusperType else 'no' for otu in amps_stats['#OTU ID']]
    amps_stats = amps_stats[amps_stats['Otus_in_amp'] == 'amp']
    totalreads_in_amp = amps_stats[sample].sum() # sum of reads of all the otus in sample mapped

    # Calculate reads and freq per otu list
    for typeEvent, listEvent in {'wt': wt_in_amp, 'dels': dels_in_amp, 'ext': extraction_in_amp}.items():
        reads_in_amp, freq_in_amp = calculateStats(amps_stats, 'Status_' + typeEvent, listEvent, totalreads_in_amp, sample)
        stats_per_amp['otus_' + typeEvent] = len(listEvent)
        stats_per_amp['reads_' + typeEvent] = reads_in_amp
        stats_per_amp['freq_' + typeEvent] = freq_in_amp

    # Stats for dsoligo
    if dsOligo_status == 'True':
        otus_with_dsoligo = [item[0] for item in nonRedundantItem0(itertools.chain(*[v for k, v in otus_dsoligo_amp.items()]))]
        dsoligo_in_amp = [otu for otu in otus_with_dsoligo if otu in otusperType] # if otu in otusperType is redundant... a dsoligo can be in wt otu
        reads_dsoligo_in_amp, freq_dsoligo_in_amp = calculateStats(amps_stats, 'Status_dsoligo', dsoligo_in_amp, totalreads_in_amp, sample)
        stats_per_amp['otus_dsoligo'] = len(dsoligo_in_amp)
        stats_per_amp['reads_dsoligo'] = reads_dsoligo_in_amp
        stats_per_amp['freq_dsoligo'] = freq_dsoligo_in_amp
    
    # Generate csv per sample
    if dsOligo_status == 'False':
        stats_amp = pd.DataFrame({'Params': ['wt_otus', 'dels_otus', 'ext_otus', 'Reads_wt', 'Reads_dels', 'Reads_ext', 'Freq_wt', 'Freq_dels', 'Freq_ext'],
                                    sample: [len(set(wt_in_amp)), len(set(dels_in_amp)), len(set(extraction_in_amp)),
                                            stats_per_amp['reads_wt'], stats_per_amp['reads_dels'], stats_per_amp['reads_ext'],
                                            stats_per_amp['freq_wt'], stats_per_amp['freq_dels'], stats_per_amp['freq_ext']]})
        stats_amp.to_csv(os.path.dirname(sam_bwa) + '/temp_stats_indels/' + outputname, sep="\t", index=False)
    else:
        stats_amp = pd.DataFrame({'Params': ['wt_otus', 'dels_otus', 'ext_otus', 'dsoligo_otus', 'Reads_wt', 'Reads_dels', 'Reads_ext', 'Reads_dsoligo', 'Freq_wt', 'Freq_dels', 'Freq_ext', 'Freq_dsoligo'],
                                    sample: [len(set(wt_in_amp)), len(set(dels_in_amp)), len(set(extraction_in_amp)), len(set(dsoligo_in_amp)),
                                            stats_per_amp['reads_wt'], stats_per_amp['reads_dels'], stats_per_amp['reads_ext'], reads_dsoligo_in_amp,
                                            stats_per_amp['freq_wt'], stats_per_amp['freq_dels'], stats_per_amp['freq_ext'], freq_dsoligo_in_amp]})
        stats_amp.to_csv(os.path.dirname(sam_bwa) + '/temp_stats_indels/' + outputname, sep="\t", index=False)