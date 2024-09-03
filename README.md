# 🚀**ampAnalysis**

`ampAnalysis` was developed for the analysis of NGS amplicon deep sequences. It comprised a comprehensive pipeline to cover the entire process from amplicon assembly to insertions/deletions (InDels) detection. This software encompasses four main steps: (i) the obtention of samples’ Amps (denoised amplicons), (ii) the trimming of non-α-gliadin Amps, (iii) the alignment of Amps to the BW208 Amps database, and (iv) the InDels/double-stranded oligodeoxynucleotide (dsODN) analysis.

![ampAnalysis pipeline](https://github.com/MiriamMarinS/ampAnalysis/blob/main/images/Pipeline.png)

For the obtention of sample Amps, raw reads were processed using the Usearch v9.2.64 (Edgar, 2010), with steps including merging (`-fastq_mergepairs`), low-quality filtering (`-fastq_filter`), dereplication (`-fastqx_uniques`), and denoising (`-unoise2`), using option values from Bayesian optimization performed in previous research (Guzmán-López et al., 2021). The resulting denoised amplicons, termed Amps, had their abundance quantified by searching and mapping raw reads to the Amps database (`-search_global`). To filter for α-gliadin Amps, an α-gliadins' sequences database was constructed with sequences from the NCBI nt database. The ultra-fast software MMSeqs2 (Steinegger and Söding, 2017) matched the Amps against this database, removing any non-matching sequences to eliminate off-target products.
Before InDels and dsODN analyses, each Amp sequence was aligned to the BW208 Amps database to identify the closest reference sequence. Two alignment strategies were used: (i) local alignment with BWA-MEM (Li, 2013), a Burrow-Wheeler Aligner for short reads, and (ii) global alignment with BBmap v39.01 (Bushnell, 2014).
| :exclamation:  Equipment used   |
|----------------------------------------------|
|All these steps were executed using the ampAnalysis software on an HPC Cluster with 24 Bull x440 nodes, with 192 GB of RAM per node, totaling 960 cores (Advanced Computing Unit (UCAS), University of Córdoba, Spain).|

## **InDels and dsODN insertion analysis**
The single-guide RNAs (sgRNAs) were searched within the WT Amps to identify their positions if the seed sequence (from protospacer adjacent motif (PAM) to 12 bp upstream of the sequence) and the PAM motif had a perfect match. The positions of InDels for each sample Amp, obtained from the alignment results, were then compared to the relative reference positions of the sgRNAs to classify these mutations as putative CRISPR editions. Given the repetitive nature of the α-gliadin genes, the initial positions of the InDels did not always overlap perfectly with the sgRNAs positions. To address this, re-alignment functions were implemented to refine these sections and reduce the number of false-negative edits. The frequencies of InDels were then calculated based on the Amps abundance data obtained in the initial step of the pipeline. The insertion of dsODN fragments was analyzed using alignment results, taking into account the positions of the sgRNAs in the reference Amp sequences.

## **Run `ampAnalysis`**
```
python ampAnlysisv2.py -i </path/to/raw_fastq/> -f -s -p <project name> -m <minampsize value, int> -a -t <path/to/metafile.txt> -n > log.out
```
> [!NOTE]  
> `-i` <string, /path/to/folder/ of raw fastq files> \
> `-f` If there is more than one sample to process, do not put a prefix in `-i` instead. \
> `-s` Run usearch per sample (parallelize). \
> `-p` <string, name of the project> \
> `-m` <int, minampsize usearch parameter value> \
> `-t` <string, /path/to/metafile.txt, example in /examples/metafile.txt> \
> `-c` <int, number of differences allowed in dsoligo search> \
> `-d` Search dsOligo in denoised amplicons. \
> `-a` Removed denoised amplicons not annotated as alpha-gliadins. \
> `-n` Search indels in denoised amplicons.

| :warning: Required softwares           |
|:----------------------------|
| Usearch     |
| MMSeqs2 |
| BWA-MEM |
| BBmap |
| python libraries: biopython, numpy, pandas, progressbar, and re |

*Published in Marín-Sanz et al. (2024): Cas9 and Cas12a-mediated excision and replacement of the celiac disease-related α-gliadin immunogenic complex in hexaploid wheat*

### **References**
*Bushnell, B. (2014) BBMap: a fast, accurate, splice-aware aligner. Lawrence Berkeley National Lab.(LBNL), Berkeley, CA (United States).* \
*Edgar, R. (2010) Usearch. Lawrence Berkeley National Lab.(LBNL), Berkeley, CA (United States).* \
*Guzmán-López, M.H., Marín-Sanz, M., Sánchez-León, S., and Barro, F. (2021a) A bioinformatic workflow for InDel analysis in the wheat multi-copy α-gliadin gene family engineered with CRISPR/Cas9. Int. J. Mol. Sci., 22, 13076.* \
*Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv Prepr. arXiv1303.3997.* \
*Steinegger, M. and Söding, J. (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat. Biotechnol., 35, 1026–1028.*
