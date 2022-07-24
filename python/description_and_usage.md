# This folder contains all custom python scripts written for this project. 

## AT_RIP_IDX.py

#### Description
Calculates RIP indices for AT rich regions in the genome. Required input is genomic sequences in fasta format, and a user-defined name for the output files (outstem).
The script will create two output files in .tsv format - one with dinucleotide frequencies of each sequence, and one with RIP indices.

#### Usage
Fasta sequences of AT rich regions can be obtained with BEDTOOLS getfasta using coordinates in the OcculterCut gff (or .bed file) and the reference genome as input. The resulting file is then run through the python script. This script can be used for GC rich regions in the same way.

`bedtools getfasta -fi [reference genome.fasta] -fo [AT_outstem.fasta] -bed [AT.gff]`

`python AT_RIP_IDX.py --seqs [AT_outstem].fna --outstem [outstem]`

## RIP_dinucs.py

#### Description
This script achieves the same outcome as above, but was tailored for TE sequences (slight changes in output formating)

#### Usage
`python RIP_dinucs.py --seqs [sequences.fasta] --outstem [outstem]`


## chrom_lens.py
#### Description
Creates a tsv file of chromosome lengths

#### Usage
`python chrom_lens.py --seq [reference_genome.fasta] --outstem outstem`

## fam_lens.py
#### Description
Creates a .tsv file of TE copies and length

#### Usage
`python fam_lens.py --seq [TE_seqs.fasta] --outstem outstem`

## fix_consensus_x.py
#### Description
Resolves nucleotide ambiguities in the new TE consensus during the final stage of manual curation. For the script to work, the new consensus must be named `CON` with ambiguous bases denoted as `X`. The tie-breaking algorithm was established following advice from TE experts.

#### Usage
`python fix_consensus_x.py --alignment [alignment.fasta]`

