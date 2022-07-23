# This folder contains all custom python scripts written for this project. 

## AT_RIP_IDX.py

#### Description
Calculates RIP indices for AT rich regions in the genome. Required input is genomic sequences in fasta format, and a user-defined name for the output files (outstem).
The script will create two output files in .tsv format - one with dinucleotide frequencies of each sequence, and one with RIP indices.

#### Usage
Fasta sequences of AT rich regions can be obtained with BEDTOOLS getfasta using a .bed file from OcculterCut coordinates and the reference genome as input.

`bedtools getfasta -fi [reference genome].fna -fo [AT_outstem].fna -bed [AT_coordinates].bed`

`python AT_RIP_IDX.py --seqs [AT_outstem].fna --outstem [outstem]`

## RIP Dinucs.py

#### Description
This script achieves the same outcome as above, but was tailored for 

#### Usage

## chrom_lens.py
#### Description
Creates a tsv file of chromosome lengths

#### Usage
`python chrom_lens.py --seq [reference_genome].fasta --outstem outstem`

## fam_lens.py
#### Description
Creates a .tsv file of TE copies and length

#### Usage
`python fam_lens.py --seq [TE_seqs].fasta --outstem outstem`

## fix_con_x.py
#### Description
Resolves nucleotide ambiguities in the final stage of manual curation

#### Usage

