from Bio import SeqIO
from collections import Counter
import argparse

#to add in RIP indices
#Rip indices = TpA/ApT and (CpA + TpG)/(ApC + GpT)
#where Rip index 1 = TpA/ApT
#and Rip index 2 = (CpA+TpG)/(ApC/GpT)

parser = argparse.ArgumentParser(
    description = "makes two files with info on dinucs and a file with RIP indices"
    )

parser.add_argument("--seqs", dest = "seqs", help = "put consensus sequences here (fasta)", required = True)
parser.add_argument("--outstem", dest = "outstem", help = "name specified for output files. Indices file will automatically be '[outstem]_indices'")

def count_dinucs(record):
    """count dinucleotides"""
    all_dinucs = []
    for i in range (len(record.seq)-1):
        all_dinucs.append(str(record.seq[i:i+2]))
    dinuc_counted = Counter(all_dinucs)
    return(dinuc_counted)

def calc_RIP_idx(dinuc_counted):
    try:
        RIP_idx_1 = round(dinuc_counted['TA']/dinuc_counted['AT'],4)
    except ZeroDivisionError:         
        RIP_idx_1 = "NaN"
    try:
        RIP_idx_2 = round((dinuc_counted['CA']+dinuc_counted['TG'])/(dinuc_counted['AC']+dinuc_counted['GT']), 4)
    except ZeroDivisionError:
        RIP_idx_2 = "NaN"
    return([RIP_idx_1, RIP_idx_2])
  

if __name__ == "__main__":
    args = parser.parse_args()
    seqs = SeqIO.parse(args.seqs, "fasta")
    R_idx_handle = open(args.outstem + "_indices.tsv", "w")
    R_idx_handle.write("family\tTE_class\tRIP_idx_1\tRIP_idx_2\n")

    with open(args.outstem + ".tsv", "w") as out:
       
        for rec in seqs:
            dn = count_dinucs(rec)
            family,TE_class = rec.id.split("#")
            idx_1, idx_2 = calc_RIP_idx(dn)

            R_idx_handle.write("{}\t{}\t{}\t{}\n".format(family, TE_class, idx_1, idx_2))
            biglist = ['AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC']
            denom = sum([n for (d, n) in dn.items() if "N" not in d]) # to ignore Ns #d= dinucleotide, n = count, take sum of n only if no N in d
            for d in biglist:
                percentage = dn[d]/denom * 100
                out.write("{}\t{}\t{}\t{}\t{:.4f}\n".format(family, TE_class, d, dn[d], percentage))



  

                 

            
