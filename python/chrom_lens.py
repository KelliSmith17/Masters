from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
	description = "makes output file with chromosome name and length"
	)
parser.add_argument("--seq", dest = "seq", help = "reference genome (fasta)",
	required = True)
parser.add_argument("--outstem", dest = "outstem", help = 
	"assigned name for output file", required = True)

def chrom_len(record):
	""""""
	for rec in seq:
		length = len(rec)
		chrom = rec.id
		yield([chrom, length]) 

if __name__ == "__main__":
	args = parser.parse_args()
	seq = SeqIO.parse(args.seq, "fasta")
	
	with open(args.outstem + ".tsv", "w") as out:
		out.write("Chromosome\tLength\n")
		for chrom,length in chrom_len(seq): 
			out.write("{}\t{}\n".format(chrom, length))
