from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
	description = "makes output file with chromosome name and length"
	)
parser.add_argument("--seq", dest = "seq", help = "TE family file (fasta)",
	required = True)
parser.add_argument("--outstem", dest = "outstem", help = 
	"assigned name for output file", required = True)

def fam_len(record):
	""""""
	for rec in seq:
		length = len(rec)
		fam = rec.id
		yield([fam, length]) 

if __name__ == "__main__":
	args = parser.parse_args()
	seq = SeqIO.parse(args.seq, "fasta")
	
	with open(args.outstem + ".tsv", "w") as out:
		out.write("Family\tLength\n")
		for fam,length in fam_len(seq): 
			out.write("{}\t{}\n".format(fam, length))
