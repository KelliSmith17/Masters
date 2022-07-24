from collections import Counter
import argparse
import sys
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


parser = argparse.ArgumentParser(
    description = ""
    )

parser.add_argument("--alignment", dest = "alignment", help = "alignment in .fa")

def is_primary_seq(rec):
    return(rec.id[0] in ["r", "_", "E", "l"])

def base_consensus(base_list):
    """ """
    
    base_c = Counter(base_list)
    max_count = max(base_c.values())
    most_common_base  = [b for b in base_c.items() if b[1] == max_count]
   
    if len(most_common_base) == 1:
        return(most_common_base[0][0])
    tied_bases = [base for base,count in most_common_base]
    if len(tied_bases) > 1:
        #print(tied_bases)
        if base_list[0] != "-":
            RM_con = base_list[0]
        else:
            RM_con = "N"
      
        tb = sorted(tied_bases)
        if tb == ['A', 'C', 'G', 'T']:
            return(RM_con)
        elif tb == ['A', 'C', 'G']:
            return(RM_con)
        elif tb == ['A', 'C', 'T']:
            return("C")
        elif tb == ['A', 'G', 'T']:
            return("G")
        elif tb == ['C', 'G', 'T']:
            return(RM_con)
        elif tb == ['A', 'C']:
            return("C")
        elif tb == ['A', 'G']:
            return("G")
        elif tb == ['A', 'T']:
            return(RM_con) 
        elif tb == ['C', 'G']:
            return(RM_con)
        elif tb == ['C', 'T']:
            return("C")
        elif tb == ['G', 'T']:
            return("G")
       
        elif tb == ['-','A', 'C', 'G', 'T']:
            return(RM_con) 
        elif tb == ['-','A', 'C', 'G']:
            return(RM_con)
        elif tb == ['-','A', 'C', 'T']:
            return("C")
        elif tb == ['-','A', 'G', 'T']:
            return("G")
        elif tb == ['-','C', 'G', 'T']:
            return(RM_con)
        elif tb == ['-','A', 'C']:
            return("C")
        elif tb == ['-','A', 'G']:
            return("G")
        elif tb == ['-','A', 'T']:
            return(RM_con) 
        elif tb == ['-','C', 'G']:
            return(RM_con)
        elif tb == ['-','C', 'T']:
            return("C")
        elif tb == ['-','G', 'T']:
            return("G")
        elif tb == ['-', 'A']:
            return("A")
        elif tb == ['-','C']:
            return("C")
        elif tb == ['-', 'G']:
            return("G")
        elif tb == ['-', 'T']:
            return("T")


if __name__ == "__main__":
    args = parser.parse_args()
    ali = SeqIO.to_dict(SeqIO.parse(args.alignment, "fasta"))
    print("working on", args.alignment)
    if not "CON" in ali.keys():
        print("There is no seq named 'CON'") 
        sys.exit(1)



    new_seq = ""
    for base_index, base in enumerate(ali["CON"]):
        if base == "X":
            #print("=====one to fix!====\n") # for testing
            #print(base_index,base)
            bases =([r[base_index] for r in ali.values() if is_primary_seq(r)])
            new_seq += base_consensus(bases)
           
        else:
            new_seq += base
            #print(base_index, base)
           

        new_consensus = SeqRecord(id = "NEW_CON", seq = Seq(new_seq))
        ali["CON"] = new_consensus
       
        SeqIO.write(ali.values(), args.alignment + "_filled.fa", "fasta")
    print("done")
    sys.exit(0)
