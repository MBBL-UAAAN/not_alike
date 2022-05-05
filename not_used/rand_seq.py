#!/usr/bin/env python3

from random import sample
from sys import argv, exit

def rand_seq(c=200, s=600):
	"""
		Generates a random sequence of nitrogen bases
		c:= number of each element in sampling population
		s:= Size of sequence
	"""

	bases = ["A", "C", "G", "T"]
	rs = "".join(sample(bases, counts=[c,c,c,c], k=s))
	return rs

def create_seqs_str(num_of_seqs=100, counts=200, size=600):
	"""
		Creates a data structure with "num_of_seqs"
		random sequences
		num_of_seqs:= number of sequences.
	"""

	seqs = {}
	
	seqs_name = [">seq_" + str(i) for i in range(num_of_seqs)]
	for i in range(num_of_seqs):
		seqs[seqs_name[i]] = rand_seq(counts, size)
	return seqs

def write_seqs_to_file(dictionary, outfile):
    """

    """
#    print(dictionary)
    outstring = ""
    for header, sequence in dictionary.items():
        outstring += header + "\n" + sequence + "\n"

    #print(outstring)
    try:
        FHIN = open(outfile, "w")
        FHIN.write(outstring)
        FHIN.close()
        return True
    except Exception:
        return False


def main():
	"""
		Main function
	"""
	usage = """###########################\nUSAGE:\n\t rand_seq.py <num_of_seqs> <counts> <size> <output file>\n###########################\n"""
	if len(argv) != 5:
		print(usage)
		exit(0)
	else:
		num_of_seqs =	int(argv[1])
		counts      =	int(argv[2])
		size    	=	int(argv[3])
		outfile     =	argv[4]

	seqs = create_seqs_str(num_of_seqs, counts, size)
	if write_seqs_to_file(seqs, outfile):
		print("Sequences stored in file")
	else:
		print("Something went wrong")


if __name__ == "__main__":
	main()

