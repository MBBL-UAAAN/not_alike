#!/usr/bin/env python3

#   This module will update the transfered file from master
#   which contains the input sequences. This module will delete
#   all sequences that had a hit in previous blast searching.

from sys import argv, exit


def store_seqs(filename):
    """
        Stores the sequences of input file
        in a dictionary data sctructure.
        Requires the file name.
    """
    seq = {}
    FHIN = open(filename, "r")
    for line in FHIN:
        line = line.strip()
        if line[0] == ">":
            head = line[1:]
            seq[head] = ""
        else:
            seq[head] += line
    FHIN.close()
    return seq

def update_seqs(seq_dicty, guid_file, guid_hit_file):
    """
        This function takes sequences that had no hit in
        previous blast searching.
    """
    
    FHIN_GUID = open(guid_file, "r")
    x = []
    for line in FHIN_GUID:
        line = line.strip()
        x.append(line)

    FHIN_GUID.close()

    FHIN_GUID_HIT = open(guid_hit_file, "r")
    y = []
    for line in FHIN_GUID_HIT:
        line = line.strip()
        y.append(line)

    FHIN_GUID_HIT.close()

    r = set(x).difference(set(y))
    r = list(r)
    outseq_dicty = {}
    for h, s in seq_dicty.items():
        if h in r:
            outseq_dicty[h] = s


    return outseq_dicty

def main():
    """
        Main function
    """

    seqs_filename       =   argv[1]
    guid_filename       =   argv[2]
    guid_hit_filename   =   argv[3]
    
    seq_dicty = store_seqs(seqs_filename)
    useq_dicty = update_seqs(seq_dicty, guid_filename, guid_hit_filename)
    for h, s in useq_dicty.items():
        print(">" + h + "\n" + s)


if __name__ == "__main__":
    main()
