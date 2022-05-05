#!/usr/bin/env python3

from sys import argv, exit

usage = """
USAGE
        split-genome <input_file> <fragment_size> <step_size> <output_file> <separate_file> <num_of_workers>
"""

if len(argv) != 7:
    print("Error at parsing argumets.\nI expect 6 arguments")
    print(usage)
    exit(0)
else:
    infile      =   argv[1]
    frg_size    =   int(argv[2])
    stp_size    =   int(argv[3])
    outfile     =   argv[4]
    sep_file    =   bool(int(argv[5]))
    workers     =   int(argv[6])

seqs = {}
with open(infile, "r") as FHIN:
    for line in FHIN:
        line = line.strip("\n")
        if line[0] == ">":
            header = line
            seqs[header] = ""
        else:
            seqs[header] += line

    FHIN.close()

frg_seqs = {}
for header, sequence in seqs.items():
    len_seq = len(sequence)
    pivot = 0
    while pivot <= len(sequence) - frg_size:
        name = header + "_" + str(pivot)
        frg_seqs[name] = sequence[pivot : (pivot + frg_size)]
        pivot += stp_size

if sep_file:
    num_of_seqs = len(frg_seqs)
    working_batch = num_of_seqs // workers
    working_batch_num = 0
    counter = 0
    tmp_string = ""
    for h, s in frg_seqs.items():

        tmp_string += h + "\n" + s + "\n"
        counter += 1

        if counter == working_batch:
            FHOUT = open(outfile + "_" + str(working_batch_num) + ".fasta", "w")
            FHOUT.write(tmp_string)
            FHOUT.close()
            working_batch_num += 1
            counter = 0
            tmp_string = ""

    if tmp_string != "":
        FHOUT = open(outfile + "_" + str(working_batch_num) + ".fasta", "w")
        FHOUT.write(tmp_string)
        FHOUT.close()
else:
    big_string = ""
    for h,s in frg_seqs.items():
        big_string += h + "\n" + s + "\n"

    with open(outfile, "w") as FHOUT:
        FHOUT.write(big_string)
        FHOUT.close()

