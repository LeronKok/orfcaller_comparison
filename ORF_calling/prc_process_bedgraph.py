#!/usr/bin/env python3

"""Update the coordinates of the P-sites assigned by PRICE

Authors: l.w.kok-15@prinsesmaximacentrum.nl
Date: 01-12-2021

This script reads a bedgraph file containing the P-sites assigned by PRICE. 
The coordinates are then updated such that P-sites are localized to the first
nucleotide of the codon instead of the whole codon.

Input arguments:
    argv[1] - the input .bedgraph file.
    argv[2] - the output .bedgraph file to which the results will be written.
    argv[3] - '+' or '-', to indicate whether the P-sites are those from the 
              positive or the negative strand.
"""

from sys import argv

if __name__ == "__main__":
    count = 0
    if argv[3] != "+" and argv[3] != "-":
        raise ValueError("Third argument should be '+' or '-'")
    with open(argv[1]) as f_open:
        with open(argv[2], "w+") as f_out:
            if argv[3] == '+':
                for line in f_open:
                    s_line = line.strip().split()
                    s_line[2] = str(int(s_line[1])+1)
                    f_out.write("\t".join(s_line) + "\n")
            else:
                for line in f_open:
                    s_line = line.strip().split()
                    s_line[1] = str(int(s_line[1])+2)
                    s_line[2] = str(int(s_line[1])+1)
                    f_out.write("\t".join(s_line) + "\n")