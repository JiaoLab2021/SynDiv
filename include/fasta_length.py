#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os
import sys

# include the path of the file_open.py
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + os.sep)
from file_open import MyOpener

# Get the length of each sequence in fasta
def getLength(fastaFilePath, outputFilePath):
    # Create a dictionary to store the length of each sequence
    length_dict = {}
    seq = ""
    name = ""

    fileOpener = MyOpener(fastaFilePath)

    # Open the fasta file
    line = [None]
    while fileOpener.read_line(line):
        info = line[0].strip()

        # skip empty lines
        if not info:
            continue

        if info.startswith('>'):
            if seq != "":
                # Count the length of the sequence and save it in the dictionary
                length_dict[name] = len(seq)
                seq = ""
            name = info[1:]
        else:
            seq += info
    # process the last sequence
    length_dict[name] = len(seq)
        
    # Output the ID and length of each sequence
    with open(outputFilePath, "w") as f:
        for seq_id, seq_len in length_dict.items():
            f.write(f"{seq_id}\t{seq_len}\n")

    return 0
