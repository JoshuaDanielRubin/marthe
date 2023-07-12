#!/usr/bin/env python

import os
import sys
import re
from tqdm import tqdm

def calculate_map_score(matches):
    #Calculate map score
    if len(matches)!=0:
        return 1/len(matches)
    else:
        return -1

def find_matches(string, seq):
    #Find matches of a given string in a sequence
    return re.findall(string, seq)

def read_reference_file(file_path):
    #Read a reference file and return its content as a string
    with open(file_path, "r") as ref_file:
        seq = "".join(line.strip() for line in ref_file)
    return seq

def process_bam_file(bam_file_path, seq):
    #Process a bam file and print map scores
    with open(bam_file_path, "r") as bam_file:
        for line in bam_file:
            pos = line.split()[3]
            string = str(line.split()[9])
            matches = find_matches(string, seq)
            print(calculate_map_score(matches))

def main(ref_file_path, bam_file_path):
    seq = read_reference_file(ref_file_path)
    process_bam_file(bam_file_path, seq)

if __name__ == "__main__":
    ref_file_path = "/home/databases/genomes/Homo_sapiens/hs37d5/hs37d5.fa"
    bam_file_path = "/home/projects/marthe/tmp/test.bam"
    main(ref_file_path, bam_file_path)
