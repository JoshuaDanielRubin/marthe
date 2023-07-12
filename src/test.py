#!/usr/bin/env python

import sys
import os
from function import *
import re
from tqdm import tqdm
import mmap

"""def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

ref_file = "/home/databases/genomes/Homo_sapiens/hs37d5/hs37d5.fa"
bam_file = "/home/projects/marthe/tmp/test.bam"
#ref_file = open(sys.argv[1], "r")
#bam_file = open(sys.argv[2], "r")
ref_file = open(ref_file, "r")
bam_file = open(bam_file, "r")
seq =[]

for line in ref_file:
    seq.append(line[:-1])
ref_file.close()
seq="".join(seq)
print("sequence joined")

for line in bam_file:
    pos = line.split()[3]
    string = str(line.split()[9])
    matches = re.findall(string, seq)
    if len(matches)!=0:
        map_score=1/len(matches)
    else: 
        map_score=-1
    print(map_score)
bam_file.close()


#print(seq)

string = r"ATTCCATTCCAATCCATTCCTTTCCTTTCGCTTGC"
matches = re.findall(string, seq)
map_score=1/len(matches)"""

#nice -19 python3 /home/projects/marthe/karyokan_ang/src/test.py /home/projects/marthe/karyokan_ang/data/hs37d5_100.genmap.bedgraph /home/projects/marthe/karyokan_ang/data/hs37d5_100.genmap.short

"""#Read the mappability file and calculate average mappability score for every window
referencefai = open("/home/databases/genomes/Homo_sapiens/hs37d5/hs37d5.fa.fai", "r") 
mapfile = open(sys.argv[1], "r")
mapfileshort =  open(sys.argv[2], "w")
#this code stores the start and end coordinates of the chromosomes found from samtools faidx 

winsize = 1000000
chrranges = []
chrnames   = {} 
for linefai in referencefai:
    faia=linefai.split("\t") 
    chrnames[ faia[0] ] = 0 
    for c in range(1, (int(faia[1])-winsize), winsize):
        rangedict =	{ "chr": faia[0], "start": c,   "end": c+(winsize-1), "ecov":1 }      
        chrranges.append(rangedict) 
referencefai.close() 


score=0
for rau in tqdm(range(len(chrranges))):
    start_fai = int(chrranges[rau]['start'])
    end_fai = int(chrranges[rau]['end'])
    for line in mapfile:
        if line.startswith("GL00") or line.startswith("hs37"):
            pass
        else:
            elements = line.split()
            start = int(elements[1])
            end = int(elements[2])
            chromosome = elements[0]
            if chromosome == chrranges[rau]['chr']:
                if start>start_fai and end<end_fai:
                    score+= 1/float(elements[3])
                else:
                    final_score = score/winsize
                    chrranges[rau]["map_score"] = final_score
                    mapfileshort.write(str(rau) + "\t"+ str(final_score) + "\n")
                    score=0
                    #print(line, final_score)
                    break
mapfileshort.close()
mapfile.close()"""

###################################################################################






    
