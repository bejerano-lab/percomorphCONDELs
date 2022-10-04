
"""
Converts bed file outputs from findDels.py from 

chr1    191    234  chr3.897.1453.balAcu1_chain13902_SS

to

chr3	897	1453	balAcu1_chain13902_SS

"""
import sys

filename = sys.argv[1]
outfilename = sys.argv[2]

forbiddenChr = ["chrUn", "random", "chrM", "alt"]

outfile = open(outfilename, "w")

with open(filename, 'r') as delFile:
    for line in delFile:

        # only forbid weird chr in target (i.e. hg38)
        if any(x in line.strip().split("\t")[3] for x in forbiddenChr):
            continue

        line = line.strip()
        entry = line.split("\t")[3]
        print(entry.replace("*", "\t"), file=outfile)

outfile.close()
