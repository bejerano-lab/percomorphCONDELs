"""

SUMMARY: Given a gzipped chain file that encodes a pairwise alignment between target assembly
and query assembly, this script prints out a BED entry for each region of contiguous alignment 
from the query to the target genome. Output is in target genome coordinates.

"""

import sys
import gzip

filename = sys.argv[1]
outfilename = sys.argv[2]
intactRegions = []  # query axt blocks given in target or query coordinates
currSpecies = filename.split("/")[len(filename.split("/"))-1].split(".")[1]  # e.g. oviAri3


outfile = open(outfilename, 'w')

with gzip.open(filename, 'rt') as chainFile:
    for line in chainFile:

        if "#" in line:  # ignore chain file header lines
            continue

        line = line.strip()

        if "chain" in line:  # each chain begins with a header line containing the word "chain"
            scrap, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, \
            qStart, qEnd, chainID = line.split(" ")

            tCurrPos = int(tStart)

        elif len(line.split()) > 1:
            alignmentSize, dt, dq = line.split("\t")

            # save off the aligning region in terms of target genome coordinates
            print("\t".join([tName, str(tCurrPos), str(tCurrPos + int(alignmentSize)), currSpecies + "_chain" + chainID + "_"]), file=outfile)

            # update target genome position
            tCurrPos += (int(alignmentSize) + int(dt))

        elif line != "":
            print("\t".join([tName, str(tCurrPos), str(tCurrPos + int(line)), currSpecies + "_chain" + chainID + "_"]), file=outfile)

outfile.close()
