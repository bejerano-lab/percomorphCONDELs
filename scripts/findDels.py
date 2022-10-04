

"""

SUMMARY: Given a gzipped chain file that encodes a pairwise alignment between target assembly
and query assembly, this script prints out a BED entry for each single- or double-sided gap 
present in all chains within the gzipped file, written in query coordinates and tagged with
a representation of reference coordinates.

"""

import sys
import gzip

filename = sys.argv[1]
outfilename = sys.argv[2]

dels_queryCoords = []  # query deletions given in target or query coordinates
delSizeThreshold = 10  # min deletion size to consider
delPaddingThreshold = 0  # padding threshold (to make sure deletions are "nowhere near" assembly gaps

currSpecies = filename.split("/")[len(filename.split("/"))-1].split(".")[1]  # e.g. oviAri3

with gzip.open(filename, 'rt') as chainFile:

    for line in chainFile:

        if "#" in line:  # ignore chain file header lines
            continue

        line = line.strip()

        if "chain" in line:  # each chain begins wtih a header line containing the word "chain"
            scrap, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, \
            qStart, qEnd, chainID = line.split(" ")

            if qStrand == "-":  # change arithmetic method for chains encoding inversions
                qCurrPos = int(qSize) - int(qStart)
            else:
                qCurrPos = int(qStart)
            tCurrPos = int(tStart)

            # tCurrPos is current position in target genome coordinates
            # qCurrPos is current position in query genome coordinates

            # dels_queryCoords.append("# Dels_qCoord: " + line)

        elif len(line.split()) > 1:
            alignmentSize, dt, dq = line.split("\t")

            if qStrand == "-":
                #qEntry = "\t".join([qName, str(qCurrPos - int(alignmentSize)), str(qCurrPos)])
                qCurrPos -= (int(alignmentSize) + int(dq))
            else:
                #qEntry = "\t".join([qName, str(qCurrPos), str(qCurrPos + int(alignmentSize))])
                qCurrPos += (int(alignmentSize) + int(dq))

            if int(dt) >= delSizeThreshold:

                kind = "SS"  # single-sided chain gap

                if int(dq) > 0:
                    kind = "DS"  # double-sided chain gap

                delID = "_".join([currSpecies+"_chain"+chainID, kind])  # BED entry ID that looks like e.g. "balAcu1_chain13902_SS"

                # save off the deletion in terms of target genome coordinates
                tDelEntry = "*".join([tName, str(tCurrPos + int(alignmentSize)), str(tCurrPos + int(alignmentSize) + int(dt)), delID])

                # save off the deletion in terms of query coordinates
                if qStrand == "-":
                    qDelEntry = "\t".join([qName, str(qCurrPos - delPaddingThreshold), str(qCurrPos + int(dq) + delPaddingThreshold)])
                else:
                    qDelEntry = "\t".join([qName, str(qCurrPos - int(dq) - delPaddingThreshold), str(qCurrPos + delPaddingThreshold)])

                # add info to list of deletions
                dels_queryCoords.append(qDelEntry + "\t" + tDelEntry)

            # update target genome position
            tCurrPos += (int(alignmentSize) + int(dt))

""" Example output line: chr1    191    234  chr3*897*1453*balAcu1_chain13902_SS """
with open(outfilename, "w") as outfile:
    for element in dels_queryCoords:
        print(element, file=outfile)

