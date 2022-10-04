#!/usr/bin/env python3

import sys
import gzip

'''
Takes a file from ensembl custom downloads like with 7 columns:
1	Gene stable ID
2       Transcript stable ID
3	Chromosome/scaffold name
4	Transcription start site (TSS)
5	Gene name
6	Transcript length (including UTRs and CDS)
7       Gene description

#Gene stable ID  Transcript stable ID    Chromosome/scaffold name        Transcription start site (TSS)   Gene name      Transcript length (including UTRs and CDS)      Gene description

'''
filename = sys.argv[1]

geneDict = {}  # key = ENSG id | value = [chr, TSSpos, name, transcriptLength]

with gzip.open(filename, 'rt') as ensgFile:
	for line in ensgFile:
		if "#" not in line:
			line = line.strip()
			ensg, enst, chrom, pos, name, length, desc = line.split("\t")
			if name == "":
				name = "[no name]"
			#if ("scaffold" not in chrom ) and ("contig" not in chrom):
			#	chrom = "chr" + chrom
			if ensg not in geneDict:
				geneDict[ensg] = [ensg, enst, chrom, pos, name, length, desc]
			else:
				if int(geneDict[ensg][5]) < int(length):  # if length of current transcript is longer, update TSS info
					geneDict[ensg] = [ensg, enst, chrom, pos, name, length, desc]

for val in geneDict.values():
	print("\t".join(val))
