import sys
import gzip

'''
Subsets a file like oryLat04.${query}.all.chain.gz to just the chains that
are predicted to contain orthologous alignments.

'''

IDsFile = sys.argv[1]  # e.g. bolPec01.chain_ids
filename = sys.argv[2]  # e.g. oryLat04.bolPec01.all.chain.gz
outFileName= sys.argv[3] # e.g. oryLat04.bolPec01.ortho.chain

# chain_ids file looks like this:
# ENST00000000233 656 ENSG00000004059
# ...

IDs= set()
with open(IDsFile, 'r') as id_file:
	for line in id_file:
		chainID = line.strip().split(" ")[1]
		IDs.add(chainID)

outfile = open(outFileName, "w")

keepChain = False

with gzip.open(filename, 'rt') as chainFile:
	for line in chainFile:
		line = line.strip()
		if "#" in line:  # chain file header lines
			print(line, file=outfile)
			keepChain = False  # shouldn't be needed but just in case

		if line[0:5] == "chain":  # each chain begins wtih a header line containing the word "chain"
			scrap, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, \
			qStart, qEnd, chainID = line.split(" ")

			if chainID in IDs:
				keepChain = True
				print(line, file=outfile)
			else:
				keepChain = False
		elif keepChain:
			print(line, file=outfile)

outfile.close()
