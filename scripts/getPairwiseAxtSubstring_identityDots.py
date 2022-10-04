import sys
import re

"""
Given an axt file, return the nucleotide axt corresponding to the specified interval coordinates.

Usage: getPairwiseSubstring.py nucleotideAlignment.axt bedStartPos bedEndPosNotInclusive outfilename

e.g. getPairwiseSubstring.py nucleotideAlignment.axt 52353197 52353207 would output...

	AGATT----TTAGA
	AGATTATATTTGGT

...to outfilename, where nucleotideAlignment.axt contained...

	40 chr14 52353197 52353297 chr14 45052736 45052861 + 4143
	AGATT----TTAGATTATAAATTATTTTTAAAAATCAAACATCCTGGGGAAATGTTTTTAATACATTAATAATATTATTGTCAAGGTCTTTTCCCACAAAACCAC
	AGATTATATTTGGTTTATAAACTGTTTTTTAAAATCATGCATCCTCTAGAAAATATTTTAGCACAT---TGACATTATTATCAAGATGTTTCCCTGCAGAAGCAC

"""



# def print_axt_using_identity_dots(refSeq, querySeq, outfile):
	# transformedQuerySeq=""
	# for i in range(0, len(refSeq)):
	# 	if refSeq[i] is querySeq[i]:
	# 		transformedQuerySeq += "."
	# 	else:
	# 		transformedQuerySeq += querySeq[i]
# 	outfile.write(refSeq)
# 	outfile.write(transformedQuerySeq)
# 	outfile.close()
# 	return

axtFilename = sys.argv[1]

intervalStart = int(sys.argv[2])
intervalEnd = int(sys.argv[3])  # noninclusive

outfile = open(sys.argv[4], "w")

axtStringStartIndex = -1
axtStringEndIndex = -1

with open(axtFilename, 'rt') as axt:
	for line in axt:
		#print(line)
		if re.search(r"[0-9]+ [a-zA-Z]+", line): # 0 CP020681.1 6103595 6103684 CCOE01000292.1 1600015 1600104 - 6036
			chrom = line.strip().split(" ")[1]
			axtStart = int(line.strip().split(" ")[2]) - 1  # switch from 1-based to 0-based coords
			axtEnd = int(line.strip().split(" ")[3]) + 1  # switch from 1-based to 0-based coords

			if (axtStart <= intervalStart) and (axtEnd >= intervalEnd):  # correct axt block has been reached
				refSeq = axt.readline().strip().upper()
				querySeq = axt.readline().strip().upper()
				#print(refSeq)
				#print(querySeq)
				genomePos = int(axtStart)

				for i in range(0,len(refSeq)):

					if (genomePos == intervalStart) and (refSeq[i] != "-"):
						axtStringStartIndex = i

					if (genomePos == intervalEnd-1) and (refSeq[i] != "-"):  # noninclusive
						axtStringEndIndex = i

					if (axtStringStartIndex >= 0) and (axtStringEndIndex >= 0):
						
						refSubSeq = refSeq[axtStringStartIndex:axtStringEndIndex+1]
						querySubSeq = querySeq[axtStringStartIndex:axtStringEndIndex+1]
						
						transformedQuerySeq = ""

						for i in range(0, len(refSubSeq)):
							if refSubSeq[i] is querySubSeq[i]:
								transformedQuerySeq += "."
							else:
								transformedQuerySeq += querySubSeq[i]
						
						print(refSubSeq, file=outfile)
						print(querySubSeq, file=outfile)
						print(transformedQuerySeq, file=outfile)

						sys.exit()

					if refSeq[i] != "-":
						genomePos += 1

			else:
				axt.readline()
				axt.readline()

