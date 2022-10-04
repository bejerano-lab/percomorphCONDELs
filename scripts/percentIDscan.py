
import sys
import gzip

"""
Compute percent ID for each window of size STEP_SIZE across chain-derived alignment (axt) files 
and output bed entry for each windw with the chain/ID it came from and its % ID calculation, if
the %ID is above some minimum threshold to be printed (otherwise output file is ginormous)

Output format:

chr   start   end   species_chain#_   %ID

e.g.

chr14	24687985	24688085	mm10_chain1_	0.609523809524

"""
filename = sys.argv[1]  # /cluster/u/hchen17/aquaticMammals/findDels/OUTGROUP_CHAINS_2BITS/hg38.bosTau8.ortho.axt.gz
STEP_SIZE = int(sys.argv[2])  # step size will be implemented in terms of the ungapped HUMAN genomic interval!!
PCT_ID_THRESHOLD_FOR_PRINTING = float(sys.argv[3])  # to save space, only print out entries/windows that have % ID greater than threshold
outfilename = sys.argv[4]

outfile = open(outfilename, "w")

with gzip.open(filename, 'rt') as axt:
	for line in axt:
		if "chain" in line:  # 0 CP020681.1 6103595 6103684 CCOE01000292.1 1600015 1600104 - 6036 ampCit01_chain1000_

			chrom, start, end = line.strip().split(" ")[1:4]
			chainID=line.strip().split(" ")[len(line.strip().split(" "))-1]  # e.g. mm10_chain100_
			hg38 = axt.readline().strip().upper()
			query = axt.readline().strip().upper()

			genomePos = int(start) - 1 # switch from 1-based to 0-based coords
			genomeAxtLength = int(end) - int(start) + 1

			gappedAxtPos = 0
			gappedAxtLength = len(hg38)

			while genomePos <= int(end) - STEP_SIZE:
				currHg38sequence = ""
				effectiveStepCounter = 0
				gappedStepCounter = 0

				while (effectiveStepCounter < STEP_SIZE):  # get effective window sequence
					currHg38Base = hg38[gappedAxtPos + gappedStepCounter]
					currHg38sequence += currHg38Base
					gappedStepCounter += 1
					if (gappedAxtPos + gappedStepCounter) == gappedAxtLength:  # hit end of chain without making entire chain
						gappedStepCounter -= 1
						break

					if currHg38Base != "-":  # don't count hg38 gaps [---] in the step count so that the step size accurately reflects x bp of the ungapped human genomic interval
						effectiveStepCounter += 1

				querySeq = ""

				if ("N" not in currHg38sequence) and (currHg38sequence[0] != "-") and (len(currHg38sequence) == STEP_SIZE):  # don't consider windows that are in hard-masked regions (i.e. contain Ns )
																				    # and don't use windows where hg38 sequence starts with "-"
					querySeq = query[gappedAxtPos:gappedAxtPos+gappedStepCounter]

					matches = sum(nt1 == nt2 for nt1, nt2 in zip(currHg38sequence, querySeq))
					pct_identity = 1.0 * matches / len(currHg38sequence)

					if pct_identity >= PCT_ID_THRESHOLD_FOR_PRINTING:
						bedEntry = "\t".join([chrom, str(genomePos), str(genomePos + effectiveStepCounter), chainID, str(pct_identity)])
						print(bedEntry, file=outfile)

				gappedAxtPos += 1
				if hg38[gappedAxtPos] != "-":
					genomePos += 1

outfile.close()