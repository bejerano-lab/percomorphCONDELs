import sys, gzip
from pathlib import Path

script_dir = str(Path( __file__ ).parent.absolute())
outdir = script_dir + "/../processedInputs/filterBEDs/"


outstreams = {}

features = {"CDS": "oryLat04_ensembl98_codingExons.bed", 
            "nonCodingExon": "oryLat04_ensembl98_nonCodingExons.bed", 
            "five_prime_utr": "oryLat04_ensembl98_5UTR.bed", 
            "three_prime_utr": "oryLat04_ensembl98_3UTR.bed"}

for f in features:
    outstream = open(outdir+features[f], "w")
    outstreams[f] = outstream


with gzip.open(sys.argv[1], "rt") as gtf:

    for line in gtf:
        if "#"  in line:
            continue

        words = line.strip().split("\t")
        chrom = "chr"+words[0]
        start = int(words[3])-1 # convert to 0-based
        end = words[4]
        geneID = words[8].split("\"")[1] + "_" + words[8].split("\"")[5]
        featureType = words[2]

        printMe = False


        if "utr" in featureType:
            outfile = outstreams[featureType]
            printMe = True
        elif featureType == "CDS":
            outfile = outstreams[featureType]
            printMe = True
        elif (featureType == "exon") and ("protein_coding" not in line):
            outfile = outstreams["nonCodingExon"]
            printMe = True

        if printMe:
            print("\t".join([chrom, str(start), end, geneID]), file=outfile)


for stream in outstreams:
    outstreams[stream].close()
