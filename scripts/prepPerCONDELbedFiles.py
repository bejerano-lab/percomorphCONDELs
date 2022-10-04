import sys

resultsDir = sys.argv[1]
uniqCONDELsbedPath = resultsDir + "/" + "oryLat04.CONDELs.lengthsOnly.bed"
intermediateOutdir = resultsDir + "/" + "annotationIntermediateFiles"
bedPrefix = "pCONDEL_"

counter = 1
with open(uniqCONDELsbedPath, 'r') as file:
    for line in file:
        line = line.strip()
        if '#' not in line:
            subfolder = str(int(counter/10))
            outfilename = "/".join([intermediateOutdir, subfolder, bedPrefix+str(counter)+".bed"])
            with open(outfilename, 'w') as outfile:
                print(line, file=outfile)
            counter+=1
