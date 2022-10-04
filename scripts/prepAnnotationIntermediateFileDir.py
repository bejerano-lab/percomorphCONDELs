import os, sys

resultsDir = sys.argv[1]
numCONDELs = int(sys.argv[2])

numFoldersToMake = int(numCONDELs / 10)

os.mkdir("/".join([resultsDir, "annotationIntermediateFiles"]))

for i in range(numFoldersToMake+1):
    os.mkdir("/".join([resultsDir, "annotationIntermediateFiles", str(i)]))

