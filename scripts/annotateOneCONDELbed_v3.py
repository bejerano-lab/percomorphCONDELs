import sys, os, re
from pybedtools import BedTool
from subprocess import run,PIPE
from pathlib import Path

resultsDir = os.path.abspath(sys.argv[1]) 
oneCONDELBedPath = os.path.abspath(sys.argv[2])
currCONDELbed = BedTool(oneCONDELBedPath)
condelID = os.path.basename(oneCONDELBedPath).split(".")[0]
folderNum = str(int(int(condelID.split("_")[1])/10))
outAnnotationName = condelID+".annot"
annotDirName = "annotationIntermediateFiles"
outAnnotationPath = "/".join([resultsDir,annotDirName,folderNum,outAnnotationName])

#print(outAnnotationPath)

script_dir = Path( __file__ ).parent.absolute()

filterBEDdirPath = str(script_dir)+"/../processedInputs/filterBEDs/"

defFile = resultsDir+"/DEF"

def grepDEF(keyword):
    line = run(["grep", keyword, defFile], stdout=PIPE).stdout.decode().strip()
    value = re.sub('\"', '', line.split("=")[1])
    return value

POSSIBLE_TARGETS = sorted(["molMol01", "takFla02", "takRub01", "tetNig2", "synSco01", "hipCom01", "hipEre01", "cynSem", "monAlb01"])
POSSIBLE_OUTGROUPS = sorted(["ampCit01", "oreNil02", "neoBri01", "hapBur01", "punNye01", "mayZeb03", "notFur02", "kryMar01", "ausLim01", 
                    "xipMac02", "xipHel01", "xipCou01", "poeRet02", "poeFor01", "funHet01", "latCal01", "parOli02", "serQui01", 
                    "serDum01", "priCar1", "gasAcu14", "punPun02", "larCro01", "miiMii01", "dicLab01", "labBer01", "bolPec01"])
TARGET_CLADES = sorted(["tetraodontiformes", "syngnathids", "sole", "eel"])
TARGET_CLADE_ASSIGNMENTS = {"molMol01": "tetraodontiformes", "takFla02": "tetraodontiformes", "takRub01": "tetraodontiformes", 
                            "tetNig2": "tetraodontiformes", "synSco01": "syngnathids", "hipCom01": "syngnathids",
                            "hipEre01": "syngnathids", "cynSem": "sole", "monAlb01": "eel"}

windowSizes = ["10", "25", "50", "100"]

ORTHO_CHAIN_DIR = str(script_dir)+"/../processedInputs/orthoChains/"

GENE_REGEX = 'ENSORLG0+'

# chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id 
CHAIN_HEADER_LINE_REGEX = "chain [0-9]\\+ \\S\\+ [0-9]\\+ [-+] [0-9]\\+ [0-9]\\+ \\S\\+ [0-9]\\+ [-+] [0-9]\\+ [0-9]\\+ "

# return a bed interval in reference coordinates spanned by the chain
def getChainCoverageBed(chainID, chainFilePath):
    chainRegEx = CHAIN_HEADER_LINE_REGEX + chainID
    chainHeader = run(["zgrep", "-we", chainRegEx, chainFilePath], stdout=PIPE).stdout.decode().strip().split(" ")
    # print(chainHeader)
    bedInterval = " ".join([chainHeader[2], chainHeader[5], chainHeader[6]])
    return BedTool(bedInterval, from_string=True)

## COLLECT SCREEN PARAMETERS ##

actualTargets = []
for clade in TARGET_CLADES:
    species = grepDEF(clade).split(" ")
    if species[0] != "":
        actualTargets.extend(species)
# print(actualTargets)

lookupDirPath = grepDEF("lookupDir")
# print(lookupDir)


delDirPath = grepDEF("delsDir")
consWindowsDirPath = grepDEF("consWindowsDir")
intactAxtsDirPath = grepDEF("intactAxtsDir")
minConSize = int(grepDEF("minConSize"))
mergeDist = int(grepDEF("mergeDist"))
minCONDELsize = int(grepDEF("minConDelSize"))

geneTaggedCONDELbed = BedTool(resultsDir+"/oryLat04.CONDELs.bed")

threePrimeUTR = BedTool(filterBEDdirPath+"/oryLat04_ensembl98_3UTR.bed")
fivePrimeUTR = BedTool(filterBEDdirPath+"/oryLat04_ensembl98_5UTR.bed")
codingExons = BedTool(filterBEDdirPath+"/oryLat04_ensembl98_codingExons.bed")
introns = BedTool(filterBEDdirPath+"/oryLat04_ensembl98_introns.bed")
nonCodingExons = BedTool(filterBEDdirPath+"/oryLat04_ensembl98_nonCodingExons.bed")

filterBeds = [codingExons, fivePrimeUTR, threePrimeUTR, introns, nonCodingExons]

canonicalTssInfo = BedTool(filterBEDdirPath+"/ASM223467v1_ensembl98_canonicalTSS_withInfo.tab")


#all_SS_targetDels_bed = BedTool("/cluster/u/hchen17/treeWAS/inputData/percomorphDels/TARGETs_20200803_SSgaps")

pDels = BedTool(resultsDir+"/oryLat04.numberedDELs.bed")

#rand = uuid.uuid1() # uniq id to avoid clashes (though not likely)


def getIntersectingElementNames(filterBed, cBed):
    result = "[none]"
    intersectingElements = set()

    for i in filterBed.intersect(cBed):
        intersectingElements.add(i.name)

    if len(intersectingElements) > 0:
        result = ",".join(sorted(list(intersectingElements)))

    return result

counter = 1
with open(outAnnotationPath, "w") as resultFile:
    # get list of orthologs/chains that were used to call this candidate
    geneTags = set()
    for i in geneTaggedCONDELbed.intersect(currCONDELbed):
        geneTags.add(i.name.strip().split("_")[1])
    genes = ",".join(sorted(list(geneTags)))

    # get list of window regimen(s) used to call this candidate
    windowsUsed = []
    windowTaggedCONDELbeds = []
    for x in range(4):
        if len(BedTool(resultsDir+"/uncleanedResultFiles/oryLat04.CONDELsCombinedUnmerged_"+windowSizes[x]+"bp").intersect(currCONDELbed)) > 0:
            windowsUsed.append("yes")
        else:
            windowsUsed.append("no")

    # get any SS target dels involved in this candidate
    #targetSSdels = []
    #for i in all_SS_targetDels_bed.intersect(candidate):
        #targetSSdels.append("|".join(list(i)))
    #ssDels = "-"
    #if len(targetSSdels) > 0:
        #ssDels = ",".join(targetSSdels)
    
    # get the larger del in which this CONDEL is located
    pDel = pDels.intersect(currCONDELbed)[0].name # should only return 1 interval

    # get gene feature annotations
    ixElements = []
    for bed in filterBeds:
        ixElements.append(getIntersectingElementNames(bed, currCONDELbed))

    firstUsTSSinfo = "\t".join(currCONDELbed.closest(canonicalTssInfo, t="first", D="ref", fu=True)[0][7:])
    firstDsTSSinfo = "\t".join(currCONDELbed.closest(canonicalTssInfo, t="first", D="ref", fd=True)[0][7:])

    # for each gene used to map a given CONDEL
    for gene in list(geneTags):

        # run through each mapped species at this region to note violators
        # use this dict to note violators
        # set all to violator initially, then flip to "ok" as appropriate
        targetStatus = {}
        targetCladeStatus = {}
        outgroupStatus = {}
        for species in POSSIBLE_TARGETS:
            targetStatus[species] = "unmapped"
        for species in POSSIBLE_OUTGROUPS:
            outgroupStatus[species] = "unmapped"
        for clade in TARGET_CLADES:
            targetCladeStatus[clade] = "unmapped"

        numMappedTargets = 0
        numMappedOutgroups = 0

        folder = str(int(int(re.sub(GENE_REGEX,'', gene))/1000))
        lookupFilePath = "/".join([lookupDirPath, folder, gene])
        
        # get chain mappings for this gene
        with open(lookupFilePath, "r") as file:
            for line in file:
                if "\t" in line: # skip the gene TSS line: "chrX\t[start]\t[end]"
                    continue
                line = line.strip()
                species, chain = line.split("_") # e.g. larCro01_chain2
                #print(species, chain)

                chainFilePath = ORTHO_CHAIN_DIR + "/" + ".".join(["oryLat04", species, "ortho.chains.gz"])
                # print(chainFilePath)
                chainID = chain[5:] # get string representation of chain ID -- i.e. '98' from 'chain98'

                chainCoverageInterval = getChainCoverageBed(chainID, chainFilePath)

                chainBedName = "_".join([species,chain,".bed"]) # e.g. larCro01_chain2_.bed

                currChainDelsPath = "/".join([delDirPath, species, chainBedName])

                if species in POSSIBLE_TARGETS:
                    numMappedTargets+=1

                    if len(currCONDELbed.intersect(chainCoverageInterval)) == 0:
                        targetStatus[species] = "inconclusive - chain too short"
                        if targetCladeStatus[TARGET_CLADE_ASSIGNMENTS[species]] == "unmapped":
                            targetCladeStatus[TARGET_CLADE_ASSIGNMENTS[species]] = "inconclusive - chain too short"
                        # print("here")
                        continue


                    if os.path.exists(currChainDelsPath):
                        delBed = BedTool(currChainDelsPath).merge(d=mergeDist)
                        ixDels = delBed.intersect(currCONDELbed, wa=True)
                        if len(ixDels) > 0:
                            for i in ixDels:
                                if i.length >= minCONDELsize:
                                    targetStatus[species] = "ok"
                                    targetCladeStatus[TARGET_CLADE_ASSIGNMENTS[species]] = "ok"
                                else:
                                    targetStatus[species] = "small del"
                        else:
                            targetStatus[species] = "violator"
                            # as long as one species in the target clade is ok, call the whole clade ok
                            if targetCladeStatus[TARGET_CLADE_ASSIGNMENTS[species]] == "unmapped":
                                targetCladeStatus[TARGET_CLADE_ASSIGNMENTS[species]] = "violator"
                elif species in POSSIBLE_OUTGROUPS:
                    numMappedOutgroups+=1

                    if len(currCONDELbed.intersect(chainCoverageInterval)) == 0:
                        outgroupStatus[species] = "inconclusive - chain too short"
                        # print("here")
                        continue



                    exitEarly = False
                    for x in range(4):
                        # each CONDEL could've been mapped using any of 4 sliding window regimens
                        # if the CONDEL intersects any CONs element (regardless of which originating 
                        # window size regime), mark this species as a non-violator/valid
                        wSize = windowSizes[x]
                        consBedPath = "/".join([consWindowsDirPath, species, wSize,chainBedName])
                        if os.path.exists(consBedPath):
                            consBed = BedTool(consBedPath).merge(d=mergeDist)
                            ixConIntervals = consBed.intersect(currCONDELbed, wa=True)
                            #print(ixConIntervals)
                            if len(ixConIntervals) > 0:
                                for i in ixConIntervals:
                                    #print(i, i.length)
                                    if i.length >= minConSize:
                                        if windowsUsed[x] == "yes":
                                            outgroupStatus[species] = "conserved"
                                            exitEarly = True
                                            break # as long as one regime is officially ok, move on
                                        else:
                                            outgroupStatus[species] = "conserved but not via officially used window regime"
                                    else:
                                        if outgroupStatus[species] == "unmapped":
                                            outgroupStatus[species] = "cons interval smaller than minConSize"
                        if exitEarly:
                            break

                    if exitEarly:
                        continue

                    if outgroupStatus[species][0:4] == "cons":
                        continue

                    currAxtBedPath = "/".join([intactAxtsDirPath, species, chainBedName])
                    if os.path.exists(currAxtBedPath):
                        axtBed = BedTool(currAxtBedPath).merge(d=mergeDist)
                        ixAxts = axtBed.intersect(currCONDELbed) # just get intersecting portion axt interval (not entire axt block)
                        if len(ixAxts) > 0:
                            for i in ixAxts:
                                if i.length >= minCONDELsize:
                                    if "con" not in outgroupStatus[species]:
                                        outgroupStatus[species] = "intact axt exists but not highly conserved"
                                    exitEarly = True
                                    break
                                else:
                                    outgroupStatus[species] = "violator"

                    if exitEarly:
                        continue

                    if os.path.exists(currChainDelsPath):
                        delBed = BedTool(currChainDelsPath).merge(d=mergeDist)
                        candidateRemainingAfterIxDels = currCONDELbed.intersect(delBed, v=True)
                        if len(candidateRemainingAfterIxDels) > 0:
                            for i in candidateRemainingAfterIxDels:
                                if i.length < minCONDELsize:
                                    outgroupStatus[species] = "violator"
                                else:
                                    outgroupStatus[species] = "minimally conserved"
                                    break
                    if outgroupStatus[species] == "unmapped":
                        outgroupStatus[species] = "violator"

                else:
                    print("invalid species encountered")
                    print("invalid species encountered", file=resultFile)
                    exit(1)

        numTargetViolators = len([i for i in targetStatus.values() if i == "violator"])
        numTargetOk = len([i for i in targetStatus.values() if i == "ok"])

        numOutgroupViolators =  len([i for i in outgroupStatus.values() if i == "violator"])
        numOutgroupOk = len([i for i in outgroupStatus.values() if i[0:9] == "conserved"])


        numTargetCladeViolators = len([i for i in targetCladeStatus.values() if i == "violator"])
        numMappedSpecies = numMappedOutgroups+numMappedTargets

        targetStatusList = "\t".join([targetStatus[x] for x in sorted(targetStatus.keys())])
        ogStatusList = "\t".join([outgroupStatus[x] for x in sorted(outgroupStatus.keys())])
        cladeStatusList = "\t".join([targetCladeStatus[x] for x in sorted(targetCladeStatus.keys())])

        print(str(currCONDELbed).strip("\n"), condelID, gene, ",".join(windowsUsed), 
            firstUsTSSinfo, firstDsTSSinfo, pDel, "\t".join(ixElements), numMappedSpecies, numTargetViolators, numTargetCladeViolators,
            numTargetOk, numMappedTargets, numOutgroupViolators, numOutgroupOk, numMappedOutgroups, targetStatusList, ogStatusList, 
            cladeStatusList, sep="\t", file=resultFile)
