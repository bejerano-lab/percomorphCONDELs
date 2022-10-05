from pathlib import Path

script_dir = str(Path( __file__ ).parent.absolute())
intronsBedPath = script_dir + "/../processedInputs/filterBEDs/oryLat04_ensembl98_introns.bed"

parentIDdict = {}
intronsDict = {}


def makeAttributeDict(attributeList):
    d = {}
    for i in attributeList:
        key, val = i.split("=")
        d[key] = val
    return d

with open("/tmp/exonsAndIntrons.gff3", "r") as gff3:
    for line in gff3:
        if "#"  in line:
            continue
        words = line.strip("\n").split()
        
        if words[2] == "mRNA":
            attributes = makeAttributeDict(words[8].split(";"))
            mRNA_id = attributes["ID"]
            transcript_id = attributes["transcript_id"]
            gene_id = attributes["gene_id"]
            parentIDdict[mRNA_id] = gene_id + "_" + transcript_id

        if words[2] == "intron":
            attributes = makeAttributeDict(words[8].split(";"))
            chrom = words[0]
            start = str(int(words[3]) - 1)
            end = words[4]
            mRNA_id = attributes["Parent"]
            intronList = intronsDict.get(mRNA_id, [])
            intronList.append([chrom, start, end])
            intronsDict[mRNA_id] = intronList

with open(intronsBedPath, "w") as intronsBed:

    for mRNA_id in intronsDict:
        for intron in intronsDict[mRNA_id]:
            ensemblId = parentIDdict[mRNA_id]
            bedInfo = intron
            bedInfo.append(ensemblId)
            print("\t".join(bedInfo), file=intronsBed)
