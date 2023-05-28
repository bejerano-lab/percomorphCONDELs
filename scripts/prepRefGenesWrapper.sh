#!/bin/bash

set -beEuo pipefail

scriptsDir=$(realpath $(dirname $0))

ensemblBiomartFile=$(readlink -f $1)
ensemblGTFgz=$(readlink -f $2)

python3 ${scriptsDir}/prepRefGeneInfo.py ${ensemblBiomartFile}

sort -k1,1 -k2,2n $scriptsDir"/../processedInputs/filterBEDs/ASM223467v1_ensembl98_canonicalTSS_withInfo.tab" > /tmp/tempSortFile
mv /tmp/tempSortFile $scriptsDir"/../processedInputs/filterBEDs/ASM223467v1_ensembl98_canonicalTSS_withInfo.tab"


## identify introns ##

# first make a gff3 file of exons from ensembl's gtf

zcat ${ensemblGTFgz} | awk '$3 == "exon"' > /tmp/exonsOnly.gtf
gt gtf_to_gff3 -force -o /tmp/exonsOnly.gff3 /tmp/exonsOnly.gtf
gt gff3 -addintrons -sort /tmp/exonsOnly.gff3 | sed 's/ /_/g' > /tmp/exonsAndIntrons.gff3 

python3 ${scriptsDir}/parseExonIntronGff3.py

rm /tmp/exons*

python3 ${scriptsDir}/parseReferenceGTF.py ${ensemblGTFgz}

