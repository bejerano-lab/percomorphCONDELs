#!/bin/bash

# This script generates an ortho chain "lookup file" for each canonical gene and any aligning species to which the gene has been mapped

# referenceGeneTSSs:
# chr1	1604082	1604083	ENSORLG00000000001
# ...

set -beEuo pipefail

scriptsDir=$(realpath $(dirname $0))

gene=$1 # e.g. ENSORLG00000000001

num=$(echo $gene | sed 's/ENSORLG//' | sed "s/^0\+//")
folderPrefix=$(($num/1000))

#echo $num $folderPrefix

outroot=${scriptsDir}/../processedInputs/orthoChainLookupFiles
inroot=${scriptsDir}/../pickOrthoChains/

mkdir -p ${outroot}/${folderPrefix}

outpath=${outroot}/${folderPrefix}/$gene

referenceGeneTSSs=${scriptsDir}/../processedInputs/filterBEDs/ASM223467v1_ensembl98_canonicalTSS.bed

grep $gene $referenceGeneTSSs | awk '{print $1 "\t" $2 "\t" $3}' > $outpath
results=$(grep "$gene" ${inroot}/*reformattedChainIds)

if [ -n "${results}" ]
then 
	echo "$results" | cut -d" " -f2 | sort -u >> $outpath
fi

