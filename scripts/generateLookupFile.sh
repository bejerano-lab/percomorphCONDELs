#!/bin/bash

# This script generates an ortho chain "lookup file" for each canonical gene and any aligning species to which the gene has been mapped

# referenceGeneTSSs:
# chr1	58832	58833	[no name]	ENSORLG00000025847	ENSORLT00000035921	[no description]
# ...

set -beEuo pipefail

scriptsDir=$(realpath $(dirname $0))

gene=$1 # e.g. ENSORLG00000025847

num=$(echo $gene | sed 's/ENSORLG//' | sed "s/^0\+//")
folderPrefix=$(($num/1000))

#echo $num $folderPrefix

outroot=${scriptsDir}/../processedInputs/orthoChainLookupFiles
inroot=${scriptsDir}/../pickOrthoChains/

mkdir -p ${outroot}/${folderPrefix}

outpath=${outroot}/${folderPrefix}/$gene

referenceGeneTSSs=${scriptsDir}/../processedInputs/filterBEDs/ASM223467v1_ensembl98_canonicalTSS_withInfo.tab

grep $gene $referenceGeneTSSs | awk -F"\t" 'OFS="\t" {print $1, $2, $3}' > $outpath
results=$(grep "$gene" ${inroot}/*reformattedChainIds)

if [ -n "${results}" ]
then 
	echo "$results" | cut -d" " -f2 | sort -u >> $outpath
fi

