#!/bin/bash -e

set -beEu -o pipefail

function usage {
    cat << EOF
Usage: script.sh oryLat04.{query}.ortho.chains.gz
EOF
    exit 1;
}

if [ $# -ne 1 ]; then
    usage;
fi




scriptsDir=$(realpath $(dirname $0))

chainFile=$(readlink -f $1)
target=$(basename $chainFile | cut -d"." -f1)
query=$(basename $chainFile | cut -d"." -f2)
originalGapTrackFile=${scriptsDir}/../processedInputs/gapTracks/${query}.gapTrack
out=${scriptsDir}/../processedInputs/DELs/
chromSizes=${scriptsDir}/../processedInputs/chromSizes/${query}.chrom.sizes

tempDirRoot="findChainDels"$(shuf -i1-10000 -n1)$(shuf -i1-10000 -n1)$(shuf -i1-10000 -n1)"_"$query

mkdir -p $out/$query 

HOMEDIR=/cluster/u/hchen17

if test -d "$HOMEDIR"; then
    mkdir -p /tmp/${tempDirRoot}
    intermediateDir=/tmp/${tempDirRoot}
else
    mkdir -p ${out}/${query}/${tempDirRoot}
    intermediateDir=${out}/${tempDirRoot}
fi

cp $(readlink -f $0) $out/$query/copyOfRunScript_$(basename $0)_$(date +%Y%m%d_%H%M)


slopBed -i $originalGapTrackFile -g $chromSizes -b 100 > ${intermediateDir}/${query}.100bpPaddedGapTrack

python3 $scriptsDir/findDels_no5bpPadding.py $chainFile ${intermediateDir}/$query".delsRelTo."$target".unfiltered.bed"

subtractBed -a ${intermediateDir}/$query".delsRelTo."$target".unfiltered.bed" -b  ${intermediateDir}/${query}.100bpPaddedGapTrack > ${intermediateDir}/$query".delsRelTo."$target".paddedGapSubtracted.bed"

python3 $scriptsDir/switchCoords.py ${intermediateDir}/$query".delsRelTo."$target".paddedGapSubtracted.bed" ${intermediateDir}/$query".delsRelTo."$target".paddedGapSubtractedCoordsSwitched.bed"

sort -k4,4 ${intermediateDir}/$query".delsRelTo."$target".paddedGapSubtractedCoordsSwitched.bed" > ${intermediateDir}/$query".delsRelTo."$target".paddedGapSubtractedSortedByChain.bed"

cut -f4 ${intermediateDir}/$query".delsRelTo."$target".paddedGapSubtractedSortedByChain.bed"  | awk -F"_" '{print $1 "_" $2 "_"}' | sort -u > ${intermediateDir}/listOfChainIDs # e.g. balAcu1_chain2049_

while read currChain; do
	grep $currChain ${intermediateDir}/$query".delsRelTo."$target".paddedGapSubtractedSortedByChain.bed" | sort -u | sort -k1,1 -k2,2n > ${intermediateDir}/$currChain".bed"
	mv ${intermediateDir}/$currChain".bed" $out/$query/$currChain".bed"
done < ${intermediateDir}/listOfChainIDs

rm -rf ${intermediateDir}
