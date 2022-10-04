#!/bin/bash -e

set -beEu -o pipefail

function usage {
    cat << EOF
Usage: script.sh oryLat04.${query}.ortho.chains.gz 
EOF
    exit 1;
}

if [ $# -ne 1 ]; then
    usage;
fi

scriptsDir=$(realpath $(dirname $0))
chainFile=$(readlink -f $1)
out=${scriptsDir}/../processedInputs/axtBlocks

target=$(basename $chainFile | cut -d"." -f1)
query=$(basename $chainFile | cut -d"." -f2)
tempDirRoot="findIntactAxts"$(shuf -i1-10000 -n1)$(shuf -i1-10000 -n1)$(shuf -i1-10000 -n1)"_"$query

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

python3 ${scriptsDir}/findIntactRegions.py $chainFile ${intermediateDir}/$query".IntactAligmntsRelTo."$target"_unsorted.bed"

sort -k4,4 ${intermediateDir}/$query".IntactAligmntsRelTo."$target"_unsorted.bed" > ${intermediateDir}/$query".IntactAligmntsRelTo."$target"_sortedByChainID.bed"

# split output by chainID *
while [ -s ${intermediateDir}/$query".IntactAligmntsRelTo."$target"_sortedByChainID.bed" ]; do  # while the file still has lines
	currChain=$(head -1 ${intermediateDir}/$query".IntactAligmntsRelTo."$target"_sortedByChainID.bed" | cut -f4 | awk -F"_" '{print $1 "_" $2 "_"}')  # e.g. balAcu1_chain2049_

	# save off all entries from the current chain
	cat ${intermediateDir}/$query".IntactAligmntsRelTo."$target"_sortedByChainID.bed" | grep $currChain | sort -k1,1 -k2,2n > ${intermediateDir}/$currChain".bed"

	# count how many entries came from current chain
	numLines=$(wc -l ${intermediateDir}/$currChain".bed" | cut -d" " -f1)

	# move saved entries from current chain to outdir
	cat ${intermediateDir}/$currChain".bed" > $out/$query/$currChain".bed"
	rm ${intermediateDir}/$currChain".bed"

	# delete entries that have already been considered
	sed -i '1,'"$numLines"'d' ${intermediateDir}/$query".IntactAligmntsRelTo."$target"_sortedByChainID.bed"
done

rm -r ${intermediateDir}
