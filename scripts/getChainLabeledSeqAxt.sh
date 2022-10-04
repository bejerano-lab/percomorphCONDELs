#!/bin/bash -e

set -beEu -o pipefail

function usage {
	cat << EOF
Usage: script.sh orthoChainIDMappings

e.g. pickOrthoChains/molMol01.chain_ids 
ENSORLG00000000002 456 ENSORLT00000000002
...

EOF
	exit 1;
}

if [ $# -ne 1 ]; then
	usage;
fi

query=$(basename $1 .chain_ids)

scriptsDir=$(realpath $(dirname $0))

orthoChains=${scriptsDir}/../processedInputs/orthoChains/oryLat04.${query}.ortho.chains.gz
mappings=$(readlink -f $1)

target2Bit=${scriptsDir}/../processedInputs/2bits/oryLat04.2bit
query2Bit=${scriptsDir}/../processedInputs/2bits/${query}.2bit
out=${scriptsDir}/../processedInputs/labeledAxts/

root=$(basename ${orthoChains} chains.gz)"chainLabeled.axts"


mkdir -p ${out}/${root}
cp $(readlink -f $0) $out/$root/copyOfRunScript_$(basename $0)


for chainID in $(cut -d" " -f2 ${mappings} | sort -u); do
	zcat ${orthoChains} | chainFilter -id=${chainID} stdin | chainToAxt stdin ${target2Bit} ${query2Bit} stdout | awk -v descriptor=$query"_chain"$chainID"_" '{if ($1 ~ /[0-9]/) $0 = $0 " " descriptor; print $0 }' > ${out}/${root}/$query"_chain"$chainID"_.axt"
	gzip ${out}/${root}/$query"_chain"$chainID"_.axt"
done
