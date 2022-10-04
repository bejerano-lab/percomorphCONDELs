#!/bin/bash

set -beuEo pipefail

species=$1
scriptsDir=$(dirname $0)
setupOutDir=${scriptsDir}/../pickOrthoChains/
orthoChainOutDir=${scriptsDir}/../processedInputs/orthoChains/
canonicalIDs=${scriptsDir}/../processedInputs/filterBEDs/ASM223467v1_ensembl98_canonicalTSS_withInfo.tab

# REFORMAT PICK CHAIN OUTPUT
join <(grep ^ENSORLG ${setupOutDir}/${1}.out | cut -d" " -f1,2 | sort) <(cut -f5-6 ${canonicalIDs} | sed 's/\t/ /' | sort) > ${setupOutDir}/${species}".chain_ids"

awk -v var=$species '{print $1 "\t" var"_chain"$2}'  ${setupOutDir}/${species}".chain_ids" > ${setupOutDir}/${species}_reformattedChainIds


# SUBSET TO ORTHO CHAINS
${scriptsDir}/subsetToOrthoChain.sh ${setupOutDir}/${species}".chain_ids" ${setupOutDir}/chains/oryLat04.${species}.all.chain.gz ${orthoChainOutDir}
