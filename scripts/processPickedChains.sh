#processPickedChains.sh 

#!/bin/bash

set -beuEo pipefail

species=$(basename $1 .out)
scriptsDir=$(dirname $0)
setupOutDir=${scriptsDir}/../pickOrthoChains/
orthoChainOutDir=${scriptsDir}/../processedInputs/orthoChains/
canonicalIDs=${scriptsDir}/../pickOrthoChains/ensembl98_ASM223467v1_canonical_IDs.txt

# REFORMAT PICK CHAIN OUTPUT
join <(grep ^ENSORLG $1 | cut -d" " -f1,2 | sort) ${canonicalIDs} > ${setupOutDir}/${species}".chain_ids"

awk -v var=$species '{print $1 "\t" var"_chain"$2}'  ${setupOutDir}/${species}".chain_ids" > ${setupOutDir}/${species}_reformattedChainIds


# SUBSET TO ORTHO CHAINS
${scriptsDir}/subsetToOrthoChain.sh ${setupOutDir}/${species}".chain_ids" ${setupOutDir}/chains/oryLat04.${species}.all.chain.gz ${orthoChainOutDir}
