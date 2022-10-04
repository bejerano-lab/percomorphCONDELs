#!/bin/bash

set -beEuo pipefail
scriptsDir=$(realpath $(dirname $0))
mappings=$(readlink -f $1)
allChains=$(readlink -f $2)
outdir=$(readlink -f $3)
target=$(basename $allChains | cut -d"." -f1)
query=$(basename $allChains | cut -d"." -f2)
orthoChainName=${target}.${query}.ortho.chains

python3 ${scriptsDir}/subsetToYTorthoChains.py $mappings $allChains $outdir/$orthoChainName

gzip $outdir/$orthoChainName
