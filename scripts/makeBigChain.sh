#!/bin/bash

set -beEuo pipefail

function usage {
    cat << EOF
Usage: script.sh chainFile targetAssemblyAbbrev chromSizesFile destinationDir

e.g. makeChainTrack.sh hg38.other.someChain.gz hg38 hg38.chrom.sizes bigChainTracksDir

Assumes that these files are accessible: 
* https://genome.ucsc.edu/goldenpath/help/examples/bigChain.as
* https://genome.ucsc.edu/goldenpath/help/examples/bigLink.as

EOF
    exit 1;
}

if [ $# -ne 4 ]; then
    usage;
fi
destination=$(readlink -f $4)
chainFile=$(readlink -f $1)
assembly=$2
chromSizes=$(readlink -f $3)
outNameRoot=$(basename $chainFile .chains.gz)
outdir=/tmp

scriptsDir=$(realpath $(dirname $0))

mkdir -p $outdir/$outNameRoot".chainTrack"
rm -rf $outdir/$outNameRoot".chainTrack"/*
cd $outdir/$outNameRoot".chainTrack"

# Use the hgLoadChain utility to generate the chain.tab and link.tab files needed to create the bigChain file
hgLoadChain -noBin -test $assembly bigChain $chainFile

# Create the bigChain file from your input chain file using a combination of sed, awk and the bedToBigBed utility
sed 's/\.000000//g' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > $outNameRoot".bigChain"
bedToBigBed -type=bed6+6 -as=${scriptsDir}/bigChain.as -tab $outNameRoot".bigChain" $chromSizes bigChain.bb

# To display your date in the Genome Browser, you must also create a binary indexed link file to accompany your bigChain file
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > bigChain.bigLink
bedToBigBed -type=bed4+1 -as=${scriptsDir}/bigLink.as -tab bigChain.bigLink $chromSizes bigChain.link.bb

mv bigChain.bb $destination/$outNameRoot".bigChain.bb"
mv bigChain.link.bb $destination/$outNameRoot".bigChain.link.bb"

rm -rf $outdir/$outNameRoot".chainTrack"
