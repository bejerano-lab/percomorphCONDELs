#!/bin/bash

set -beEuo pipefail

function usage {
	cat << EOF
Usage: script.sh DEF_file
where DEF file defines these parameters, e.g....
	lookupDir=/cluster/u/hchen17/treeWAS/inputData_medium/percomorphOrthoChainLookupFiles/allGenesAllSpecies
	
	tetraodontiformes="molMol01 takFla02 takRub01 tetNig2"
	eel="monAlb01"
	sole="cynSem"
	syngnathids="synSco01 hipCom01 hipEre01"
	targetCladeReq=2
	requireMonAlb=False
	requireGasAcu=True
	outgroups="ampCit01 dicLab01 funHet01 gasAcu14 hapBur01 kryMar01 larCro01 latCal01 mayZeb03 miiMii01 neoBri01 notFur02 oreNil02 punNye01 serDum01 serQui01 xipMac02"
	numReqOutgroups=17
	mergeDist=20
	tssPadding=500000
	delsDir=/cluster/u/hchen17/treeWAS/inputData/percomorphDels
	consWindowsDir=/cluster/u/hchen17/treeWAS/inputData/percomorphSlidingWindowCons
	intactAxtsDir=/cluster/u/hchen17/treeWAS/inputData/percomorphIntactAxtBlocks
	minConSize=20
	minConDelSize=20
	outDir=/cluster/u/hchen17/treeWAS/output/percomorph/
	maxViolations=0
	numSygnathids=2
	numTetraodontiformes=3
EOF
	exit 1;
}

if [ $# -ne 1 ]; then
	usage;
fi

scriptsDir=$(realpath $(dirname $0))

DEF=$(readlink -f $1)

syngnathids=$(cat $DEF | grep "syngnathids=" | cut -d"=" -f2 | sed 's/"//g')
tetraodontiformes=$(cat $DEF | grep "tetraodontiformes=" | cut -d"=" -f2 | sed 's/"//g')

axtDir=${scriptsDir}/../processedInputs/axtBlocks

oryLat04ChromSizes=${scriptsDir}/../processedInputs/chromSizes/oryLat04.chrom.sizes

tmpDir=/tmp/$(date +%Y%m%d_%H%M)_findTwoOrMoreSynAndTetIntactAxtBlocks
mkdir -p $tmpDir

# combine all axt blocks
for species in $syngnathids $tetraodontiformes; do
	cat ${axtDir}/${species}/${species}*.bed >> ${tmpDir}/simplyCombinedAxtBlocks.bed
done

sort -k1,1 ${tmpDir}/simplyCombinedAxtBlocks.bed > ${tmpDir}/sortedAxtBlocks.bed

genomeCoverageBed -bg -i ${tmpDir}/sortedAxtBlocks.bed -g $oryLat04ChromSizes  | awk '$4 >= 2' > ${scriptsDir}/../processedInputs/axtBlocks/twoOrMoreSynAndTetIntact.bed

rm -r $tmpDir

