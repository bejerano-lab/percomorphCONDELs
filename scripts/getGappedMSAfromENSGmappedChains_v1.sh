#!/bin/bash
set -beEuo pipefail

function usage {
    cat << EOF
Usage: $0 ensgIDmappedChainsFille refAssemblyAbbrev refChr refStart refEnd outdir

e.g. $0 ../processedInputs/orthoChainLookupFiles/0/ENSORLG00000000858 oryLat04 chr1 10063080 10063799 ../MSA


where ensgIDmappedChainsFile looks like

chr1  10063000 10063001
ampCit01_chain249
cynSem_chain327
dicLab01_chain9
funHet01_chain292
hapBur01_chain144
hipCom01_chain510
...

EOF
    exit 1;
}

if [ $# -ne 6 ]; then
    usage;
fi

scriptsDir=$(realpath $(dirname $0))
ensgIDchainLookupFile=$(readlink -f $1)
refAssemblyAbbrev=$2
refChr=$3
refStart=$4
refEnd=$5
outdir=$(readlink -f $6)

inputDataRoot=$(readlink -f ${scriptsDir}/../processedInputs/)

outFile=${outdir}/gappedMSAfromENSGmappedChains_${refAssemblyAbbrev}_${refChr}_${refStart}_${refEnd}_$(basename $ensgIDchainLookupFile)

rm -f $outFile

random="getGappedMSAfromENSGmappedChains"$(shuf -i1-10000 -n1)$(shuf -i1-10000 -n1)$(shuf -i1-10000 -n1)"_"

ref2bit=${inputDataRoot}/2bits/${refAssemblyAbbrev}.2bit

twoBitToFa $ref2bit -seq=$refChr -start=$refStart -end=$refEnd /tmp/${random}${refAssemblyAbbrev}
head -1 /tmp/$random$refAssemblyAbbrev >> $outFile
echo -e ${refAssemblyAbbrev}"\t""\t"$(tail -n +2 /tmp/${random}${refAssemblyAbbrev} | tr '[:lower:]' '[:upper:]' | sed -r 's/\s//g')  >> $outFile
echo -e ${refAssemblyAbbrev}"\t""\t"$(tail -n +2 /tmp/${random}${refAssemblyAbbrev} | tr '[:lower:]' '[:upper:]' | sed -r 's/\s//g') 

for speciesInfo in $(tail -n +2 $ensgIDchainLookupFile); do
    chainID=$(echo $speciesInfo | cut -d"_" -f2 | sed 's/chain//')  # e.g. "balAcu1_chain2" becomes "2"
    species=$(echo $speciesInfo | cut -d"_" -f1)   # e.g. "balAcu1_chain2" becomes "balAcu1"
    species2bit=${inputDataRoot}/2bits/${species}.2bit
    currChainSet=${inputDataRoot}/orthoChains/${refAssemblyAbbrev}.${species}*chains.gz
    ${scriptsDir}/filterChainToAxt.sh $chainID $currChainSet $ref2bit $species2bit > /tmp/${random}${species}_chain${chainID}
    python3 ${scriptsDir}/getPairwiseAxtSubstring_identityDots.py /tmp/${random}${species}_chain${chainID} $refStart $refEnd /tmp/${random}${species}
    echo $species >> $outFile
    cat /tmp/${random}${species} >> $outFile
    echo >> $outFile

    echo $species
    cat /tmp/${random}${species}
    echo
    #echo -e $species"\t""\t"$(tail -n 1 /tmp/$random$species) >> $outFile
    #echo -e $species"\t""\t"$(cat /tmp/$random$species)
done

rm /tmp/$random*
