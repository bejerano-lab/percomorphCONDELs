#!/bin/bash

set -beEuo pipefail

function usage {
    cat << EOF


Usage: $0 ensgIDmappedChainsFille refAssemblyAbbrev refChr refStart[0-based] refEnd[not included] outdir

e.g. $0 ../processedInputs/orthoChainLookupFiles/0/ENSORLG00000000858 oryLat04 chr1 10063080[0-based] 10063799[not included] ../orthoSequences

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

outFile=${outdir}/consAnchoredOrthoFastas_${refAssemblyAbbrev}_${refChr}_${refStart}_${refEnd}_$(basename $ensgIDchainLookupFile)

rm -f $outFile

random="getConsAnchoredOrthologousFastas_"$(basename $ensgIDchainLookupFile)_$(date "+%Y%m%d-%H%M%S")"-"$(shuf -i1-10000 -n1)
mkdir -p /tmp/${random}
echo -e $refChr"\t"$refStart"\t"$((refStart+1)) > /tmp/${random}/refStartCoord.bed
echo -e $refChr"\t"$((refEnd-1))"\t"$((refEnd)) > /tmp/${random}/refEndCoord.bed

ref2bit=${inputDataRoot}/2bits/${refAssemblyAbbrev}.2bit

#twoBitToFa $ref2bit -seq=$refChr -start=$refStart -end=$refEnd /tmp/$random$refAssemblyAbbrev
#head -1 /tmp/$random$refAssemblyAbbrev >> $outFile
#echo -e $refAssemblyAbbrev"\t""\t"$(tail -n +2 /tmp/$random$refAssemblyAbbrev | tr '[:lower:]' '[:upper:]' | sed -r 's/\s//g')  >> $outFile
#echo -e $refAssemblyAbbrev"\t""\t"$(tail -n +2 /tmp/$random$refAssemblyAbbrev | tr '[:lower:]' '[:upper:]' | sed -r 's/\s//g') 

echo ">"${refAssemblyAbbrev}"."${refChr}":"$((refStart+1))"-"${refEnd} >> $outFile
twoBitToFa $ref2bit -seq=${refChr} -start=${refStart} -end=${refEnd} /tmp/${random}/tmpout
cat /tmp/${random}/tmpout >> $outFile
echo >> $outFile

for speciesInfo in $(tail -n +2 $ensgIDchainLookupFile); do
    chainID=$(echo $speciesInfo | cut -d"_" -f2 | sed 's/chain//')  # e.g. "balAcu1_chain2" becomes "2"
    species=$(echo $speciesInfo | cut -d"_" -f1)   # e.g. "balAcu1_chain2" becomes "balAcu1"
    species2bit=${inputDataRoot}/2bits/${species}.2bit
    currChainSet=${inputDataRoot}/orthoChains/${refAssemblyAbbrev}.${species}*chains.gz
    chainFilter -id=${chainID} $currChainSet > /tmp/${random}/${species}_currChainForLifting
    liftOver /tmp/${random}/refStartCoord.bed /tmp/${random}/${species}_currChainForLifting /tmp/${random}/${species}_liftedRefStartCoord.bed /tmp/${random}/${species}_failedLiftedRefStartCoord.bed
    liftOver /tmp/${random}/refEndCoord.bed /tmp/${random}/${species}_currChainForLifting /tmp/${random}/${species}_liftedRefEndCoord.bed /tmp/${random}/${species}_failedLiftedRefEndCoord.bed

    if [[ -s /tmp/${random}/${species}_liftedRefStartCoord.bed && -s /tmp/${random}/${species}_liftedRefEndCoord.bed ]]; then

        cat /tmp/${random}/${species}_liftedRef*Coord.bed | sort -k2,2n > /tmp/${random}/${species}_queryCoords.bed
        queryStart=$(head -1 /tmp/${random}/${species}_queryCoords.bed | cut -f2)
        queryEnd=$(tail -n 1 /tmp/${random}/${species}_queryCoords.bed | cut -f3)
        queryChr=$(cut -f1 /tmp/${random}/${species}_queryCoords.bed | sort -u)

        echo ">"${species}.${queryChr}:${queryStart}-${queryEnd} >> $outFile
        twoBitToFa -seq=${queryChr} -start=${queryStart} -end=${queryEnd} ${species2bit} /tmp/${random}/tmpout
	cat /tmp/${random}/tmpout >> $outFile
        echo >> $outFile
    else
	echo ">"${species} >> $outFile
	echo >> $outFile
    fi

done

rm -rf /tmp/${random}

