#!/bin/bash -e

set -beEuo pipefail

function usage {
	cat << EOF
Usage: script.sh windowsFile[not_gzipped] e.g. ../processedInputs/windowPctIDs/bolPec01/bolPec01_10bpWindowsSortedByPctID

EOF
	exit 1;
}

if [ $# -ne 1 ]; then
	usage;
fi

scriptsDir=$(realpath $(dirname $0))

windows=$(readlink -f $1)
numWindowsUsed=12000000 # seed number of windows must be less than 5% of oryLat04
covTarget=0.05
outdir=${scriptsDir}/../processedInputs/slidingWindowCons
root=$(basename $windows)
windowSize=$(echo $root | cut -d"_" -f2 | cut -d"b" -f1)
refChromSizes=${scriptsDir}/../processedInputs/chromSizes/oryLat04.chrom.sizes

REF_ASSEMBLY_SIZE=$(cat $refChromSizes | datamash sum 2)  # sum up all chromosome/contig lengths in column 2 to get reference assembly size

cp $0 $outdir/copyOfRunScript_$(basename $0 .sh)_$(date +%Y%m%d_%H%M)

coarse="1000000 500000 100000 10000 1000 100 10 1"
extraCoarse="5000000 1000000 500000 100000 10000 1000 100 10 1"
extraFine="1000 100 10 1"
extraExtraFine="1"
stepSizes="100000 10000 1000 100 10 1"


stepSizes=${coarse}

tempDirRoot="getWindowsCovering5pctOfRef"$(shuf -i1-10000 -n1)$(shuf -i1-10000 -n1)$(shuf -i1-10000 -n1)$root

HOMEDIR=/cluster/u/hchen17

if test -d "$HOMEDIR"; then
	mkdir -p /tmp/$tempDirRoot
	intermediateDir=/tmp/$tempDirRoot
else
	mkdir -p ${outdir}/$tempDirRoot
	intermediateDir=${outdir}/$tempDirRoot
fi

for step in $stepSizes; do
	
	# get top n windows
	head -$numWindowsUsed ${windows} > ${intermediateDir}/$root"_Nwindows"
	
	# save working list of windows from which windows that have already been considered are deleted
	sed '1,'"$numWindowsUsed"'d' ${windows} > ${intermediateDir}/$root"_minusTopN"
	
	# bedSort top n windows
	sort -k1,1 -k2,2g ${intermediateDir}/$root"_Nwindows" > ${intermediateDir}/$root"_NwindowsSorted"
	
	# merge top n windows and print length of each element in column 5
	mergeBed -i ${intermediateDir}/$root"_NwindowsSorted" -c 4 -o distinct | awk '{print $0 "\t" $3-$2}' > ${intermediateDir}/$root"_NwindowsSortedMerged"
	# some merged regions will come from different chains! known caveat of this method

	# determine what fraction of the reference genome is covered by the n windows (when flattened to 1x coverage)
	currBasePairCov=$(cat ${intermediateDir}/$root"_NwindowsSortedMerged"| datamash sum 5)
	currCov=$(echo $currBasePairCov $REF_ASSEMBLY_SIZE | awk '{print $1/$2}')

	# generate boolean-like flag to determine whether more windows need to be added to approach target coverage fraction (i.e. 5% or 10%)
	needMoreCov=$(echo $covTarget $currCov | awk '{ if ($1 > $2) { print "1" } else { print "0" } }')
	echo "NUM WINDOWS USED:" $numWindowsUsed "CURRENT COVERAGE:" $currCov "NEED MORE?:" $needMoreCov "STEP SIZE:" $step >> $outdir/$root"_breadcrumbs"
	echo "REF_ASSEMBLY_SIZE:" $REF_ASSEMBLY_SIZE >> $outdir/$root"_breadcrumbs"

	while [ $needMoreCov -eq 1 ]; do

		# print out salient values for easier debugging and QA
		echo "NUM WINDOWS USED:" $numWindowsUsed "CURRENT COVERAGE:" $currCov "NEED MORE?:" $needMoreCov "STEP SIZE:" $step >> $outdir/$root"_breadcrumbs"

		# combine current merged windows with stepSize more windows; then sort and merge/flatten
		cat ${intermediateDir}/$root"_NwindowsSortedMerged" > ${intermediateDir}/$root"_Nwindows"
		head -$step ${intermediateDir}/$root"_minusTopN" >> ${intermediateDir}/$root"_Nwindows"
		sort -k1,1 -k2,2g ${intermediateDir}/$root"_Nwindows" > ${intermediateDir}/$root"_NwindowsSorted"

		# merge top n windows and print length of each element in column 5
		mergeBed -i ${intermediateDir}/$root"_NwindowsSorted" -c 4 -o distinct | awk '{print $0 "\t" $3-$2}' > ${intermediateDir}/$root"_NwindowsSortedMerged"
		# some merged regions will come from different chains! known caveat of this method

		# determine what fraction of the reference genome is covered by the n windows (when flattened to 1x coverage)
		currBasePairCov=$(cat ${intermediateDir}/$root"_NwindowsSortedMerged"| datamash sum 5)
		currCov=$(echo $currBasePairCov $REF_ASSEMBLY_SIZE | awk '{print $1/$2}')

		if [[ $currCov =~ ^[[:space:]] ]]; then  # error because disk space is full
			echo "*****PROBABLY DISK SPACE FULL*****" >> $outdir/$root"_breadcrumbs"
			exit 1
		fi

		# generate boolean-like flag to determine whether more windows need to be added to approach target coverage fraction (i.e. 5% or 10%)
		needMoreCov=$(echo $covTarget $currCov | awk '{ if ($1 > $2) { print "1" } else { print "0" } }')

		# in place, delete the most recent set of windows that have just been considered

		if [[ $needMoreCov -eq 1 ]]; then  # i.e. don't add stepSize and don't delete top window if overshot coverage target 
			numWindowsUsed=$[numWindowsUsed + $step]
			sed -i '1,'"$step"'d' ${intermediateDir}/$root"_minusTopN"
		fi
	done
done

# print out salient values for easier debugging and QA
echo "NUM WINDOWS USED:" $numWindowsUsed "CURRENT COVERAGE:" $currCov "NEED MORE?:" $needMoreCov "STEP SIZE:" $step >> $outdir/$root"_breadcrumbs"

# get ties
pctID=$(head -1 ${intermediateDir}/$root"_minusTopN" | awk '{print $5}') # save pctID for current top window
awk -v num="$pctID" '$5 == num' ${intermediateDir}/$root"_minusTopN" > ${intermediateDir}/$root"_ties"
numTies=$(wc -l ${intermediateDir}/$root"_ties" | awk '{print $1}')

numWindowsUsed=$[numWindowsUsed + $numTies]

cat ${intermediateDir}/$root"_ties" > ${intermediateDir}/$root"_Nwindows"
cat ${intermediateDir}/$root"_NwindowsSortedMerged" >> ${intermediateDir}/$root"_Nwindows"
sort -k1,1 -k2,2g ${intermediateDir}/$root"_Nwindows" > ${intermediateDir}/$root"_NwindowsSorted"

# merge top n windows and print length of each element in column 4
	mergeBed -i ${intermediateDir}/$root"_NwindowsSorted" -c 4 -o distinct | awk '{print $0 "\t" $3-$2}' > ${intermediateDir}/$root"_NwindowsSortedMerged"
# some merged regions will come from different chains! known caveat of this method

# determine what fraction of the reference genome is covered by the new set of windows (when flattened to 1x coverage)
currBasePairCov=$(cat ${intermediateDir}/$root"_NwindowsSortedMerged"| datamash sum 5)
currCov=$(echo $currBasePairCov $REF_ASSEMBLY_SIZE | awk '{print $1/$2}')

if [[ $currCov =~ ^[[:space:]] ]]; then  # error because disk space is full
	echo "*****PROBABLY DISK SPACE FULL*****" >> $outdir/$root"_breadcrumbs"
	exit 1
fi

# print out salient values for easier debugging and QA
echo "--------------------" >> $outdir/$root"_breadcrumbs"
echo "numTiesAdded:" $numTies >> $outdir/$root"_breadcrumbs"
echo "lastWindowPctID:" $pctID >> $outdir/$root"_breadcrumbs"
echo "numWindowsUsed:" $numWindowsUsed >> $outdir/$root"_breadcrumbs"
echo "currCov:" $currCov >> $outdir/$root"_breadcrumbs"

#newRoot=$root"_"$covTarget"coverage_"$numWindowsUsed"windows"

#awk 'OFS="\t" {print $1,$2,$3,$4}' ${intermediateDir}/$root"_NwindowsSortedMerged" > $outdir/$newRoot


### NOW COLLAPSE WINDOWS BY CHAIN ###
species=$(basename $windows | cut -d"_" -f1)

out=$outdir/$species/$windowSize
mkdir -p $out

head -$numWindowsUsed $windows > ${intermediateDir}/$root"_"$numWindowsUsed"windows.unsorted"

# sort windows by chain ID and save in working file
sort -k4,4 ${intermediateDir}/$root"_"$numWindowsUsed"windows.unsorted" > ${intermediateDir}/windowsSortedByChainID
rm ${intermediateDir}/$root"_"$numWindowsUsed"windows.unsorted"

while [ -s ${intermediateDir}/windowsSortedByChainID ]; do # while there are still windows

	currentChainID=$(head -1 ${intermediateDir}/windowsSortedByChainID | awk '{print $4}')  # save current chainID
	cat ${intermediateDir}/windowsSortedByChainID | grep "$currentChainID"  > ${intermediateDir}/$currentChainID"unmerged.bed"  # save all windows on that chain
	sort -k1,1 -k2,2g ${intermediateDir}/$currentChainID"unmerged.bed" > ${intermediateDir}/$currentChainID"unmergedSorted.bed"  # sort in prep for merging
	mergeBed -i ${intermediateDir}/$currentChainID"unmergedSorted.bed" > ${intermediateDir}/$currentChainID".bed"  # flatten all windows on that chain
	numLines=$(wc -l ${intermediateDir}/$currentChainID"unmergedSorted.bed" | cut -d" " -f1)  # determine how many windows were just saved
	cat ${intermediateDir}/$currentChainID".bed" > $out/$currentChainID".bed"  # copy final output to /cluster
	sed -i '1,'"$numLines"'d' ${intermediateDir}/windowsSortedByChainID  # delete already saved windows from working file
	rm ${intermediateDir}/$currentChainID*  # remove intermediate files
done

rm -rf ${intermediateDir}
