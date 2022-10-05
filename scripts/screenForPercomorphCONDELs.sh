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

reauth

HOMEDIR=/cluster/u/hchen17

#NUM_SPECIES=38 # INCLUDING oryLat04

scriptsDir=$(realpath $(dirname $0))

FILTER_BEDS=$(readlink -f ${scriptsDir}/../processedInputs/filterBEDs)

nonOrthoTargetAxtsBed=$(readlink -f ${scriptsDir}/../processedInputs/axtBlocks/twoOrMoreTargetsIntact.bed)

oryLat04ChromSizes=$(readlink -f i${scriptsDir}/../processedInputs/chromSizes/oryLat04.chrom.sizes)
oryLat04ChromIDregEx="chr[1-9]+"
oryLat04_codingExons=${FILTER_BEDS}/oryLat04_ensembl98_codingExons.bed
oryLat04_nonCodingExons=${FILTER_BEDS}/oryLat04_ensembl98_nonCodingExons.bed
oryLat04_3UTR=${FILTER_BEDS}/oryLat04_ensembl98_3UTR.bed
oryLat04_5UTR=${FILTER_BEDS}/oryLat04_ensembl98_5UTR.bed
oryLat04_introns=${FILTER_BEDS}/oryLat04_ensembl98_introns.bed

windowSizes="10 25 50 100"
out=$(readlink -f $(cat $1 | grep "outDir=" | cut -d"=" -f2))/$(date +%Y%m%d_%H%M)"_RESULTS_"$(basename $1 .txt | cut -d"_" -f2-)
mkdir -p $out/par $out/uncleanedResultFiles $out/logFiles

cp $1 $out/DEF
cp $(readlink -f $0) $out/logFiles/copyOf$(basename $0)  # save a copy of this running script

# screen parameters definition file
DEF=$(readlink -f $out/DEF)

lookupDir=$(readlink -f $(cat $DEF | grep "lookupDir=" | cut -d"=" -f2))
tssPadding=$(cat $DEF | grep "tssPadding=" | cut -d"=" -f2)
delsDir=$(readlink -f $(cat $DEF | grep "delsDir=" | cut -d"=" -f2))
consDir=$(readlink -f $(cat $DEF | grep "consWindowsDir=" | cut -d"=" -f2))
axtDir=$(readlink -f $(cat $DEF | grep "intactAxtsDir=" | cut -d"=" -f2))
minConSizeInBp=$(cat $DEF | grep "minConSize=" | cut -d"=" -f2)
minConDelSizeInBp=$(cat $DEF | grep "minConDelSize=" | cut -d"=" -f2)

numReqOutgroups=$(cat $DEF | grep "numReqOutgroups=" | cut -d"=" -f2)
outgroups=$(cat $DEF | grep "outgroups=" | cut -d"=" -f2)

mergeDist=$(cat $DEF | grep "mergeDist=" | cut -d"=" -f2)

targetCladeReq=$(cat $DEF | grep "targetCladeReq=" | cut -d"=" -f2)
syngnathids=$(cat $DEF | grep "syngnathids=" | cut -d"=" -f2)
eel=$(cat $DEF | grep "eel=" | cut -d"=" -f2)
tetraodontiformes=$(cat $DEF | grep "tetraodontiformes=" | cut -d"=" -f2)

numTetraodontiformes=$(cat $DEF | grep "numTetraodontiformes=" | cut -d"=" -f2)
numSyngnathids=$(cat $DEF | grep "numSyngnathids=" | cut -d"=" -f2)

sole=$(cat $DEF | grep "sole=" | cut -d"=" -f2)

maxViolations=$(cat $DEF | grep "maxViolations=" | cut -d"=" -f2)
requireMonAlb=$(cat $DEF | grep "requireMonAlb=" | cut -d"=" -f2)
requireGasAcu=$(cat $DEF | grep "requireGasAcu=" | cut -d"=" -f2)

runningHost=$(hostname)

echo $runningHost $lookupDir $tssPadding $delsDir $consDir $axtDir $minConDelSizeInBp $minConSizeInBp $numReqOutgroups $outgroups $mergeDist $targetCladeReq $syngnathids $tetraodontiformes $numTetraodontiformes $eel $sole $maxViolations $requireMonAlb

cat << EOF > $out/logFiles/parametersReceivedFromDEF
host 					$runningHost

lookupDir				$lookupDir
tssPadding				$tssPadding
delsDir 				$delsDir
consDir 				$consDir
axtDir 					$axtDir
minConDelSizeInBp		$minConDelSizeInBp
minConSizeInBp			$minConSizeInBp

numReqOutgroups			$numReqOutgroups
outgroups 				$outgroups

mergeDist				$mergeDist

targetCladeReq 			$targetCladeReq
syngnathids				$syngnathids
tetraodontiformes		$tetraodontiformes
numTetraodontiformes 	$numTetraodontiformes
numSyngnathids			$numSyngnathids
eel						$eel
sole 					$sole

maxViolations			$maxViolations
requireMonAlb 			$requireMonAlb
requireGasAcu			$requireGasAcu


EOF

for windowSize in $windowSizes; do
	num=0
	while [[ $num -lt 31 ]]; do 
		mkdir -p $out/$windowSize/$num #$out/TEMP${windowSize}
		num=$(( num + 1 ))
	done
done

# write out 4 scripts to $out/par -- one script per window size

for windowSize in $windowSizes; do

TEMP_DIR=/tmp

echo -e "TEMP_DIR\t" $TEMP_DIR >> $out/logFiles/parametersReceivedFromDEF

cat << EOF > $out/par/"analyzeENSG_"$windowSize"_Window.sh" # gene delDir conDir outDir timestamp
#/!/bin/bash -e

set -beEu -o pipefail

ensg=\$(basename \$1)
outDir=\$2
uniqID=\$3

requireMonAlb=$requireMonAlb
requireGasAcu=$requireGasAcu

rm -rf $TEMP_DIR/"analyzeENSG_"\$ensg"_"\$uniqID/
mkdir -p $TEMP_DIR/"analyzeENSG_"\$ensg"_"\$uniqID/

tmpFileRoot=$TEMP_DIR/"analyzeENSG_"\$ensg"_"\$uniqID/\$ensg"_"\$uniqID

lookup_file=\$(readlink -f \$1)

cat \$lookup_file | egrep $oryLat04ChromIDregEx > \$tmpFileRoot"_tss.bed"
slopBed -i \$tmpFileRoot"_tss.bed" -g $oryLat04ChromSizes -b $tssPadding > \$tmpFileRoot"_interval.bed"


#touch \$tmpFileRoot".chainGaps"
touch \$tmpFileRoot".OGchainGaps"
touch \$tmpFileRoot".cons"
touch \$tmpFileRoot".axtBlocks"


canProceed=False

if [ \$requireGasAcu = True ]; then 
	speciesInfo=\$(cat \$lookup_file | grep gasAcu || true) # || true prevents exit code of 1 when grep returns no match


	if [[ \$speciesInfo != "" ]]; then

		# save off curr species name
		species=\$(echo \$speciesInfo | cut -d"_" -f1)

		# save off ID for the chain that has the curr species ortho alignment for this gene
		chainID="_"\$(echo \$speciesInfo | cut -d"_" -f2)"_"

		axtElementsFile=$axtDir/\$species/\$speciesInfo"_.bed"

		if [ -e \$axtElementsFile ]; then
			intersectBed -wa -u -a \$axtElementsFile -b \$tmpFileRoot"_interval.bed" | mergeBed -d $mergeDist -i stdin | awk -v var=\$speciesInfo '{print \$0 "\t" var"_"}' > \$tmpFileRoot".GasAcuAxtBlocks"
			if [[ -s \$tmpFileRoot".GasAcuAxtBlocks" ]]; then
				canProceed=True
			fi
		fi

	fi

else
	canProceed=True
fi


if [ \$canProceed = False ]; then
	touch \$tmpFileRoot".CONDELs"
	# even if no ConDels are found, move to outdir as record of job completion
	mv \$tmpFileRoot".CONDELs" \$outDir/
	exit
fi



## CONSIDER TARGETS BY CLADE ##

for speciesInfo in \$(cat \$lookup_file | egrep -e \$(echo $tetraodontiformes | sed -e 's/[[:space:]]/|/g')); do
	
	# save off curr species name
	species=\$(echo \$speciesInfo | cut -d"_" -f1)

	# save off ID for the chain that has the curr species ortho alignment for this gene
	chainID="_"\$(echo \$speciesInfo | cut -d"_" -f2)"_"

	# get file with all dels for this species' gene's ortho chains
	delFile=$delsDir/\$species/\$speciesInfo"_.bed"

	if [ -e \$delFile ]; then 

		# then subset to only those dels within the x Mb interval ± ENSG TSS on the proper ortho chain
		intersectBed -wa -u -a \$delFile -b \$tmpFileRoot"_interval.bed" >> \$tmpFileRoot"_tetraodontiformes.unmergedChainGaps"
	fi
done

for speciesInfo in \$(cat \$lookup_file | egrep -e \$(echo $syngnathids | sed -e 's/[[:space:]]/|/g')); do
	
	# save off curr species name
	species=\$(echo \$speciesInfo | cut -d"_" -f1)

	# save off ID for the chain that has the curr species ortho alignment for this gene
	chainID="_"\$(echo \$speciesInfo | cut -d"_" -f2)"_"

	# get file with all dels for this species' gene's ortho chains
	delFile=$delsDir/\$species/\$speciesInfo"_.bed"

	if [ -e \$delFile ]; then 

		# then subset to only those dels within the x Mb interval ± ENSG TSS on the proper ortho chain
		intersectBed -wa -u -a \$delFile -b \$tmpFileRoot"_interval.bed" >> \$tmpFileRoot"_syngnathids.unmergedChainGaps"
	fi
done

if [ $eel != "" ]; then
	for speciesInfo in \$(cat \$lookup_file | egrep -e \$(echo $eel | sed -e 's/[[:space:]]/|/g')); do
		
		# save off curr species name
		species=\$(echo \$speciesInfo | cut -d"_" -f1)

		# save off ID for the chain that has the curr species ortho alignment for this gene
		chainID="_"\$(echo \$speciesInfo | cut -d"_" -f2)"_"

		# get file with all dels for this species' gene's ortho chains
		delFile=$delsDir/\$species/\$speciesInfo"_.bed"

		if [ -e \$delFile ]; then 
			# then subset to only those dels within the x Mb interval ± ENSG TSS on the proper ortho chain
			intersectBed -wa -u -a \$delFile -b \$tmpFileRoot"_interval.bed" >> \$tmpFileRoot".monAlb01chainGaps"
		fi
	done
fi

if [ $sole != "" ]; then

	for speciesInfo in \$(cat \$lookup_file | egrep -e \$(echo $sole | sed -e 's/[[:space:]]/|/g')); do
		
		# save off curr species name
		species=\$(echo \$speciesInfo | cut -d"_" -f1)

		# save off ID for the chain that has the curr species ortho alignment for this gene
		chainID="_"\$(echo \$speciesInfo | cut -d"_" -f2)"_"

		# get file with all dels for this species' gene's ortho chains
		delFile=$delsDir/\$species/\$speciesInfo"_.bed"

		if [ -e \$delFile ]; then 
			# then subset to only those dels within the x Mb interval ± ENSG TSS on the proper ortho chain
			intersectBed -wa -u -a \$delFile -b \$tmpFileRoot"_interval.bed" >> \$tmpFileRoot".cynSemChainGaps"
		fi
	done

fi

existingTargetClades=0

if [ -e \$tmpFileRoot".monAlb01chainGaps" ]; then
	cat \$tmpFileRoot".monAlb01chainGaps" >> \$tmpFileRoot".chainGapsByTargetClade"
	let "existingTargetClades=existingTargetClades+1"
fi

if [ -e \$tmpFileRoot".cynSemChainGaps" ]; then
	cat \$tmpFileRoot".cynSemChainGaps" >> \$tmpFileRoot".chainGapsByTargetClade"
	let "existingTargetClades=existingTargetClades+1"
fi

if [ -e \$tmpFileRoot"_tetraodontiformes.unmergedChainGaps" ]; then
	genomeCoverageBed -bg -i \$tmpFileRoot"_tetraodontiformes.unmergedChainGaps" -g $oryLat04ChromSizes  | awk '\$4 >= $numTetraodontiformes' > \$tmpFileRoot"_tetraodontiformes.unmergedChainGaps.bedGraph"
	
	if [ -s \$tmpFileRoot"_tetraodontiformes.unmergedChainGaps.bedGraph" ]; then
		mergeBed -c 4 -o max -i \$tmpFileRoot"_tetraodontiformes.unmergedChainGaps.bedGraph" >> \$tmpFileRoot".chainGapsByTargetClade"
		let "existingTargetClades=existingTargetClades+1"

	fi
fi

if [ -e \$tmpFileRoot"_syngnathids.unmergedChainGaps" ]; then
	genomeCoverageBed -bg -i \$tmpFileRoot"_syngnathids.unmergedChainGaps" -g $oryLat04ChromSizes  | awk '\$4 >= $numSyngnathids' > \$tmpFileRoot"_syngnathids.unmergedChainGaps.bedGraph"
	if [ -s  \$tmpFileRoot"_syngnathids.unmergedChainGaps.bedGraph" ]; then
		let "existingTargetClades=existingTargetClades+1"
		mergeBed -c 4 -o max -i \$tmpFileRoot"_syngnathids.unmergedChainGaps.bedGraph" >> \$tmpFileRoot".chainGapsByTargetClade"
	fi
fi


## PROCEED ONLY IF THERE ARE ENOUGH TARGET CLADES REPRESENTED and if required gasAcu files exist ##

if [ \$existingTargetClades -ge $targetCladeReq ] && [ \$canProceed = True ]; then

	## CONSIDER OUTGROUPS ##
	for speciesInfo in \$(cat \$lookup_file | egrep \$(echo $outgroups | sed -e 's/[[:space:]]/|/g') | egrep -v $oryLat04ChromIDregEx); do

		# save off curr species name
		species=\$(echo \$speciesInfo | cut -d"_" -f1)

		# save off ID for the chain that has the curr species ortho alignment for this gene
		chainID="_"\$(echo \$speciesInfo | cut -d"_" -f2)"_"

		# get file with all dels for this species' gene's ortho chains
		delFile=$delsDir/\$species/\$speciesInfo"_.bed"

		# get file with all conserved windows for this species' ortho chain
		conFile=$consDir/\$species/$windowSize/\$speciesInfo"_.bed"

		# get file with all intact alignment regions for this species' ortho chain
		axtElementsFile=$axtDir/\$species/\$speciesInfo"_.bed"

		# then subset to only those conserved windows INTERSECTING! the x Mb interval ± ENSG TSS on the proper ortho chain
		# won't trim conserved region that intersects but extends beyond the x Mb interval

		if [ -e \$delFile ] && [ -e \$conFile ] && [ -e \$axtElementsFile ]; then
			intersectBed -wa -u -a \$delFile -b \$tmpFileRoot"_interval.bed" >> \$tmpFileRoot".OGchainGaps"
			intersectBed -wa -u -a \$conFile -b \$tmpFileRoot"_interval.bed" | awk -v var=\$speciesInfo '{print \$0 "\t" var"_"}' >> \$tmpFileRoot".cons"
			intersectBed -wa -u -a \$axtElementsFile -b \$tmpFileRoot"_interval.bed" >> \$tmpFileRoot".axtBlocks"
		fi

	done

	# genomeCoverageBed requires file to be grouped by chromosome so no need for further sorting; output from above should all be from same chain/chromosome

	# get all regions that have highly conserved sequence in ≥ numReqOutgroups that spans more than minConSizeInBp
	genomeCoverageBed -bg -i \$tmpFileRoot".cons" -g $oryLat04ChromSizes | awk '\$4 >= $numReqOutgroups' | mergeBed -d $mergeDist -i stdin | sort -k1,1 -k2,2n | awk '\$3-\$2 >= $minConSizeInBp' > \$tmpFileRoot".positionsOfSufConsCoveragePreFilter"

	# get all regions that are deleted in > n=maxViolations outgroups
	genomeCoverageBed -bg -i \$tmpFileRoot".OGchainGaps" -g $oryLat04ChromSizes | awk '\$4 > $maxViolations' | mergeBed -i stdin > \$tmpFileRoot".ogDels"

	subtractBed -a \$tmpFileRoot".positionsOfSufConsCoveragePreFilter" -b \$tmpFileRoot".ogDels" > \$tmpFileRoot".positionsOfSufConsCoverage"

	# generate a file that should have bed entry for every region in x Mb window in which at least one member from each of n=numReqOutgroups possible outgroup clades is intact
	genomeCoverageBed -bg -i \$tmpFileRoot".axtBlocks" -g $oryLat04ChromSizes  | awk '\$4 >= $numReqOutgroups' | mergeBed -i stdin | awk '\$3-\$2 >= $minConSizeInBp' > \$tmpFileRoot".positionsOfSufAxtRegionCoverage"

	# get del regions that are missing in {≥targetCladeReq} target clades 
	genomeCoverageBed -bg -i \$tmpFileRoot".chainGapsByTargetClade" -g $oryLat04ChromSizes | awk '\$4 >= $targetCladeReq' | mergeBed -i stdin > \$tmpFileRoot".targetDels"

	if [ \$requireMonAlb = True ]; then 
		# require that dels exist in rice eel
		intersectBed -a \$tmpFileRoot".targetDels" -b \$tmpFileRoot".monAlb01chainGaps" > \$tmpFileRoot".targetDelsTEMPORARY"
		mv  \$tmpFileRoot".targetDelsTEMPORARY"  \$tmpFileRoot".targetDels"
	fi

	# get shared TargetDels (Target-specific) regions that are missing in all target groups but intact in ≥ n=numReqOutgroups subset of outgroup clades
	# TODO: CONSIDER SIZE OF OVERLAP
	intersectBed -a \$tmpFileRoot".targetDels" -b \$tmpFileRoot".positionsOfSufAxtRegionCoverage" > \$tmpFileRoot".dels"

	# get raw ConDels; don't filter ConDels by size right now	TODO: CONSIDER SIZE OF OVERLAP
	intersectBed -a \$tmpFileRoot".dels" -b \$tmpFileRoot".positionsOfSufConsCoverage" > \$tmpFileRoot".CONDELsAnySizePreGasAcuFilter"

	if [ \$requireGasAcu = True ]; then
		# require that axt block exists in gasAcu
		intersectBed -a \$tmpFileRoot".CONDELsAnySizePreGasAcuFilter" -b \$tmpFileRoot".GasAcuAxtBlocks" > \$tmpFileRoot".CONDELsAnySizeTmp"
		mv   \$tmpFileRoot".CONDELsAnySizeTmp"  \$tmpFileRoot".CONDELsAnySizePreGasAcuFilter"
	fi

	cat \$tmpFileRoot".CONDELsAnySizePreGasAcuFilter" | sort -k1,1 -k2,2n | mergeBed -d $mergeDist -i stdin | sort -k1,1 -k2,2n | awk -F"\t" -v var="\$ensg" 'OFS="\t" {print \$1,\$2,\$3,\$3-\$2"bp_"var}' > \$tmpFileRoot".CONDELsAnySize"

	touch \$tmpFileRoot".CONDELs"

	while read rawCONDEL; do

		size=\$(echo \$rawCONDEL | awk '{print \$3-\$2}')

		if [ \$size -ge $minConDelSizeInBp ]; then
			echo "\$rawCONDEL" >> \$tmpFileRoot".CONDELs"
			# echo \$rawCONDEL >> \$outDir/\$ensg"_"\$uniqID".CONDELs"
		# else
			# echo "\$rawCONDEL" >> \$tmpFileRoot".smallCONDELs"
		fi

	done < \$tmpFileRoot".CONDELsAnySize"

	if [[ -s \$tmpFileRoot".CONDELs" ]]; then # if >0 CONDELs are predicted, then...
		# use big enough CONDELs (if they exist) to get entire original lesion
		intersectBed -wa -u -a \$tmpFileRoot".targetDels"  -b \$tmpFileRoot".CONDELs" | awk -F"\t" -v var="\$ensg" 'OFS="\t" {print \$1,\$2,\$3,\$3-\$2"bp_"var}'> \$outDir/\$ensg"_"\$uniqID".DELs"
	fi

	# even if no ConDels are found, move to outdir as record of job completion
	mv \$tmpFileRoot".CONDELs" \$outDir/

	# if [[ -s \$tmpFileRoot".smallCONDELs" ]]; then 
	# 	# keep any too-small CONDELs that intersect coding regions
	# 	intersectBed -wa -u -a \$tmpFileRoot".smallCONDELs" -b $oryLat04_5UTR $oryLat04_3UTR $oryLat04_codingExons $oryLat04_introns  > \$tmpFileRoot".smallGenicCONDELs"

	# 	if [[ -s \$tmpFileRoot".smallGenicCONDELs" ]]; then
	# 		mv \$tmpFileRoot".smallGenicCONDELs" \$outDir/
	# 	fi
	# fi

else
	touch \$tmpFileRoot".CONDELs"
	# even if no ConDels are found, move to outdir as record of job completion
	mv \$tmpFileRoot".CONDELs" \$outDir/
fi

rm -rf $TEMP_DIR/"analyzeENSG_"\$ensg"_"\$uniqID/
EOF

done

# make scripts executable
chmod +x $out/par/analyzeENSG*.sh

counter=$(shuf -i 1-10000 -n 1)
minReqChains=$(( $numReqOutgroups + $targetCladeReq ))

for lookupFile in $(find $lookupDir -name ENSORLG*); do
	numLines=$(wc -l $lookupFile | cut -d" " -f1)
	if [ $numLines -gt $minReqChains ]; then

		if [[ ($requireMonAlb = True  &&  $(grep monAlb01 $lookupFile | wc -l) -gt 0)  ||  ($requireMonAlb = False) ]]; then
			echo $(basename $lookupFile) >> $out/logFiles/genesConsidered.txt
			gene=$(basename $lookupFile .lookup)
			number=$(echo $gene | sed 's/ENSORLG//' | sed "s/^0\+//")
			folder=$(($number/1000))

			for windowSize in $windowSizes; do
				echo $out/par/analyzeENSG_$windowSize"_Window.sh" $lookupFile $out/$windowSize/$folder $windowSize"bp_"$counter "{check out exists" $out/$windowSize/$folder/$gene"_"$windowSize"bp_"$counter".CONDELs}">> $out/par/jobList
				let counter=counter+1
			done
		else
			echo $(basename $lookupFile) >> $out/logFiles/genesNOTConsidered.txt
		fi

	else
		echo $(basename $lookupFile) >> $out/logFiles/genesNOTConsidered.txt
	fi
done


mv $out/par/jobsOfInterest $out/par/jobList


#read -rsp $'Press any key to start jobs...\n' -n1 key

if test -d "$HOMEDIR"; then
	#must use para make to wait for jobs to finish!
	ssh hoxa "cd $out/par && para make jobList" 2>&1 | tee -a $out/logFiles/parasolErrOut   # append and print stderr and stdout from para make to the error log file
else
	cat $out/par/jobList | parallel
fi
#read -rsp $'Press any key to continue with collating results...\n' -n1 key

numOutFiles=$(ls $out/10/* $out/25/* $out/50/* $out/100/* | grep "\.CONDELs" | wc -l)
numJobs=$(cat $out/par/jobList* | wc -l)

echo "Num outfiles:" $numOutFiles
echo "Num jobs:" $numJobs

if [[ $numJobs -eq $numOutFiles ]]; then
	echo "Collating results..."
	for size in $windowSizes; do
		num=0
		while [[ $num -le 30 ]]; do 
			cat $out/$size/$num/*".CONDELs" >> $out/uncleanedResultFiles/oryLat04.CONDELsCombinedUnmerged_$size"bp"

			numlines=$(find $out/$size/$num/ -name "*.DELs" | wc -l)
			if [[ $numlines -gt "0" ]]; then
				cat $out/$size/$num/*".DELs" >> $out/uncleanedResultFiles/oryLat04.DELsCombinedUnmerged_$size"bp"
			fi

			# numlines=$(find $out/$size/$num/ -name "*.smallGenicCONDELs" | wc -l)
			# if [[ $numlines -gt "0" ]]; then
			# 	cat $out/$size/$num/*".smallGenicCONDELs" >> $out/uncleanedResultFiles/oryLat04.smallGenicCONDELsCombinedUnmerged_$size"bp"
			# fi

			num=$(( num + 1 ))
		done
	done

	cat $out/uncleanedResultFiles/oryLat04.CONDELsCombinedUnmerged_* | sort -k1,1 -k2,2n > $out/uncleanedResultFiles/oryLat04.CONDELs_preIntactTargetFilter.bed
	subtractBed -A -a $out/uncleanedResultFiles/oryLat04.CONDELs_preIntactTargetFilter.bed -b $nonOrthoTargetAxtsBed > $out/oryLat04.CONDELs.bed
	
	mergeBed -d $mergeDist -i $out/oryLat04.CONDELs.bed | sort -k1,1 -k2,2n | awk '{print $0 "\t" $3-$2}' > $out/oryLat04.CONDELs.lengthsOnly.bed
	cut -f4 $out/oryLat04.CONDELs.bed | cut -f2 -d"_" | sort -u > $out/uniqGenes
 
	intersectBed -wa -u -a $out/oryLat04.CONDELs.lengthsOnly.bed -b $oryLat04_codingExons > $out/oryLat04.coding.CONDELs.bed
	intersectBed -wa -u -a $out/oryLat04.CONDELs.lengthsOnly.bed -b $oryLat04_5UTR $oryLat04_3UTR > $out/oryLat04.UTR.CONDELs.bed
	intersectBed -wa -u -a $out/oryLat04.CONDELs.lengthsOnly.bed -b $oryLat04_nonCodingExons > $out/oryLat04.nonCodingExons.CONDELs.bed
	intersectBed -wa -u -a $out/oryLat04.CONDELs.lengthsOnly.bed -b $oryLat04_introns > $out/oryLat04.introns.CONDELs.bed

	numlines=$(find $out/uncleanedResultFiles/ -name "oryLat04.DELsCombinedUnmerged_*" | wc -l)

	if [[ $numlines -gt 0 ]]; then
		cat $out/uncleanedResultFiles/oryLat04.DELsCombinedUnmerged_* | sort -k1,1 -k2,2n | mergeBed -i stdin  > $out/oryLat04.tempDELs.bed
		intersectBed -wa -u -a $out/oryLat04.tempDELs.bed -b $out/oryLat04.CONDELs.lengthsOnly.bed | awk 'OFS="\t" {print $0, "pDEL_"NR}' > $out/oryLat04.numberedDELs.bed
	# 	intersectBed -wa -u -a $out/oryLat04.DELs.bed -b $oryLat04_codingExons > $out/oryLat04.coding.DELs.bed
	#	intersectBed -wa -u -a $out/oryLat04.DELs.bed -b $oryLat04_5UTR $oryLat04_3UTR > $out/oryLat04.UTR.DELs.bed
	#	intersectBed -wa -u -a $out/oryLat04.DELs.bed -b $oryLat04_nonCodingExons > $out/oryLat04.nonCodingExons.DELs.bed
	#	intersectBed -wa -u -a $out/oryLat04.DELs.bed -b $oryLat04_introns > $out/oryLat04.introns.DELs.bed
	fi

	# numlines=$(find $out/uncleanedResultFiles/ -name "oryLat04.smallGenicCONDELsCombinedUnmerged_*" | wc -l)
	# if [[ $numlines -gt 0 ]]; then

	# 	cat $out/uncleanedResultFiles/oryLat04.smallGenicCONDELsCombinedUnmerged_* | sort -k1,1 -k2,2n > $out/oryLat04.smallGenicCONDELs.bed
	# fi
	wc -l $out/*bed > $out/stats
	wc -l $out/logFiles/*txt >> $out/stats
	wc -l $out/uniqGenes* >> $out/stats
	
	${scriptsDir}/annotateCONDELsWrapper.sh $out

else
	echo "WARNING: Some jobs did not finish"
fi

if [[ -s $out/oryLat04.coding.CONDELs.bed ]]; then
	rm $out/oryLat04.tempDELs.bed
	for window in $windowSizes; do
		rm -r $out/$window
	done
fi
