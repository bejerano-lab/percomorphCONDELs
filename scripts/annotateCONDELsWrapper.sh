#!/bin/bash

set -beEuo pipefail

# annotateCONDELsWrapper.sh

SCRIPTS_DIR=$(realpath $(dirname $0))

resultsDir=$(readlink -f $1)

cp $0 ${resultsDir}/logFiles/
cp ${SCRIPTS_DIR}/annotateOneCONDELbed_v3.py  ${resultsDir}/logFiles/

numResults=$(wc -l ${resultsDir}/oryLat04.CONDELs.lengthsOnly.bed)
python3 ${SCRIPTS_DIR}/prepAnnotationIntermediateFileDir.py $resultsDir $numResults
python3 ${SCRIPTS_DIR}/prepPerCONDELbedFiles.py $resultsDir

# generate job list
mkdir -p ${resultsDir}/annotationIntermediateFiles/par
for file in $(find ${resultsDir}/annotationIntermediateFiles/*/*bed); do
    echo python3 ${SCRIPTS_DIR}/annotateOneCONDELbed_v3.py ${resultsDir} ${file} >> ${resultsDir}/annotationIntermediateFiles/par/jobList
done


HOMEDIR=/cluster/u/hchen17
if test -d "$HOMEDIR"; then
	ssh hoxa "cd ${resultsDir}/annotationIntermediateFiles/par && para make jobList" 2>&1 | tee -a ${resultsDir}/logFiles/parasolErrOut   # append and print stderr and stdout from para make to the error log file

else
	cat ${resultsDir}/annotationIntermediateFiles/par/jobList | parallel
fi

python3 $SCRIPTS_DIR/printAnnotationHeader.py > ${resultsDir}/perGeneCONDELannotations.tab

cat ${resultsDir}/annotationIntermediateFiles/*/*annot | sort -k1,1 -k2,2n >> ${resultsDir}/perGeneCONDELannotations.tab

if [[ -s ${resultsDir}/perGeneCONDELannotations.tab ]]; then
	rm -rf ${resultsDir}/annotationIntermediateFiles
fi

