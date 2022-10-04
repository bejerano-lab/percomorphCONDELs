#!/bin/bash -e

set -beEu -o pipefail


scriptsDir=$(realpath $(dirname $0))
dirWithSpeciesWindows=${scriptsDir}/../processedInputs/windowPctIDs

out=${dirWithSpeciesWindows}/par/catSortPar_$(date +%Y%m%d_%H%M)
mkdir -p $out
#-----------------------#

cat << EOF > ${out}/catAndSort.sh

#/!/bin/bash -e

set -beEu -o pipefail

windowsFolder=\$(readlink -f \$1)

root=\$(echo \$windowsFolder | awk -F"/" '{print \$(NF-1)"_"\$NF}')

random="catAndSortWindows"\$(shuf -i1-10000 -n1)\$(shuf -i1-10000 -n1)\$(shuf -i1-10000 -n1)\$root
mkdir -p /tmp/\$random


cat \$windowsFolder/*bpWindowPctIDs > /tmp/\$random/\$root"bpWindowsWithPctIDcombined"

sort -k5,5nr /tmp/\$random/\$root"bpWindowsWithPctIDcombined" > \$windowsFolder/../\$root"bpWindowsSortedByPctID"

rm -rf /tmp/\$random

EOF

#-----------------------#

chmod +x ${out}/catAndSort.sh

cp $(readlink -f $0) ${out}/copyOfRunScript_$(basename $0 .sh)_$(date +%Y%m%d_%H%M)

windowSizes="10 25 50 100"

for step in $windowSizes; do
	for folder in $(find ${dirWithSpeciesWindows}/*/${step} -type d); do
		echo ${out}/catAndSort.sh $(readlink -f $folder) >> ${out}/catAndSortWindowsJoblist
	done
done


HOMEDIR=/cluster/u/hchen17

if test -d "$HOMEDIR"; then
	ssh hoxa "cd ${out} && para make catAndSortWindowsJoblist" 2>&1 | tee -a ${out}/parasolErrOut   # append and print stderr and stdout from para make to the error log file
else
	cat ${out}/catAndSortWindowsJoblist | parallel
fi
