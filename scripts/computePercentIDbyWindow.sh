#!/bin/bash -e

set -beEu -o pipefail

scriptsDir=$(realpath $(dirname $0))
#echo $scriptsDir
axtDir=${scriptsDir}/../processedInputs/labeledAxts
out=${scriptsDir}/../processedInputs/windowPctIDs
threshold=0.0

windowSizes="10 25 50 100"

mkdir -p $out/par
cp $0 $out/par/copyOfRunScript_$(basename $0 .sh)_$(date +%Y%m%d_%H%M)

pythonScript=${scriptsDir}/percentIDscan.py


#-----------------------#

cat << EOF > $out/par/computePercentIDperWindow.sh # axt.gz windowSize outdir
#/!/bin/bash -e

set -beEu -o pipefail

axtFile=\$(readlink -f \$1)  # e.g. ausLim01_chain1000_.axt.gz
windowSize=\$2
outdir=\$(readlink -f \$3)
outname=\$(basename \$axtFile .axt.gz)\$windowSize"bpWindowPctIDs"  # e.g. ausLim01_chain1000_10bpWindowPctIDs

# don't sort -- will have to concatenate a bunch later and sort by pct id then
python3 $pythonScript \$axtFile \$windowSize $threshold \$outdir/\$outname

EOF

#-----------------------#

chmod +x $out/par/computePercentIDperWindow.sh

for axtSetFolder in ${axtDir}/*axts; do
	query=$(basename ${axtSetFolder} | cut -d"." -f2)
	for size in ${windowSizes}; do
		mkdir -p ${out}/${query}/${size}
		for axt in ${axtSetFolder}/*axt.gz; do
			echo $out/par/computePercentIDperWindow.sh $(readlink -f $axt) $size $out/$query/$size "{check out exists+" $out/$query/$size/$(basename $axt .axt.gz)$size"bpWindowPctIDs}">> $out/par/joblist
		done
	done
done


HOMEDIR=/cluster/u/hchen17

if test -d "$HOMEDIR"; then
	ssh hoxa "cd ${out}/par && para make joblist" 2>&1 | tee -a $out/parasolErrOut   # append and print stderr and stdout from para make to the error log file
else
	cat ${out}/par/joblist | parallel
fi

bash ${scriptsDir}/catAndSortWindowsByPctID.sh



