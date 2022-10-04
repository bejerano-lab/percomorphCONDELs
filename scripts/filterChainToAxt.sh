#!/bin/bash

function usage {
	cat << EOF

Example usage: $0 125 hg38.mm10.all.chain(.gz) hg38.2bit mm10.2bit

$0 chainID target.query.chain.gz target2bit query2bit

The script takes in chainId (=125) for chain file (hg38.mm10.all.chain[.gz]), the target (hg38.2bit) and query 2bit (mm10.2bit) files

EOF
	exit 1;
}

if [ $# -ne 4 ]; then
	usage;
fi

chainId=$1
chainFile=$2
target2Bit=$3
query2Bit=$4

chainFilter -id=$chainId $chainFile | chainToAxt stdin $target2Bit $query2Bit stdout
