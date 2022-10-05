#!/bin/bash

set -beEuo pipefail

scriptsDir=$(realpath $(dirname $0))

ensemblBiomartFile=$(readlink -f $1)
ensemblGTFgz=$(readlink -f $2)

python3 ${scriptsDir}/prepRefGeneInfo.py ${ensemblBiomartFile}

## identify introns ##

# first make a gff3 file of exons from ensembl's gtf

zcat ${ensemblGTFgz} | awk '$3 == "exon"' > /tmp/exonsOnly.gtf
gt gtf_to_gff3 -force -o /tmp/exonsOnly.gff3 /tmp/exonsOnly.gtf
gt gff3 -force yes  -addintrons yes -sort yes /tmp/exonsOnly.gff3 | sed 's/ /_/g' > /tmp/exonsAndIntrons.gff3 



