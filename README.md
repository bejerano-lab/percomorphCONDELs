# percomorphCONDELs

This repository provides the code for identifying genomic deletions that are associated with pelvic and caudal fin reduction in independent lineages of fish, as described in the following manuscript:

[coming soon]


## System requirements
- Install UCSC Genome Browser utils from http://hgdownload.soe.ucsc.edu/admin/exe/ to be available in your $PATH
- Configure available python packages according to the file `environment.yaml`

## Prep and run the screen

**1. Populate `processedInputs/2bits` with `2bit` files for all genome assemblies**

> Download assembly fasta files from the sources listed in Table __, perform repeat masking (if needed), and convert the fasta files to `2bit` format using `faToTwoBit`
>
> *OR*
>
> Download and decompress a tarball of 36 pre-generated `2bit` files from https://doi.org/10.5281/zenodo.7140839

**2. Populate `processedInputs/gapTracks` with gapTrack files for all genome assemblies**

> Generate assembly “gapTrack” files from each `2bit` using `twoBitInfo -nBed` and filter to only assembly gaps of 5 Ns or longer
>
> *OR*
>
> Download and decompress a tarball of 36 pre-generated gapTrack files from https://doi.org/10.5281/zenodo.7140839

**3. Populate `processedInputs/chromSizes`**

> Generate assembly “chromosome sizes” files from each `2bit` using `twoBitInfo`
>
> *OR*
>
> Download and decompress a tarball of 36 pre-generated chrom.sizes files from https://doi.org/10.5281/zenodo.7140839


**4. Populate `processedInputs/filterBEDs` and `pickOrthoChains` with files based on ensembl98 gene annotations for ASM223467v1**

> From [ensembl biomart](https://www.ensembl.org/biomart/martview), download a gzipped, tab-delimited file with the following columns of information
> 1. Gene stable ID
> 2. Transcript stable ID
> 3. Chromosome/scaffold name
> 4. Transcription start site (TSS)
> 5. Gene name
> 6. Exon start
> 7. Exon end
> 8. Gene description
> 
> From [ensembl FTP](http://ftp.ensembl.org/pub/), download `Oryzias_latipes.ASM223467v1.98.chr.gtf.gz` 
> 
> Then run 
> ```
> scripts/prepRefGenesWrapper.sh {pathToBiomartTabFile.gz} Oryzias_latipes.ASM223467v1.98.chr.gtf.gz
> ```
> 
> *OR*
>
> Download and decompress a tarball of pre-generated files from https://doi.org/10.5281/zenodo.7140839; place these files in `processedInputs/filterBEDs`
> Download and decompress a tarball of pre-generated files from https://doi.org/10.5281/zenodo.7140839; place these files in `pickOrthoChains`

**5. Acquire whole-genome pairwise alignments to the reference assembly**

> *Skip to step 7 to use pre-generated orthologous alignment files*
> 
> Generate the alignments from scratch
> - Use [`doBlastzChainNet.pl`](https://github.com/ENCODE-DCC/kentUtils/blob/master/src/hg/utils/automation/doBlastzChainNet.pl), the parameters provided in `blastz_chain_DEF_files`, and `chainLinearGap=medium`
> - Place the resulting `{reference}.{query}.all.chain.gz` files in `pickOrthoChains/chains`
> 
> *OR*
>
> Download and decompress a tarball of pre-generated `{reference}.{query}.all.chain.gz` files from https://doi.org/10.5281/zenodo.7140839; place the files in `pickOrthoChains/chains`

**6. Map reference gene orthologs to identify orthologous alignment chains**

> For each query species, run `python2 scripts/pickChains_0.1threshold.py {query assembly name}`
> 
> e.g. `python2 scripts/pickChains_0.1threshold.py hipCom01`
> 
> This step may require 5-10g of RAM for each query assembly, as the `{reference}.{query}.all.chain.gz` alignments file will be decompressed into memory. It follows the method described in *Turakhia et al. Nucleic Acids Res. 2020 Sep 18;48(16):e91. https://doi.org/10.1093/nar/gkaa550.*

**7. Process orthologous chains**

> For each query species, run `scripts/processPickedChains.sh {query assembly name}`
> 
> e.g. `./scripts/processPickedChains.sh hipCom01`
> 
> *OR* (if skipping here from step 4)
> 
> Download and decompress a tarball of pre-generated `{reference}.{query}.ortho.chains.gz` files from https://doi.org/10.5281/zenodo.7140839; place the files in `processedInputs/orthoChains`
> 
> Download and decompress a tarball of pre-generated orthologous gene mapping info from https://doi.org/10.5281/zenodo.7140839; place the files in `pickOrthoChains`

**8. Generate lookup files for each reference gene**

> For each canonical reference gene, run `scripts/generateLookupFile.sh {ENSG ID}`
> 
> For example,
> 
> ```
> for gene in $(cut -f5 processedInputs/filterBEDs/ASM223467v1_ensembl98_canonicalTSS_withInfo.tab); do
>   echo scripts/generateLookupFile.sh ${gene} >> joblist
> done
> 
> cat joblist | parallel
> ```
> 
**9. Identify chain gaps**

> For each set of chains in `processedInputs/orthoChains/`, run `scripts/getChainDelsMinusPaddedAssemblyGaps.sh oryLat04.{query}.ortho.chains.gz`
> 
> *OR*
> 
> Download and decompress a tarball of pre-generated files from https://doi.org/10.5281/zenodo.7140839 and place the files in `processedInputs/DELs`

**10. Identify intact alignment (axt) blocks**

> For each set of chains in `processedInputs/orthoChains/`, run `scripts/findIntactAxts.sh  oryLat04.{query}.ortho.chains.gz`
> 
> Then run `scripts/getAxtBlocksIntactIn2TetraSyngnath.sh screenDefinitionParameters`
> 
> *OR*
> 
> Download and decompress a tarball of pre-generated files from https://doi.org/10.5281/zenodo.7140839 and place the files in `processedInputs/axtBlocks`
> 
> Then run `scripts/getAxtBlocksIntactIn2TetraSyngnath.sh screenDefinitionParameters`

**11. Identify sliding window conserved sequence elements**

> *Extract chain ID-labeled fasta alignments:*
>
> For each `{query}.chain_ids` file in `pickOrthoChains/`, run `scripts/getChainLabeledSeqAxt.sh pickOrthoChains/{query}.chain_ids`
>
> *Compute percent sequence identity between the reference and query assemblies in sliding windows across the genome:*
> 
> If all previous steps have been performed, simply run `scripts/computePercentIDbyWindow.sh`
> 
> *Identify the most conserved windows:*
> 
> For each collated windows file `processedInputs/windowPctIDs/{query}/{query}_*bpWindowsSortedByPctID`, run `scripts/getWindowsCovering5pctOfREFandCollapseByChain.sh processedInputs/windowPctIDs/{query}/{query}_{size}bpWindowsSortedByPctID`
> 
> *OR*
> 
> Download and decompress a tarball of pre-generated files from https://doi.org/10.5281/zenodo.7140839 and place the files in `processedInputs/slidingWindowCons`

**12. Run the screen**

> Run `scripts/screenForPercomorphCONDELs.sh screenDefinitionParameters` 
> 
> Screen output will be saved in `results`
