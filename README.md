# percomorphCONDELs

This repository provides the code for identifying genomic deletions that are associated with pelvic and caudal fin reduction in independent lineages of fish, as described in the following manuscript:

[coming soon]


**1. Populate `processedInputs/2bits` with `2bit` files for all genome assemblies**

> Download assembly fasta files from the sources listed in Table __ and convert them to `2bit` format using `faToTwoBit`
>
> *OR*
>
> Download and decompress a tarball of 36 pre-generated `2bit` files from [zenodo tbd]

**2. Populate `processedInputs/gapTracks` with gapTrack files for all genome assemblies**

> Generate assembly “gapTrack” files from each `2bit` using `twoBitInfo -nBed` and filter to only assembly gaps of 5 Ns or longer
>
> *OR*
>
> Download and decompress a tarball of 36 pre-generated gapTrack files from [zenodo tbd]

**3. Populate `processedInputs/filterBEDs` and `pickOrthoChains` with files based on ensembl98 gene annotations for ASM223467v1**

> Download ... using BioMart
>
> *OR*
>
> Download and decompress a tarball of pre-generated files from  [zenodo tbd]; place these files in `processedInputs/filterBEDs`
> Download and decompress a tarball of pre-generated files from  [zenodo tbd]; place these files in `pickOrthoChains`

**4. Acquire whole-genome pairwise alignments to the reference assembly**

> *Skip to step 6 to use pre-generated orthologous alignment files*
> 
> Generate the alignments from scratch
> - Use doBlastzChainNet.pl, the parameters provided in `blastz_chain_DEF_files`, and chainLinearGap=medium
> - Place the resulting {reference}.{query}.all.chain.gz files in `pickOrthoChains/chains`
> 
> *OR*
>
> Download and decompress a tarball of pre-generated {reference}.{query}.all.chain.gz from  [zenodo tbd]; place the files in `pickOrthoChains/chains`

**5. Map reference gene orthologs to identify orthologous alignment chains**

> For each query species, run `python2 scripts/pickChains_0.1threshold.py {query assembly name}`
> 
> e.g. `python2 scripts/pickChains_0.1threshold.py hipCom01`
> 
> This step may require 5-10g of RAM, as the {reference}.{query}.all.chain.gz alignments file will be decompressed into memory. It follows the method described in *Turakhia et al. Nucleic Acids Res. 2020 Sep 18;48(16):e91. https://doi.org/10.1093/nar/gkaa550.*

**6. Process orthologous chains**
>

**7. Generate lookup files for each reference gene**
>
