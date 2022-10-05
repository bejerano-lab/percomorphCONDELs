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

**3. Populate `processedInputs/chromSizes`**

> Generate assembly “chromosome sizes” files from each `2bit` using `twoBitInfo`
>
> *OR*
>
> Download and decompress a tarball of 36 pre-generated chrom.sizes files from [zenodo tbd]


**4. Populate `processedInputs/filterBEDs` and `pickOrthoChains` with files based on ensembl98 gene annotations for ASM223467v1**

> From [ensembl biomart](https://www.ensembl.org/biomart/martview), download a gzipped tab-delimted file with the following columns of information
> 1. Gene stable ID
> 2. Transcript stable ID
> 3. Chromosome/scaffold name
> 4. Transcription start site (TSS)
> 5. Gene name
> 6. Exon start
> 7. Exon end
> 8. Gene description
>
> Then run `python3 scripts/prepRefGeneInfo.py {pathToBiomartTabFile.gz}`
> 
> From [ensembl FTP](http://ftp.ensembl.org/pub/), download `Oryzias_latipes.ASM223467v1.98.gtf.gz` 
> 
> Then run 
> 
> *OR*
>
> Download and decompress a tarball of pre-generated files from [zenodo tbd]; place these files in `processedInputs/filterBEDs`
> Download and decompress a tarball of pre-generated files from [zenodo tbd]; place these files in `pickOrthoChains`

**5. Acquire whole-genome pairwise alignments to the reference assembly**

> *Skip to step 7 to use pre-generated orthologous alignment files*
> 
> Generate the alignments from scratch
> - Use doBlastzChainNet.pl, the parameters provided in `blastz_chain_DEF_files`, and chainLinearGap=medium
> - Place the resulting {reference}.{query}.all.chain.gz files in `pickOrthoChains/chains`
> 
> *OR*
>
> Download and decompress a tarball of pre-generated {reference}.{query}.all.chain.gz from [zenodo tbd]; place the files in `pickOrthoChains/chains`

**6. Map reference gene orthologs to identify orthologous alignment chains**

> For each query species, run `python2 scripts/pickChains_0.1threshold.py {query assembly name}`
> 
> e.g. `python2 scripts/pickChains_0.1threshold.py hipCom01`
> 
> This step may require 5-10g of RAM for each query assembly, as the {reference}.{query}.all.chain.gz alignments file will be decompressed into memory. It follows the method described in *Turakhia et al. Nucleic Acids Res. 2020 Sep 18;48(16):e91. https://doi.org/10.1093/nar/gkaa550.*

**7. Process orthologous chains**

> For each query species, run `scripts/processPickedChains.sh {query assembly name}`
> 
> e.g. `./scripts/processPickedChains.sh hipCom01`
> 
> *OR* (if skipping here from step 4)
> 
> Download and decompress a tarball of pre-generated {reference}.{query}.ortho.chains.gz from [zenodo tbd]; place the files in `processedInputs/orthoChains`
> Download and decompress a tarball of pre-generated orthologous gene mapping info from [zenodo tbd]; place the files in `pickOrthoChains`

**8. Generate lookup files for each reference gene**

> For each canonical reference gene, run `scripts/generateLookupFile.sh {ENSG ID}`


**9. Identify chain gaps**

**10. Identify intact alignment (axt) blocks**

**11. Identify sliding window conserved sequence elements**

*Extract chain ID-labeled alignment sequences*

*Compute percent sequence identity between the reference and query assemblies in sliding windows across the genome*

*Identify the most conserved windows*

**12. Run the screen**
