# percomorphCONDELs

This repository provides the code for identifying genomic deletions that are associated with pelvic and caudal fin reduction in independent lineages of fish. Please refer to the following manuscript for more details:

[coming soon]


**1. Populate `processedInputs/2bits` with `2bit` files for all genome assemblies**

> Download assembly fasta files from the sources listed in Table __ and convert them to `2bit` format using `faToTwoBit`
>
> *OR*
>
> Download and decompress a tarball of 36 pre-generated `2bit` files from [zenodo tbd]

**2. Populate `processedInputs/gapTracks` with gapTrack files for all genome assemblies**

>Generate assembly “gapTrack” files from each `2bit` using `twoBitInfo -nBed` and filter to only assembly gaps of 5 Ns or longer
>
> *OR*
>
> Download and decompress a tarball of 36 pre-generated gapTrack files from [zenodo tbd]

