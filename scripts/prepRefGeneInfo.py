import sys
import gzip
import pickle
from pathlib import Path

'''
Takes in a gzipped file from ensembl biomart  with 8 columns:
1. Gene stable ID
2. Transcript stable ID
3. Chromosome/scaffold name
4. Transcription start site (TSS)
5. Gene name
6. Exon start
7. Exon end
8. Gene description
'''
biomartFilePath = sys.argv[1]
#biomartFilePath = "../ensembl98_ASM223467v1_geneInfo_originalFile.tab.gz"
script_dir = str(Path( __file__ ).parent.absolute())
pickOrthoChainsDir = script_dir + "/../pickOrthoChains/"

genes_chr = {}
genes_coords = {}

all_genes_transcripts = {} # key = gene_id | values = list of isoform transcript IDs
canonical_genes_transcript_ids = {} # gene_id | values = canonical transcript
all_transcripts_exon_coords = {} # key = transcript_id | values = list of exon coord tuples
all_TSSs = {} # key = transcript_id | values = [gene_id, chrom, tssStart, tssStart+1, name, desc]

count = 0
with gzip.open(biomartFilePath, 'rt') as geneInfo:
    for line in geneInfo:
        if "Chromosome/scaffold" in line:
            continue
        line = line.strip("\n")
        
        
        words = line.split("\t")
        gene_id = words[0]
        transcript_id = words[1]
        chrom = "chr"+words[2]
        tssStart = int(words[3]) - 1
        name = words[4]
        exon_start = int(words[5]) - 1
        exon_end = int(words[6])
        desc = words[7]
        
        if name == "":
            name = "[no_name]"
        if desc == "":
            desc = "[no_description]"
        
        if transcript_id not in all_TSSs:
            all_TSSs[transcript_id] = [chrom, str(tssStart), str(tssStart+1), name, gene_id, transcript_id, desc]
        
        exons = all_transcripts_exon_coords.get(transcript_id, [])
        exons.append((exon_start, exon_end))
        all_transcripts_exon_coords[transcript_id] = exons
        genes_chr[gene_id] = chrom
        transcripts = all_genes_transcripts.get(gene_id, [])
        if transcript_id not in transcripts:
            transcripts.append(transcript_id)
            all_genes_transcripts[gene_id] = transcripts

# get exons for longest (i.e. canonical) transcript
for gene_id in genes_chr.keys():
    max_bases = 0
    for transcript_id in all_genes_transcripts[gene_id]:
        bases = sum([c[1] - c[0] for c in all_transcripts_exon_coords[transcript_id]])
        if bases > max_bases:
            max_bases = bases
            genes_coords[gene_id] = all_transcripts_exon_coords[transcript_id]
            canonical_genes_transcript_ids[gene_id] = transcript_id
    count += len(genes_coords[gene_id])
    
with open(pickOrthoChainsDir+'genes_coords.p', 'wb') as genesCoords:
    pickle.dump(genes_coords, genesCoords)
with open(pickOrthoChainsDir+'genes_chr.p', 'wb') as genesChr:
    pickle.dump(genes_chr, genesChr)

with open(script_dir+"/../processedInputs/filterBEDs/ASM223467v1_ensembl98_canonicalTSS_withInfo.tab", "w") as canonicalTSSinfo:
    for transcript_id in sorted(canonical_genes_transcript_ids.values()):
        print("\t".join(all_TSSs[transcript_id]), file=canonicalTSSinfo)

#print(count)
#print(len(genes_chr.keys()))

#for k, v in canonical_genes_transcript_ids.items():
#    print(k, v)


