import sys
import pickle

genes_chr = {}
genes_coords = {}

genes_transcripts = {}
genes_transcript_ids = {}
transcripts_coords = {}


count = 0
for line in open(sys.argv[1]):
	if "#" not in line:
		words = line.split()
		gene_id = words[0]
		transcript_id = words[1]
		chrom = words[2]
		exon_start = int(words[3])
		exon_end = int(words[4])

		exons = transcripts_coords.get(transcript_id, [])
		exons.append((exon_start, exon_end))
		transcripts_coords[transcript_id] = exons
		genes_chr[gene_id] = chrom
		transcripts = genes_transcripts.get(gene_id, [])
		if transcript_id not in transcripts:
			transcripts.append(transcript_id)
			genes_transcripts[gene_id] = transcripts


# get exons for longest (i.e. canonical) transcript
for gene_id in genes_chr.keys():
	max_bases = 0
	for transcript_id in genes_transcripts[gene_id]:
		bases = sum([c[1] - c[0] for c in transcripts_coords[transcript_id]])
		if bases > max_bases:
			max_bases = bases
			genes_coords[gene_id] = transcripts_coords[transcript_id]
			genes_transcript_ids[gene_id] = transcript_id
	count += len(genes_coords[gene_id])

pickle.dump(genes_coords, open('genes_coords.p', 'w'))
pickle.dump(genes_chr, open('genes_chr.p', 'w'))

print count
print len(genes_chr.keys())

for k, v in genes_transcript_ids.items():
	print k, v



