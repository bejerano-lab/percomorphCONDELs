import gzip
import numpy as np
import os.path
import math
import pickle
import sys

threshold = 0.1

min_chain_score = 10000
interval_len = 10000

args = sys.argv[1:]
query = args[0]
#ROOT='/cluster/u/hchen17/treeWAS/setup/percomorphs/ensembl98_ASM223467v1/'
ROOT= args[1]
filename = ROOT + 'chains/oryLat04.' + query + '.all.chain.gz'
outfile = open(ROOT + query + '.out', 'w')

f = gzip.open(filename, 'rb')
file_content = f.read()
lines = file_content.splitlines()
f.close()

print >> outfile,  'Reading chains'

chains = {} # dict.fromkeys(range(1, last_chain_id+1))
chain_blocks = {} # dict.fromkeys(range(1, last_chain_id+1), [])
chain_t_gaps = {} # dict.fromkeys(range(1, last_chain_id+1), [])
chain_q_gaps = {} # dict.fromkeys(range(1, last_chain_id+1), [])
chain_id = 0
chain_intervals = {}

assigned = {}
not_assigned = []
selected_chains = set([])

for line in lines:
    words = line.split()
    if (line[0:5] == 'chain'):
        chain_chr = words[2]
        if chain_chr not in chain_intervals.keys():
            chain_intervals[chain_chr] = [[]] * (1 + (int(words[3]) / interval_len))
        chain_start = int(words[5])
        chain_end = int(words[6])
        chain_q_chr = words[7]
        chain_q_chr_len = int(words[8])
        if chain_q_chr not in assigned.keys():
            assigned[chain_q_chr] = []
        chain_q_strand = words[9]
        chain_q_start = int(words[10])
        chain_q_end = int(words[11])
        chain_score = int(words[1])
        chain_id = int(words[12])
        if chain_score < min_chain_score:
#            print chain_id
            break
        chains[chain_id] = [chain_chr, chain_start, chain_end, chain_q_strand, chain_q_chr, chain_q_chr_len, chain_q_start, chain_q_end, chain_score]
        chain_blocks[chain_id] = []
        chain_t_gaps[chain_id] = []
        chain_q_gaps[chain_id] = []
        interval_start = (chain_start / interval_len)
        interval_end = 1 + (chain_end / interval_len)
        for i in range(interval_start,   interval_end):
            chain_intervals[chain_chr][i].append(chain_id)
    elif (chain_id > 0):
        if (len(words) > 0):
            block_length = int(words[0])
            chain_blocks[chain_id].append(block_length)
        if (len(words) > 1):
            t_gap_length = int(words[1])
            q_gap_length = int(words[2])
            chain_t_gaps[chain_id].append(t_gap_length)
            chain_q_gaps[chain_id].append(q_gap_length)

chain_lens = {}
for chain_id in chain_blocks.keys():
    chain_lens[chain_id] = sum(chain_blocks[chain_id])

print 'Done parsing chains.'
print 'Number of chains ', len(chain_t_gaps.keys())



def find_query_coords(chain_id, t_start, t_end, strand):
    chain_start = chains[chain_id][1]
    chain_end = chains[chain_id][2]
    if (chain_end < t_start) or (chain_start > t_end):
        return ['', 0, 0, 0]
    else:
        q_chr = chains[chain_id][4]
        t_curr = chain_start
        q_chr_len = chains[chain_id][5]
        q_curr = chains[chain_id][6]
        q_start = 0
        q_end = 0
        overlap_start = 0
        overlap_end = 0
        overlap = 0
        gaps_len = len(chain_t_gaps[chain_id])
        for (k, block) in enumerate(chain_blocks[chain_id]):
            overlap += max(0, min(t_curr + block, t_end) - max(t_curr, t_start))
#            if (overlap_start == 1) and (overlap_end == 0):
#                overlap += block
            if (q_start == 0) and (t_curr + block >= t_start):
                q_start = q_curr + max(0, (t_start - t_curr))
#                overlap += t_curr + block - t_start
#                overlap_start = 1
            if (q_end == 0) and (t_curr + block >= t_end):
                q_end = q_curr + max(0, t_end - t_curr)
#                overlap += t_end - t_curr
#                overlap_end = 1
            t_curr += block
            q_curr += block
            if (k < gaps_len):
                t_curr += chain_t_gaps[chain_id][k]
                q_curr += chain_q_gaps[chain_id][k]
        if (q_end == 0):
            q_end = q_curr
#        print [q_chr, q_start, q_end, overlap]
        if strand == '+':
            return [q_chr, q_start, q_end, overlap]
        else:
            return [q_chr, q_chr_len - q_end, q_chr_len - q_start, overlap]

def find_orthologous_gene(t_chr, t_exon_starts, t_exon_ends, gene_id):
    gene_len = sum([(t_exon_ends[k] - t_exon_starts[k]) for k in range(len(t_exon_starts))])
    best_chain_q_chr = ''
    best_chain_q_strand = ''
    best_chain_q_starts = [0] * len(t_exon_starts)
    best_chain_q_ends = [0] * len(t_exon_ends)
    exon_avail = [1] * len(t_exon_ends)
    best_chain_len = 1
    level = 1
    best_chain_level = 1
    highest_overlap = 1
    second_highest_overlap = 0
    highest_chain_score = 1
    second_highest_chain_score = 0
    overlapping_chains = []
    best_chain_id = -1
    second_best_chain_id = -1
    highest_overlap_score = 0
    level1_score = 1
    level1_chr = ''
    level1_overlap = 0
    level1_chain_id = -1
    level1_gene_in_synteny = 1
    level1_q_starts = [0] * len(t_exon_starts)
    level1_q_ends = [0] * len(t_exon_ends)
#    second_highest_chain_overlap_score = 0
    overlapping_exon_highest_chain_score = [0] * len(t_exon_starts)
    coding_start = t_exon_starts[0]
    coding_end = t_exon_ends[-1]
    if (coding_end <= coding_start):
        coding_start = t_exon_starts[-1]
        coding_end = t_exon_ends[0]
    chain_ids = []
    if t_chr in chain_intervals.keys():
        interval_start = (t_exon_starts[0] / interval_len)
        interval_end = 1 + (t_exon_ends[-1] / interval_len)
        if (interval_end <= interval_start):
            interval_start = (t_exon_starts[-1] / interval_len)
            interval_end = 1 + (t_exon_ends[0] / interval_len)
#        print (interval_start, interval_end)
        for i in range(interval_start, interval_end):
            if (i < len(chain_intervals[t_chr])):
                chain_ids.extend(chain_intervals[t_chr][i])
    chains_keys = sorted(set(chain_ids))
#    print (chains_keys)
    for chain_id in chains_keys:
        chain = chains[chain_id]
        chain_score = chain[-1]
        if (chain[0] == t_chr):
            chain_start = chain[1]
            chain_end = chain[2]
            chain_strand = chain[3]
            chain_overlap = min(coding_end, chain_end) - max(chain_start, coding_start)
            if (chain_overlap > 0):
                q_chr = chain[5]
                q_starts = []
                q_ends = []
                overlap_score = []
                overlaps = []
                for (k, exon_start) in enumerate(t_exon_starts):
                    exon_end = t_exon_ends[k]
                    coords = find_query_coords(chain_id, exon_start, exon_end, chain_strand)
                    overlap = coords[3]
                    if (exon_avail[k] > 0):
                        if (coords[3] > (t_exon_ends[k] - t_exon_starts[k])/4):
                            exon_avail[k] = 0
                    else:
                        overlap = 0
                    q_starts.append(coords[1])
                    q_ends.append(coords[2])
                    overlap_score.append(overlap)
                    overlaps.append(coords[3])
                total_overlap_score = sum(overlap_score)
                total_overlap = sum(overlaps)
                overlapping_chains.append([chain_id, total_overlap, chain_score])
#                chain_overlap_score = math.log10(chain_score) * total_overlap 
#                if (chain_overlap_score > highest_chain_overlap_score):
                if (total_overlap_score > highest_overlap_score):
                    highest_overlap_score = total_overlap_score
                    best_chain_id = chain_id
                    best_chain_q_chr = chain[4]
                    best_chain_q_strand = chain[3]
                    best_chain_q_starts = q_starts
                    best_chain_q_ends = q_ends
                    best_chain_level = level
                    highest_overlap = total_overlap
                    highest_chain_score = chain_score
                    best_chain_len = chain_end - chain_start
                if (chain_overlap == coding_end - coding_start):
                    if (level == 1):
                        level1_score = chain_score
                        level1_overlap = float(total_overlap) / gene_len
                        level1_chain_id = chain_id
                        level1_chr = chain[4]
                        level1_q_starts = q_starts
                        level1_q_ends = q_ends
                    level = level+1
    for overlapping_chain in overlapping_chains:
        chain_id = overlapping_chain[0]
        total_overlap = overlapping_chain[1]
        chain_score = overlapping_chain[2]
        if ((2*total_overlap > highest_overlap) and (chain_id != best_chain_id) and (chain_score > second_highest_chain_score)):
            second_highest_chain_score = chain_score
            second_highest_overlap = total_overlap
            second_best_chain_id = chain_id
    gene_in_synteny = float(gene_len) / chain_lens.get(best_chain_id, 1)
    level1_gene_in_synteny = float(gene_len) / chain_lens.get(level1_chain_id, 1)
    second_best_score_ratio = (float(second_highest_chain_score) / highest_chain_score) 
    second_best_overlap_ratio = (float(second_highest_overlap) / highest_overlap) 
    level1_score_ratio = float(highest_chain_score) / level1_score
    return (gene_id, best_chain_id, second_best_chain_id, best_chain_q_chr, best_chain_q_starts, best_chain_q_ends, best_chain_level, gene_in_synteny, second_best_score_ratio, second_best_overlap_ratio, best_chain_q_strand, level1_chain_id, level1_chr, level1_q_starts, level1_q_ends, level1_gene_in_synteny, level1_score_ratio, level1_overlap)

def is_assigned(best_chain_id, best_chain_q_chr, best_chain_q_starts, best_chain_q_ends, assigned):
    overlap = 0
    for coords in assigned[best_chain_q_chr]:
        for (k, q_start) in enumerate(best_chain_q_starts):
            q_end = best_chain_q_ends[k]
            overlap = max(0, min(coords[1], q_end) - max(coords[0], q_start))
            if ((overlap > 0) and (best_chain_id != coords[2])):
                print >> outfile,  overlap, best_chain_id, best_chain_q_chr, q_start, q_end, coords
                return True
    return False

#if 'ENSEMBL_ACCOUNT' in os.environ:
#    host, username, password = os.environ['ENSEMBL_ACCOUNT'].split()
#    account = HostAccount(host, username, password)
#else:
#    account = None
#
#human = Genome(Species=target, Release=Release, account=account)
#
#coding_genes = human.getFeatures(CoordName='19', feature_types='gene')

genes_chr = pickle.load(open(ROOT + 'genes_chr.p', 'rb'))
genes_coords = pickle.load(open(ROOT + 'genes_coords.p', 'rb'))

count = 0

for transcript_id in genes_chr.keys():
#  count c= 1
#  if (gene.BioType == 'protein_coding'):   
#     transcript = gene.CanonicalTranscript
#     transcript_id = transcript.StableId
#     count += 1
#     if count == 100:
#         break
     chr =  genes_chr[transcript_id]
     exon_starts = [genes_coords[transcript_id][k][0] for k in range(len(genes_coords[transcript_id]))]
     exon_ends = [genes_coords[transcript_id][k][1] for k in range(len(genes_coords[transcript_id]))]
#     for exon in transcript.Exons:
#         print >> outfile,  exon.Location.CoordName, exon.Location.Start, exon.Location.End
#         chr =  exon.Location.CoordName
#         exon_starts.append(int(exon.Location.Start))
#         exon_ends.append(int(exon.Location.End))
     (gene_id, best_chain_id, second_best_chain_id, best_chain_q_chr, best_chain_q_starts, best_chain_q_ends, best_chain_level, gene_in_synteny, second_best_score_ratio, second_best_overlap_ratio, best_chain_q_strand, level1_chain_id, level1_chr, level1_q_starts, level1_q_ends, level1_gene_in_synteny, level1_score_ratio, level1_overlap) = (find_orthologous_gene(chr, exon_starts, exon_ends, transcript_id))
     if ((gene_in_synteny < threshold) and (second_best_score_ratio < threshold) and (best_chain_level == 1)):
         print >> outfile,  gene_id, best_chain_id, second_best_chain_id, best_chain_q_chr, best_chain_q_starts, best_chain_q_ends, best_chain_level, gene_in_synteny, second_best_score_ratio, second_best_overlap_ratio, best_chain_q_strand, level1_chain_id, level1_chr, level1_score_ratio, level1_overlap
         print gene_id, best_chain_id, second_best_chain_id, best_chain_q_chr, best_chain_q_starts, best_chain_q_ends, best_chain_level, gene_in_synteny, second_best_score_ratio, second_best_overlap_ratio, best_chain_q_strand, level1_chain_id, level1_chr, level1_score_ratio, level1_overlap
         selected_chains.update(set([best_chain_id]))
         for (k, q_start) in enumerate(best_chain_q_starts):
             q_end = best_chain_q_ends[k]
             if (q_end > q_start):
                 assigned[best_chain_q_chr].append((q_start, q_end, best_chain_id))
             else: 
                 assigned[best_chain_q_chr].append((q_end, q_start, best_chain_id))
     else:
        not_assigned.append(transcript_id)

for transcript_id in not_assigned:
     chr =  genes_chr[transcript_id]
     exon_starts = [genes_coords[transcript_id][k][0] for k in range(len(genes_coords[transcript_id]))]
     exon_ends = [genes_coords[transcript_id][k][1] for k in range(len(genes_coords[transcript_id]))]
#     transcript = human.getTranscriptByStableId(transcript_id)
#     chr = ''
#     exon_starts = []
#     exon_ends = []
#     for exon in transcript.Exons:
#         print exon.Location.CoordName, exon.Location.Start, exon.Location.End
#         chr =  exon.Location.CoordName
#         exon_starts.append(int(exon.Location.Start))
#         exon_ends.append(int(exon.Location.End))
     (gene_id, best_chain_id, second_best_chain_id, best_chain_q_chr, best_chain_q_starts, best_chain_q_ends, best_chain_level, gene_in_synteny, second_best_score_ratio, second_best_overlap_ratio, best_chain_q_strand, level1_chain_id, level1_chr, level1_q_starts, level1_q_ends, level1_gene_in_synteny, level1_score_ratio, level1_overlap) = (find_orthologous_gene(chr, exon_starts, exon_ends, transcript_id)) 
     if (((gene_in_synteny < threshold) or ((best_chain_level==1) and(best_chain_id in selected_chains))) and (second_best_score_ratio < threshold) and (level1_score_ratio > 0.1 or level1_overlap <= 0.3)):
         if not(is_assigned(best_chain_id, best_chain_q_chr, best_chain_q_starts, best_chain_q_ends, assigned)):
                selected_chains.update(set([best_chain_id]))
                print >> outfile, gene_id, best_chain_id, second_best_chain_id, best_chain_q_chr, best_chain_q_starts, best_chain_q_ends, best_chain_level, gene_in_synteny, second_best_score_ratio, second_best_overlap_ratio, best_chain_q_strand, level1_chain_id, level1_chr, level1_score_ratio, level1_overlap
                for (k, q_start) in enumerate(best_chain_q_starts):
                    q_end = best_chain_q_ends[k]
                    if (q_end > q_start):
                        assigned[best_chain_q_chr].append((q_start, q_end, best_chain_id))
                    else: 
                        assigned[best_chain_q_chr].append((q_end, q_start, best_chain_id))
         else:
             print >> outfile,  'Already assigned: ', transcript_id
     elif ((level1_overlap > 0.3) and (level1_gene_in_synteny < threshold)):
         if not(is_assigned(level1_chain_id, level1_chr, level1_q_starts, level1_q_ends, assigned)):
                selected_chains.update(set([level1_chain_id]))
                print >> outfile,  gene_id, level1_chain_id, best_chain_id, level1_chr, level1_q_starts, level1_q_ends, 1, level1_gene_in_synteny, second_best_score_ratio, second_best_overlap_ratio, best_chain_q_strand, level1_chain_id, level1_chr, level1_score_ratio, level1_overlap
                for (k, q_start) in enumerate(level1_q_starts):
                    q_end = level1_q_ends[k]
                    if (q_end > q_start):
                        assigned[level1_chr].append((q_start, q_end, level1_chain_id))
                    else: 
                        assigned[level1_chr].append((q_end, q_start, level1_chain_id))
     else:
         print >> outfile,  'No match found: ', transcript_id, gene_in_synteny, second_best_score_ratio, level1_score_ratio
