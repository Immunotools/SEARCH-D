import os
import sys
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

import numpy as np

class MotifMatrix:
    def __init__(self, fname):
        self.nucls = ['A', 'C', 'G', 'T']
        self.matrix = []
        lines = open(fname).readlines()
        for i in range(1, len(lines)):
            splits = lines[i].strip().split()
            self.matrix.append([float(s) for s in splits])
        self.motif_length = len(self.matrix[0])

    def ComputeProbability(self, seq, pos):
        prob = 1
        for i in range(self.motif_length):
            if seq[pos + i] not in self.nucls:
                return 0
            prob *= self.matrix[i][self.nucls.index(seq[pos + i])]
        return prob

    def ComputeNonamerProbabiilty(self, seq):
        prob = 1
        for i in range(self.motif_length):
            if seq[i] not in self.nucls:
                return 0
            prob *= self.matrix[i][self.nucls.index(seq[i])]
        return prob

    def MotifLength(self):
        return self.motif_length

def ComputeProbabilityList(seq, motif_matrix):
    prob_list = []
    for i in range(len(seq) - motif_matrix.MotifLength() + 1):
        prob_list.append(motif_matrix.ComputeProbability(seq, i))
    return prob_list

def FindRSSOccurrences(seq, rss_df):
    matches_dict = {'HEPTAMER' : [], 'NONAMER' : []}
    for i in range(len(rss_df)):
        matches_dict[rss_df['Type'][i]].extend([m.start() for m in re.finditer(rss_df['Sequence'][i], seq)])
#    for hpos in matches_dict['HEPTAMER']:
#        for npos in matches_dict['NONAMER']:
#            if npos - hpos == 30:
#                print npos, hpos
    return matches_dict

def GetThreshold(prob_list, quan):
    return np.percentile(prob_list, quan)

def ProbabilityIsGood(prob, ref_probs, perc):
    return prob >= np.percentile(ref_probs, perc)

class DGene:
    def __init__(self, seq, l7, l9, r7, r9, pos):
        self.seq = seq
        self.l7 = l7
        self.l9 = l9
        self.r7 = r7
        self.r9 = r9
        self.pos = pos

def HeptamerIsGood(h, all_hs):
    return h in all_hs

#############################################
def main(ighd_fasta, output_fasta):
    motif_matrix = 'motif_matrices' #sys.argv[3]
    quantile = 30 #float(sys.argv[4])
    min_length = 10
    max_length = 50

    left7_matrix = MotifMatrix(os.path.join(motif_matrix, 'ighd_left_heptamer_matrix.txt'))
    right7_matrix = MotifMatrix(os.path.join(motif_matrix, 'ighd_right_heptamer_matrix.txt'))

    left9_matrix = MotifMatrix(os.path.join(motif_matrix, 'ighd_left_nonamer_matrix.txt'))
    right9_matrix = MotifMatrix(os.path.join(motif_matrix, 'ighd_right_nonamer_matrix.txt'))

    left_rss_df = pd.read_csv(os.path.join(motif_matrix, 'ighd_left_rss.txt'), delim_whitespace = True)
    left7_seqs = set(left_rss_df.loc[left_rss_df['Type'] == 'HEPTAMER']['Sequence'])
    left7_probs = [left7_matrix.ComputeNonamerProbabiilty(s) for s in left7_seqs]
    left9_seqs = set(left_rss_df.loc[left_rss_df['Type'] == 'NONAMER']['Sequence'])
    left9_probs = [left9_matrix.ComputeNonamerProbabiilty(s) for s in left9_seqs]

    right_rss_df = pd.read_csv(os.path.join(motif_matrix, 'ighd_right_rss.txt'), delim_whitespace = True)
    right7_seqs = set(right_rss_df.loc[right_rss_df['Type'] == 'HEPTAMER']['Sequence'])
    right7_probs = [right7_matrix.ComputeNonamerProbabiilty(s) for s in right7_seqs]
    right9_seqs = set(right_rss_df.loc[right_rss_df['Type'] == 'NONAMER']['Sequence'])
    right9_probs = [right9_matrix.ComputeNonamerProbabiilty(s) for s in right9_seqs]

    ighd_locus = ''
    for r in SeqIO.parse(ighd_fasta, 'fasta'):
        if r.id.find('REVERSE') != -1:
            splits = r.id.split('|')
            if splits[5] == 'REVERSE:True':
                r.seq = r.seq.reverse_complement()
        ighd_locus = str(r.seq).upper()
    print('IGHD locus (' + str(len(ighd_locus)) + ' bp) was extracted from ' + ighd_fasta)

    # finding heptamer matches
    left7_positions = [pos for pos in range(len(ighd_locus) - 7 + 1) if HeptamerIsGood(ighd_locus[pos : pos + 7],
                                                                                       left7_seqs)]
    right7_positions = [pos for pos in range(len(ighd_locus) - 7 + 1) if HeptamerIsGood(ighd_locus[pos : pos + 7],
                                                                                        right7_seqs)]

    additional_pos_l = [pos.start() for pos in re.finditer('CAC\wGTG', ighd_locus)]
    additional_pos_r = [pos.start() for pos in re.finditer('CAC\wGTG', ighd_locus)]

    left7_positions.extend(additional_pos_l)
    left7_positions = set(left7_positions)

    right7_positions.extend(additional_pos_r)
    right7_positions = set(right7_positions)

    # finding good pairs of heptamers
    d_gene_candidates = []
    for l7_pos in left7_positions:
        for r7_pos in right7_positions:
            if r7_pos > l7_pos and r7_pos - l7_pos <= max_length + 7 and r7_pos - l7_pos >= min_length + 7:
                d_gene_candidates.append((l7_pos, r7_pos, ighd_locus[l7_pos + 7 : r7_pos]))
    print(str(len(d_gene_candidates)) + ' IGHD candidated were found')

    # collect nonamers
    d_genes = []
    for l7_pos, r7_pos, d_seq in d_gene_candidates:
        l9_pos = l7_pos - 12 - 9
        r9_pos = r7_pos + 12 + 7
        if l9_pos < 0 or r9_pos >= len(ighd_locus):
            continue
        l9_prob = left9_matrix.ComputeProbability(ighd_locus, l9_pos)
        r9_prob = right9_matrix.ComputeProbability(ighd_locus, r9_pos)
        l7_seq = ighd_locus[l7_pos : l7_pos + 7]
        r7_seq = ighd_locus[r7_pos : r7_pos + 7]
        gene_pos = l7_pos + 7
        if ProbabilityIsGood(l9_prob, left9_probs, quantile) and ProbabilityIsGood(r9_prob, right9_probs, quantile):
            d_genes.append(DGene(d_seq, l7_seq, ighd_locus[l9_pos : l9_pos + 9],
                                 r7_seq, ighd_locus[r9_pos : r9_pos + 9], gene_pos))
    print(str(len(d_gene_candidates) - len(d_genes)) + ' candidated were filtered as false-positives')

    fh = open(output_fasta, 'w')
    ind = 0
    for d in sorted(d_genes, key = lambda x : x.pos):
        fh.write('>INDEX:' + str(ind + 1) + '|IGHD' + str(ind + 1) + '|L7:' + d.l7 + '|L9:' + d.l9 + '|R7:' + d.r7 + '|R9:' + d.r9 + '\n')
        fh.write(d.seq + '\n')
        ind += 1
    fh.close()
    print(str(len(d_genes)) + ' IGHD genes were written to ' + output_fasta)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('python search_d.py IGHD_locus.fasta output.fasta')
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
