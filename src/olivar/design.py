#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Functions for primer design
Architecture:
design.py
    primer_generator(object)
        WAlignScore()
    PrimerSetBadnessFast()
'''


__author__ = 'Michael X. Wang'


import os
CURRENT_DIR = os.path.dirname(__file__)


import random
import json
from collections import Counter

import numpy as np
import pandas as pd
from Bio import SeqIO

import basic
from basic import seq2num
from basic import num2seq
from basic import revcomp


GC_LETTER = ['C', 'G', 'c', 'g']


def WAlignScore(seq1, seq2=None):
    '''
    Calculate primer self badness
    '''
    PENALTY_OFFSET = -0.75 # Formulat is (PO + 1) / (PO + d+1)
    # if =1, then d=1 gives 2/3 socre, if =0, then d=1 gives 1/2 score, if =-0.5, then d=1 gives 1/3 score
    indel_penalty = 3
    replace_penalty = 2
    GC_score = 2
    AT_score = 1
    
    seq1 = seq1.lower()
    len1 = len(seq1)
    
    if seq2 is not None:
        seq2 = seq2.lower()
        len2 = len(seq2)
    else:
        seq2 = basic.revcomp(seq1)
        len2 = len1
        
    if len1 == 0:
        max_score = 0
    else:
        score = np.zeros([len1, len2])
        if seq1[0] == seq2[0]:
            if seq1[0] == 'c' or seq1[0] == 'g':
                score[0][0] = GC_score
            else:
                score[0][0] = AT_score
                
        for j in range(1, len2):
            if seq1[0] == seq2[j]:
                if seq1[0] == 'c' or seq1[0] == 'g':
                    score[0][j] = GC_score
                else:
                    score[0][j] = AT_score
            else:
                score[0][j] = max(score[0][j-1] - indel_penalty, 0)
                
        for i in range(1, len1):
            if seq1[i] == seq2[0]:
                if seq1[i] == 'c' or seq1[i] == 'g':
                    score[i][0] = GC_score
                else:
                    score[i][0] = AT_score
            else:
                score[i][0] = max(score[i-1][0] - indel_penalty, 0)
                
            for j in range(1, len2):
                if seq1[i] == seq2[j]:
                    if seq1[i] == 'c' or seq1[i] == 'g':
                        score[i][j] = GC_score + score[i-1][j-1]
                    else:
                        score[i][j] = AT_score + score[i-1][j-1]
                else:
                    score[i][j] = max([score[i-1][j] - indel_penalty, 
                                       score[i][j-1] - indel_penalty, 
                                       score[i-1][j-1] - replace_penalty, 0])
                    
        # weight by offset
        for i in range(len1):
            for j in range(len2):
                distance = len1 - i - 1
                if seq1[i] == seq2[j]:
                    score[i][j] -= 3*distance
                else:
                    score[i][j] = 0
                    
        max_score = 1.5**(score.max() - 6)
    return max_score


class primer_generator(object):
    '''
    Usage:
        g = primer_generator(temperature, salinity)
        retult_1 = g.get(seq1, dG_max=-11.8, min_GC=0.2, max_GC=0.75, min_complexity=0.4, max_len=36, check_SNP=True, prefix='')
        retult_2 = g.get(seq2, dG_max=-11.8, min_GC=0.1, max_GC=0.9, min_complexity=0.3, max_len=36, check_SNP=True)
    '''
    def __init__(self, temperature=60, salinity=0.18, end_pool_path=None):
        # Parameter = np.array([[11, -7.6, -21.3],
        #                     [12, -7.2, -20.4],
        #                     [13, -8.4, -22.4],
        #                     [14, -7.8, -21.0],
        #                     [21, -7.2, -21.3],
        #                     [22, -7.6, -21.3],
        #                     [23, -8.2, -22.2],
        #                     [24, -8.5, -22.7],
        #                     [31, -8.5, -22.7],
        #                     [32, -7.8, -21.0],
        #                     [33, -8.0, -19.9],
        #                     [34, -10.6, -27.2],
        #                     [41, -8.2, -22.2],
        #                     [42, -8.4, -22.4],
        #                     [43, -9.8, -24.4],
        #                     [44, -8.0, -19.9],])
        paraH = np.array([[-7.6, -7.2, -8.4, -7.8],
                          [-7.2, -7.6, -8.2, -8.5],
                          [-8.5, -7.8, -8.0, -10.6],
                          [-8.2, -8.4, -9.8, -8.0]])
        paraS = np.array([[-21.3, -20.4, -22.4, -21.0],
                          [-21.3, -21.3, -22.2, -22.7],
                          [-22.7, -21.0, -19.9, -27.2],
                          [-22.2, -22.4, -24.4, -19.9]])
        self.paraS = paraS + 0.368 * np.log(salinity)
        self.paraG = paraH - (temperature + 273.15) * self.paraS / 1000
        self.dG_init = 0.2 - (temperature + 273.15) * (-5.7)/1000

        # load pre-generated pool of primer end sequences
        if not end_pool_path is None:
            with open(end_pool_path) as f:
                end_pool = json.load(f)
            self.check_end = True
            self.l_end = len(end_pool[0])
            self.end_pool = [s.lower() for s in end_pool]
        else:
            self.check_end = False
            self.l_end = None
            self.end_pool = None
        
    def StacksDG(self, seq):
        rep = basic.seq2arr(seq)
        dG_stacks = 0
        for i, _ in enumerate(rep[:-1]):
            dG_stacks += self.paraG[rep[i]][rep[i+1]]
        return dG_stacks
    
    def get(self, seq, prefix='', dG_max=-11.8, min_GC=0.2, max_GC=0.75, 
            min_complexity=0.4, min_len=15, max_len=36, check_SNP=True, check_BLAST=False):
        '''
        Input:
            Sequence to generate primers from; restrictions (dG, GC, etc.)
        Output:
            DataFrame of primers; dictionary of failure cases
        '''
        INTR_SCORE = 2

        # input should be lowercase
        if seq == seq.upper(): # all uppercase
            print('input sequence should be in lowercase, with SNPs in uppercase. ')
            seq = seq.lower()
            if check_SNP:
                print('ignore check_SNP')
                check_SNP = False

        primer_seq = [] # primer sequence
        primer_arr = [] # array format
        primer_len = [] # primer length without prefix
        primer_dist = [] # distance to mother sequence end
        badness = []
        primer_dG = []
        BLAST_hits = []
        end_fail = 0 # number of primers failed end check
        SNP_fail = 0 # number of primers failed SNP check
        GC_fail = 0 # number of primers failed GC ratio check
        complexity_fail = 0 # number of primers failed complexity check
        for end_pos in range(len(seq)-1, 8, -1): # end_pos is 0-based
            # check end
            if self.check_end:
                if not seq[end_pos-self.l_end+1:end_pos+1] in self.end_pool:
                    end_fail += 1
                    continue
            # check SNP in last 5 bases
            if check_SNP:
                if seq[end_pos-4:end_pos+1] != seq[end_pos-4:end_pos+1].lower():
                    SNP_fail += 1
                    continue
            start_pos = end_pos - 1
            curr_dG = self.dG_init + self.StacksDG(seq[start_pos:end_pos+1])
            # satisfy minimum length
            while (start_pos > end_pos - min_len + 1) and (start_pos > 0):
                start_pos -= 1
                curr_dG += self.StacksDG(seq[start_pos:(start_pos+2)])
            # long enough to satisfy energy but not exceeding max_len
            while (curr_dG > dG_max) and (start_pos > end_pos - max_len + 1) and (start_pos > 0):
                start_pos -= 1
                curr_dG += self.StacksDG(seq[start_pos:(start_pos+2)])
            # check restrictions
            if (curr_dG <= dG_max) or (end_pos - start_pos + 1 == max_len):
                GC_ratio = basic.get_GC(seq[start_pos:end_pos+1])
                if min_GC <= GC_ratio <= max_GC:
                    complexity = basic.get_complexity(seq[start_pos:end_pos+1])
                    if complexity >= min_complexity:
                        s_raw = seq[start_pos:end_pos+1]
                        s = prefix + s_raw
                        primer_seq.append(s)
                        primer_arr.append(basic.seq2arr(s))
                        primer_len.append(end_pos - start_pos + 1)
                        primer_dist.append(len(seq) - end_pos - 1)
                        primer_dG.append(curr_dG)
                        if check_BLAST:
                            #count = BLAST_seq(s_raw, **BLAST_config)
                            if count > 0:
                                b = INTR_SCORE * WAlignScore(s) + 1000 * np.log10(count)**4
                            else:
                                b = INTR_SCORE * WAlignScore(s)
                        else:
                            count = -1
                            b = INTR_SCORE * WAlignScore(s)
                        BLAST_hits.append(count)
                        badness.append(b)
                    else:
                        complexity_fail += 1
                else:
                    GC_fail += 1
        primers = pd.DataFrame({'seq': primer_seq, 
                                'arr': primer_arr, 
                                'primer_len': primer_len, 
                                'dist': primer_dist, 
                                'badness': badness, 
                                'dG': primer_dG, 
                                'BLAST_hits': BLAST_hits})
        return primers, {'end_fail': end_fail, 'SNP_fail': SNP_fail, 
        'GC_fail': GC_fail, 'complexity_fail': complexity_fail}


def PrimerSetBadnessFast(all_fP: list, all_rP=[], existing=[]):
    '''
    Input:
        all_fP: list of strings
        all_rP: list of strings
        existing: list of strings
    Output:
        Badness: total Loss of the primer set (all fP and rP)
        Badness_component: Badness of each primer
    '''
    endhash4 = Counter()
    endhash5 = Counter()
    endhash6 = Counter()
    middlehash7 = Counter()
    middlehash8 = Counter()

    PENALTY_OFFSET = 0 # if=1, then d=1 gives 2/3 socre, if =0, then d=1 gives 1/2 score, if =-0.5, then d=1 gives 1/3 score

    END4 = 1
    END5 = 4
    END6 = 20
    MIDDLE7 = 100
    MIDDLE8 = 500

    # NOTE: Intramolecular checks done at GeneratePrimerCandidates to be efficient

    # set all sequences as lowercase
    all_fP = [p.lower() for p in all_fP]
    all_rP = [p.lower() for p in all_rP]
    existing = [p.lower() for p in existing]

    # Set up end hash tables (check end 4 nt to middle binding)
    for p in all_fP+all_rP+existing:
        endhash4[p[-4:]] += 1
        endhash5[p[-5:]] += 1
        endhash6[p[-6:]] += 1
        
    # Set up middle hash table
    # Middlehash penalizes based on closeness to 3' end, Badness = 2 / (distance to 3' + 2);
    # So absolute last 7 is worth 1, 1nt away is worth 0.67, 2nt away is worth 0.5, etc.
    for p in all_fP+all_rP+existing:
        l = len(p)
        for j in range(l-6):
            middlehash7[p[j:j+7]] += (PENALTY_OFFSET + 1) / (l - j - 6 + PENALTY_OFFSET)
        for j in range(l-7):
            middlehash8[p[j:j+8]] += (PENALTY_OFFSET + 1) / (l - j - 7 + PENALTY_OFFSET)

    # Run through each sequence's reverse complement to add up badness
    Badness = 0
    Badness_component = []
    for one_side_primer in [all_fP, all_rP]:
        one_side_badness = []
        for p in one_side_primer:
            p_badness = 0
            c = basic.revcomp(p)
            l = len(c)
            for j in range(l-3):
                k = c[j:j+4]
                try:
                    endscore4 = endhash4[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (endscore4 * END4 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)
                except KeyError:
                    pass
            for j in range(l-4):
                k = c[j:j+5]
                try:
                    endscore5 = endhash5[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (endscore5 * END5 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)
                except KeyError:
                    pass
            for j in range(l-5):
                k = c[j:j+6]
                try:
                    endscore6 = endhash6[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (endscore6 * END6 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)
                except KeyError:
                    pass
            for j in range(l-6):
                k = c[j:j+7]
                try:
                    midscore7 = middlehash7[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (midscore7 * MIDDLE7 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)
                except KeyError:
                    pass
            for j in range(l-7):
                k = c[j:j+8]
                try:
                    midscore8 = middlehash8[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (midscore8 * MIDDLE8 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)
                except KeyError:
                    pass
            Badness += p_badness
            one_side_badness.append(p_badness)
        Badness_component.append(one_side_badness)

    return Badness, Badness_component
