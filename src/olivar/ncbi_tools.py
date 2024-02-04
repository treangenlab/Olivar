#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Wrappers for NCBI BLAST+
Architecture:
ncbi_tools.py
    ns_simulation()
        BLAST_batch_short()
'''


__author__ = 'Michael X. Wang'


import os
CURRENT_DIR = os.path.dirname(__file__)


import copy
from time import time

import pandas as pd
import subprocess
from Bio.Blast import NCBIXML
#from Bio.Blast.Applications import NcbiblastnCommandline
# The Bio.Application modules and modules relying on it have been deprecated
from tqdm import tqdm


def BLAST_batch_short(seq_list: list, db: str, n_cpu=1, seq_names=None, mode='rough'):
    '''
    Run local BLAST on a list of short sequences and return a corresponding list of hits. 
    Input:
        seq_list: list of sequences
        db: path to BLAST database
        n_cpu: number of threads to use
        seq_names: Default None. If specified, should be a list of sequence names corresponding to seq_list
        mode: 'precise' or 'rough'
    '''
    BATCH_SIZE = 128
    N = len(seq_list)

    BLAST_input = os.path.join(CURRENT_DIR, 'BLAST_query.fasta')
    if mode == 'precise':
        BLAST_output = os.path.join(CURRENT_DIR, 'BLAST_result.xml')
    elif mode == 'rough':
        BLAST_output = os.path.join(CURRENT_DIR, 'BLAST_result.tsv')
    else:
        raise ValueError('Invalid mode. mode should be either "precise" or "rough"')

    if seq_names is None:
        seq_names = ['query_%d' % i for i in range(N)]
    if N != len(seq_names):
        raise ValueError('length of seq_names should be the same as seq_list.')
    
    print(f'Running blastn-short with {n_cpu} CPU(s), batch size {BATCH_SIZE}...')
    tik = time()
    pbar = tqdm(total=N, ascii=' >') # progress bar

    n_batch = 0
    if mode == 'precise':
        all_count_batch = []
    elif mode == 'rough':
        all_count_batch = [0] * N
    all_loc_batch = []
    while BATCH_SIZE*n_batch < N:
        batch_start = BATCH_SIZE*n_batch
        batch_stop = min(BATCH_SIZE*(n_batch+1), N)
        seq_batch = seq_list[batch_start: batch_stop]
        name_batch = seq_names[batch_start: batch_stop]
        n_batch += 1
        #print('%d/%d' % (batch_start, N))

        # prepare fasta file for query
        with open(BLAST_input, 'w') as f:
            for name, seq in zip(name_batch, seq_batch):
                f.write('>%s\n%s\n' % (name, seq))
        
        # run blast
        if mode == 'precise':
            # referred to Primer-BLAST paper
            config = {
                'outfmt': 5, # xml
                'evalue': 5000, 
                'reward': 1, 
                'penalty': -1, 
                'gapopen': 2, 
                'gapextend': 1
            }
        elif mode == 'rough':
            # default of blastn-short
            config = {
                'outfmt': 6, # tabular
                'evalue': 10, 
                'reward': 1, 
                'penalty': -3, 
                'gapopen': 5, 
                'gapextend': 2
            }
        # blastn_cline = NcbiblastnCommandline(query=BLAST_input, 
        #     db=db, 
        #     out=BLAST_output, 
        #     task='blastn-short', 
        #     num_threads=n_cpu, 
        #     **config)
        # stdout, stderr = blastn_cline()
        cmd_result = subprocess.run([
            'blastn', 
            '-query', BLAST_input, 
            '-db', db, 
            '-out', BLAST_output, 
            '-task', 'blastn-short', 
            '-num_threads', str(n_cpu), 
            '-outfmt', str(config['outfmt']), 
            '-evalue', str(config['evalue']), 
            '-reward', str(config['reward']), 
            '-penalty', str(config['penalty']), 
            '-gapopen', str(config['gapopen']), 
            '-gapextend', str(config['gapextend'])
        ])
        
        # parse result file
        with open(BLAST_output, 'r') as handle:
            if mode == 'precise':
                all_records = NCBIXML.parse(handle)
                all_count = []
                all_loc = []
                for record in all_records:
                    # single BLAST query
                    count = 0
                    loc_dict = {}
                    qlen = record.query_letters
                    for alignment in record.alignments:
                        loc_list = []
                        for hsp in alignment.hsps:
                            if (hsp.query_end == qlen) and \
                            (hsp.match[-5:].count(' ')<2) and \
                            (hsp.identities>qlen-5):
                                # no 3' overhang
                                # no more than 1 mismatch within 5nt at 3'end
                                # no more than 4 total mismatches
                                loc_list.append((hsp.sbjct_start, hsp.sbjct_end, record.query))
                                count += 1
                        if loc_list:
                            loc_dict[alignment.accession] = loc_list
                    all_count.append(count)
                    all_loc.append(loc_dict)
                all_count_batch.extend(all_count)
                all_loc_batch.extend(all_loc)
            elif mode == 'rough':
                '''
                1. qseqid      query or source (e.g., gene) sequence id
                2. sseqid      subject  or target (e.g., reference genome) sequence id
                3. pident      percentage of identical matches
                4. length      alignment length (sequence overlap)
                5. mismatch    number of mismatches
                6. gapopen     number of gap openings
                7. qstart      start of alignment in query
                8. qend        end of alignment in query
                9. sstart      start of alignment in subject
                10. send        end of alignment in subject
                11. evalue      expect value
                12. bitscore    bit score
                '''
                # get first hsp
                hsp = handle.readline()[:-1].split('\t')
                query_name = hsp[0]
                count = 1
                for line in handle:
                    hsp = line[:-1].split('\t')
                    if hsp[0] == query_name:
                        count += 1
                    else:
                        # next query sequence
                        all_count_batch[int(query_name.split('_')[1])] = count
                        query_name = hsp[0]
                        count = 1
                all_count_batch[int(query_name.split('_')[1])] = count # save count of last record

        pbar.update(batch_stop-batch_start)
    pbar.close()

    print('finished in %.2fs' % (time()-tik))
    return all_count_batch, all_loc_batch


def ns_simulation(BLAST_db: str, seq_list: list, seq_names: list, max_amp_len=1500, n_cpu=1):
    '''
    Predict all possible non-specific amplification for primers in seq_list
    Input:
        BLAST_db: path to local BLAST database
        seq_list: list of primers
        seq_names: list of primer names
        max_amp_len: maximum amplicon length
        n_cpu: number of threads to use
    return: 
        dataframe of ns amplicons
        dataframe of ns primer pairs
    '''
    all_hits, all_loc = BLAST_batch_short(seq_list, db=BLAST_db, n_cpu=n_cpu, seq_names=seq_names, mode='precise')

    # merge hits within each alignment (chromosome)
    loc_merged = {}
    for loc_dict in all_loc:
        for acc, loc_list in loc_dict.items():
            try:
                loc_merged[acc].extend(copy.copy(loc_list))
            except KeyError:
                # copy to prevent modifying all_loc
                loc_merged[acc] = copy.copy(loc_list)

    # count legitimate primer pairs within each alignemnt (chromosome)
    pair_dict = {}
    for acc, loc_list in loc_merged.items():
        loc_list.sort(key=lambda x: x[0])
        list_len = len(loc_list)
        for i, loc in enumerate(loc_list[:-1]):
            j = i + 1
            while (j < list_len) and (loc_list[j][0] - loc[0] < max_amp_len):
                # amplicon length within range
                if (loc[0] < loc[1]) and (loc_list[j][1] < loc_list[j][0]):
                    # primer orientation is correct
                    pair = (loc[2], loc_list[j][2]) # primer names
                    try:
                        pair_dict[pair].append((acc, loc, loc_list[j]))
                    except KeyError:
                        pair_dict[pair] = [(acc, loc, loc_list[j])]
                j += 1
    
    # output each non-specific amplicon
    df_ns = {
        'fP': [], 
        'rP': [], 
        'acc': [], 
        'fP_start': [], 
        'fP_stop': [], 
        'rP_start': [], 
        'rP_stop': []
    }
    df_count = {
        'fP': [], 
        'rP': [], 
        'hits': []
    }
    for pair, loc_list in pair_dict.items():
        df_count['fP'].append(pair[0])
        df_count['rP'].append(pair[1])
        df_count['hits'].append(len(loc_list))
        for loc in loc_list:
            df_ns['fP'].append(pair[0])
            df_ns['rP'].append(pair[1])
            df_ns['acc'].append(loc[0])
            df_ns['fP_start'].append(loc[1][0])
            df_ns['fP_stop'].append(loc[1][1])
            df_ns['rP_start'].append(loc[2][0])
            df_ns['rP_stop'].append(loc[2][1])
    df_ns = pd.DataFrame(df_ns)
    df_ns['amp_len'] = df_ns['rP_start'] - df_ns['fP_start'] + 1
    return df_ns, pd.DataFrame(df_count), all_hits
