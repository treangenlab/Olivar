#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Main workflow of Olivar tiling.
Architecture:
main.py
    build()
    tiling()
        design_context_seq()
            generate_context()
                find_min_loc()
        get_primer()
        optimize()
        to_df()
        tiling_save()
    validate()
'''


__author__ = 'Michael X. Wang'


import os
CURRENT_DIR = os.path.dirname(__file__)

import json
import pickle # needs python>=3.8
import random
import copy
import multiprocessing
from time import time
from math import ceil, floor
from copy import deepcopy

import pandas as pd
import numpy as np
from numpy.random import rand, default_rng
from tqdm import tqdm

from Bio import SeqIO

import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# internal modules
import basic
import design # SADDLE
from ncbi_tools import BLAST_batch_short, ns_simulation
from design import PrimerSetBadnessFast

N_POOLS = 2 # number of primer pool (not to be changed)

# algorithmic parameters
PRIMER_DESIGN_LEN = 40 # length of primer design region [40]
UPPER_GC = 0.75 # extreme GC upper bond [0.75]
LOWER_GC = 0.25 # extreme GC lower bond [0.75]
LOW_COMPLEXITY = 0.4 # low complexity lower bond [0.4]
CHOICE_RATE = 0.3 # bottom CHOICE_RATE lowest risk primer design regions are choosen from candidate region
RISK_TH = 0.1 # top RISK_TH highest risk primer design regions are considered as loss

REFEXT = '.olvr' # extension for Olivar reference file


def build(fasta_path: str, var_path: str, BLAST_db: str, out_path: str, title: str, threads: int):
    '''
    Build the Olivar reference file for tiled amplicon design
    Input:
        fasta_path: Path to the fasta reference sequence.
        var_path: Optional, path to the csv file of SNP coordinates and frequencies. 
            Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.
        BLAST_db: Optional, path to the BLAST database. 
            Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        out_path: Output directory [./]. 
        title: Name of the Olivar reference file [olivar-ref]. 
        threads: Number of threads [1]. 
    '''
    n_cpu = threads

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # all coordinates are 1-based and inclusive
    word_size = 28
    offset = 14 # word_size should be divisible by offset
    print('Building Olivar reference with word_size=%d, offset=%d' % (word_size, offset))

    n_cycle = word_size/offset
    if n_cycle != int(n_cycle):
        raise ValueError('word_size should be divisible by offset.')
    else:
        n_cycle = int(n_cycle)

    for record in SeqIO.parse(fasta_path, 'fasta'):
        seq_raw = str(record.seq).lower() # SNPs in upper case
        break # read first record
    
    # check ambiguous bases
    ambiguous_bases = set(seq_raw) - {'a', 't', 'c', 'g'}
    if ambiguous_bases:
        raise ValueError(f'Ambiguous bases are not supported. {ambiguous_bases} found in reference fasta file.')

    # load SNP/iSNV coordinates
    var_arr = np.zeros(len(seq_raw))
    if var_path:
        seq_raw = list(seq_raw)
        df = pd.read_csv(var_path, sep=',', index_col=False)

        # variant csv file sanity check
        required_col = {
            'START': 'int64', 
            'STOP': 'int64'
        }
        optional_col = {
            'FREQ': 'float64'
        }
        for cname, expected_dtype in required_col.items():
            try:
                received_dtype = df[cname].dtype
            except KeyError:
                raise ValueError(f'Column {cname} not found in {var_path}')
            if received_dtype != expected_dtype:
                raise ValueError(f'Invalid value found in column {cname} of {var_path}, only {expected_dtype} is accepted (received {received_dtype})')
        for cname, expected_dtype in optional_col.items():
            try:
                received_dtype = df[cname].dtype
            except KeyError:
                # column not provided
                print("Warning: 'FREQ' is not provided")
                df[cname] = np.nan
                continue
            if received_dtype != expected_dtype:
                raise ValueError(f'Invalid value found in column {cname} of {var_path}, only {expected_dtype} is accepted (received {received_dtype})')

        # add variant frequency to var_arr
        for i, row in df.iterrows():
            freq = row['FREQ']
            if np.isnan(freq):
                # FREQ is not provided
                freq = 1.0
            var_arr[int(row['START'])-1:int(row['STOP'])] += freq
            # capitalize SNPs
            for pos in range(int(row['START'])-1, int(row['STOP'])):
                seq_raw[pos] = seq_raw[pos].upper()
        var_arr = var_arr**0.5 # amplify low frequency variants
        seq_raw = ''.join(seq_raw)
    else:
        print('No sequence variation coordinates provided, skipped.')

    # calculate risk of GC, non-specificity, complexity
    seq_len_temp = (len(seq_raw)//offset) * offset # length of temporary sequence
    start_temp = (len(seq_raw) - seq_len_temp)//2 + 1 # start coordinate of sequence to be processed as words
    stop_temp = start_temp + seq_len_temp - 1

    # fetch all words
    all_word = []
    for pos in range(0, seq_len_temp-word_size+1, offset):
        word_start = start_temp+pos-1
        word = seq_raw[word_start : word_start+word_size]
        all_word.append(word)

    # get GC, complexity, BLAST hits of each word
    with multiprocessing.Pool(processes=n_cpu) as pool:
        print('Calculating GC content and sequence complexity...')
        tik = time()
        all_gc = np.array(pool.map(basic.get_GC, all_word))
        all_complexity = np.array(pool.map(basic.get_complexity, all_word))
        print('Finished in %.3fs' % (time()-tik))
    if BLAST_db:
        print('Calculating non-specificity...')
        all_hits, _ = BLAST_batch_short(all_word, db=BLAST_db, n_cpu=n_cpu, mode='rough') # tabular output
    else:
        print('No BLAST database provided, skipped.')
        all_hits = np.zeros(len(all_word))

    # calculate score array
    start = start_temp + (n_cycle-1)*offset # start coordinate of sequence with score array
    stop = stop_temp - (n_cycle-1)*offset
    seq_len = stop - start + 1 # length of final sequence
    gc_arr = np.zeros(seq_len)
    comp_arr = np.zeros(seq_len)
    hits_arr = np.zeros(seq_len)
    for pos in range(seq_len):
        n = pos//offset
        gc_arr[pos] = np.sum(all_gc[n:n+n_cycle])
        comp_arr[pos] = np.sum(all_complexity[n:n+n_cycle])
        hits_arr[pos] = np.sum(all_hits[n:n+n_cycle])
    gc_arr = gc_arr/n_cycle
    hits_arr = hits_arr/n_cycle

    olv_ref = {
        'seq': seq_raw, 
        'seq_record': record, 
        'start': start, 
        'stop': stop, 
        'gc_arr': gc_arr, 
        'comp_arr': comp_arr, 
        'var_arr': var_arr[start-1:stop], 
        'hits_arr': hits_arr, 
        'all_hits': all_hits
    }
    save_path = os.path.join(out_path, f'{title}{REFEXT}')
    with open(save_path, 'wb') as f:
        pickle.dump(olv_ref, f, protocol=5) # protocol 5 needs python>=3.8
        print('Reference file saved as %s' % save_path)


def find_min_loc(risk_arr, start, stop, rng):
    '''
    find the next primer design region within start and stop.
    '''
    score = []
    loc = np.array(range(start, stop-PRIMER_DESIGN_LEN+2))
    for i in loc:
        score.append(sum(risk_arr[i-1: i-1+PRIMER_DESIGN_LEN]))
    if len(loc) == 1:
        idx_rand = 0
    elif CHOICE_RATE < 1:
        k = int(len(loc)*CHOICE_RATE)
        if k == 0:
            k = 1
        idx_k_lowest = np.argpartition(score, k)[:k]
        idx_rand = rng.choice(idx_k_lowest)
    else:
        idx_rand = rng.choice(range(len(loc)))
    primer_start = loc[idx_rand]
    return primer_start, primer_start+PRIMER_DESIGN_LEN-1, score[idx_rand]


def generate_context(SOI):
    '''
    find a valid set of PDRs that cover the whole sequence of interest (SOI).
    '''
    risk_arr, start, stop, max_amp_len, min_amp_len, seed = SOI
    SOI_len = stop - start + 1
    rng = default_rng(seed)
    
    all_context_seq = [] # coordinates of all context sequences
    all_risk = []

    # generate the first pair of primer design region
    fp_start, fp_stop, fp_risk = find_min_loc(risk_arr, start, start+3*PRIMER_DESIGN_LEN-1, rng)
    rp_stop, rp_start, rp_risk = find_min_loc(risk_arr, fp_start+min_amp_len-PRIMER_DESIGN_LEN, fp_start+max_amp_len-1, rng)
    all_context_seq.append((fp_start, fp_stop, rp_stop, rp_start))
    all_risk.append((fp_risk, rp_risk))

    # generate the second pair of primer design region
    left_lim = max(fp_stop+1, rp_start+2*PRIMER_DESIGN_LEN-max_amp_len+1) # leave space for the next fP design region
    fp_start, fp_stop, fp_risk = find_min_loc(risk_arr, left_lim, rp_stop-1, rng)
    left_lim = max(fp_start+min_amp_len-PRIMER_DESIGN_LEN, rp_start+PRIMER_DESIGN_LEN+1)
    right_lim = fp_start+max_amp_len-1
    if right_lim > stop:
        raise ValueError('sequence provided is too short.')
    rp_stop, rp_start, rp_risk = find_min_loc(risk_arr, left_lim, right_lim, rng)
    all_context_seq.append((fp_start, fp_stop, rp_stop, rp_start))
    all_risk.append((fp_risk, rp_risk))

    # design the rest context sequences
    cs_id = 2 # id of the next context sequence to be designed (0-based)
    while True:
        # fp
        # all_context_seq[cs_id-2][3]: rp_start of the last last context sequence
        left_lim = max(all_context_seq[cs_id-2][3]+1, rp_start+2*PRIMER_DESIGN_LEN-max_amp_len+1)
        fp_start, fp_stop, fp_risk = find_min_loc(risk_arr, left_lim, rp_stop-1, rng)
        # rp
        left_lim = max(fp_start+min_amp_len-PRIMER_DESIGN_LEN, rp_start+PRIMER_DESIGN_LEN+1)
        right_lim = fp_start+max_amp_len-1
        if right_lim > stop:
            break
        rp_stop, rp_start, rp_risk = find_min_loc(risk_arr, left_lim, right_lim, rng)

        all_context_seq.append((fp_start, fp_stop, rp_stop, rp_start))
        all_risk.append((fp_risk, rp_risk))
        cs_id += 1

    all_context_seq = np.array(all_context_seq)
    all_risk = np.array(all_risk)
    all_risk_flatten = all_risk.flatten()
    k = int(len(all_risk_flatten) * (1-RISK_TH))

    cover_start = all_context_seq[0][1]+1
    cover_stop = all_context_seq[-1][2]-1
    coverage = (cover_stop-cover_start+1)/SOI_len

    loss = sum(np.partition(all_risk_flatten, k)[k:]**2)/(coverage**2)
    #loss = sum(all_risk_flatten[all_risk_flatten>=RISK_TH])
    return all_context_seq, all_risk, loss


def design_context_seq(config):
    '''
    PDR optimization. 
    Input:
        config: design configurations
    Output:
        all_plex_info
    '''
    ref_path = config['ref_path']
    max_amp_len = config['max_amp_len']
    w_gc = config['w_egc']
    w_complexity = config['w_lc']
    w_hits = config['w_ns']
    w_var = config['w_var']
    seed = config['seed']
    n_cpu = config['threads']
    
    # set random number generator
    rng_parent = default_rng(seed)
    
    with open(ref_path,  'rb') as f:
        olv_ref = pickle.load(f)

    seq_raw = olv_ref['seq'] # SNPs are capitalized in build()
    start = olv_ref['start'] # start of design region, 1-based, closed
    stop = olv_ref['stop'] # stop of design region, 1-based, closed
    gc_arr = olv_ref['gc_arr']
    comp_arr = olv_ref['comp_arr']
    var_arr = olv_ref['var_arr']
    hits_arr = olv_ref['hits_arr']
    print('Successfully loaded reference file %s' % ref_path)
    
    if not config['min_amp_len']:
        min_amp_len = int(max_amp_len*0.9)
    else:
        min_amp_len = config['min_amp_len']
    if min_amp_len < 3*PRIMER_DESIGN_LEN:
        raise ValueError('min_amp_len is too small')

    gc_arr = np.logical_or(gc_arr<LOWER_GC, gc_arr>UPPER_GC).astype(int)
    comp_arr = (comp_arr < LOW_COMPLEXITY).astype(int)
    if max(hits_arr) != 0:
        hits_arr = hits_arr/max(hits_arr) # normalize

    # make arrays the same length as seq_raw
    gc_arr = w_gc*0.5*np.concatenate((np.zeros(start-1), gc_arr, np.zeros(len(seq_raw)-stop)))
    comp_arr = w_complexity*0.5*np.concatenate((np.zeros(start-1), comp_arr, np.zeros(len(seq_raw)-stop)))
    hits_arr = w_hits*np.concatenate((np.zeros(start-1), hits_arr, np.zeros(len(seq_raw)-stop)))
    var_arr = w_var*10*np.concatenate((np.zeros(start-1), var_arr, np.zeros(len(seq_raw)-stop)))

    # construct risk array (first row is risk, second row is coordinate on seq_rawy)
    risk_arr = gc_arr + comp_arr + hits_arr + var_arr

    N = 500*len(risk_arr)//max_amp_len # number of primer sets to generate
    #N = 100
    rand_int = rng_parent.integers(2**32, size=N) # random seeds for each iteration

    # single thread
    # design = [generate_context((risk_arr, start, stop, max_amp_len, min_amp_len, rand_int[i])) for i in tqdm(range(N))]

    # multi threads
    print('reference sequence length: %d' % len(seq_raw))
    print('design region: %d:%d' % (start, stop))
    print('region length: %d' % (stop-start+1))
    print('number of PDR sets to be tested: %d' % N)
    print('Designing PDRs...')
    tik = time()
    with multiprocessing.Pool(processes=n_cpu) as pool:
        batch = [(risk_arr, start, stop, max_amp_len, min_amp_len, rand_int[i]) for i in range(N)]
        design = list(
            tqdm(
                pool.imap(generate_context, batch, chunksize=max(N//n_cpu//10, 1)), # chunksize=1 is slow
                total=N, 
                ascii=' >'
            )
        )
    print('Finished in %.3fs' % (time()-tik))

    # Loss of all iterations
    all_loss = [d[2] for d in design]
    all_loss.sort(reverse=True)

    # find the best arangement of primer design regions
    best_design = sorted(design, key=lambda x:x[2])
    all_context_seq, all_risk, loss = best_design[0]
    print('Loss of the best design: %.3f' % loss)
    
    # prepare output
    all_plex_info = {}
    for i, (context_seq, risk) in enumerate(zip(all_context_seq, all_risk)):
        cs = seq_raw[context_seq[0]-1:context_seq[3]] # actual context sequence
        plex_info = {
            'reference': config['title'], 
            'tube': i%2 + 1, 
            'context_seq_coords': (context_seq[0], context_seq[3]), 
            'primers_coords': tuple(context_seq), 
            'insert_coords': (context_seq[1]+1, context_seq[2]-1), 
            'context_seq': cs, 
            'context_seq_revcomp': basic.revcomp(cs), 
            'fP_design': seq_raw[context_seq[0]-1:context_seq[1]], 
            'rP_design': basic.revcomp(seq_raw[context_seq[2]-1:context_seq[3]]), 
            'risk': tuple(risk)
        }
        all_plex_info['%s_%d' % (config['title'], i+1)] = plex_info
    print('total amplicons: %d' % (i+1))
    cover_start = all_context_seq[0][1]+1
    cover_stop = all_context_seq[-1][2]-1
    print('covered region: %d:%d' % (cover_start, cover_stop))
    print('coverage of reference sequence: %.3f%%' % (100*(cover_stop-cover_start+1)/len(seq_raw)))
    return all_plex_info, risk_arr, gc_arr, comp_arr, hits_arr, var_arr, all_loss, olv_ref['seq_record']


def get_primer(all_plex_info, config):
    '''
    Generate candidates of fP and rP for each plex
    Input:
        all_plex_info: output of design_context_seq()
        config: design configurations
    Output:
        all_plex_info: add primer candidates to input all_plex_info
    '''
    print('Generating primer candidates...')
    tik = time()

    temperature_fP = config['temperature']
    temperature_rP = config['temperature']
    salinity = config['salinity']
    fP_adapter = config['fP_prefix']
    rP_adapter = config['rP_prefix']
    default_fP_setting = {
        'dG_max': config['dG_max'], 
        'min_GC': config['min_GC'], 
        'max_GC': config['max_GC'], 
        'min_complexity': config['min_complexity'], 
        'max_len': config['max_len'], 
        'check_SNP': config['check_SNP']
    }
    default_rP_setting = deepcopy(default_fP_setting)
    check_BLAST = False

    # generate primer candidates
    fP_generator = design.primer_generator(temperature_fP, salinity)
    rP_generator = design.primer_generator(temperature_rP, salinity)
    n = 0
    for plex_id, plex_info in tqdm(all_plex_info.items(), ascii=' >'):
        # set primer prefix
        fP_prefix = fP_adapter
        rP_prefix = rP_adapter
        n += 1
        
        # generate fP
        fP_design = plex_info['fP_design'] # sequence to generate fP candidates
        fP_setting = deepcopy(default_fP_setting)
        fP, fail = fP_generator.get(fP_design, fP_prefix, check_BLAST=check_BLAST, **fP_setting)
        while len(fP) == 0:
            print('Fail to generate fP, update setting: %s' % plex_id)
            # detect failure mode with the most cases
            fail_mode = max(fail, key=fail.get)
            # update setting
            if fail_mode == 'SNP_fail' and fP_setting['check_SNP'] == True:
                fP_setting['check_SNP'] = False
            elif fail_mode == 'GC_fail':
                fP_setting['min_GC'] -= 0.05
                fP_setting['max_GC'] += 0.05
            elif fail_mode == 'complexity_fail':
                fP_setting['min_complexity'] -= 0.05
            # generate primers again
            fP, fail = fP_generator.get(fP_design, fP_prefix, check_BLAST=check_BLAST, **fP_setting)
        plex_info['fP_candidate'] = fP
        plex_info['fP_setting'] = fP_setting
        
        # generate rP
        rP_design = plex_info['rP_design'] # sequence to generate rP candidates
        rP_setting = deepcopy(default_rP_setting)
        rP, fail = rP_generator.get(rP_design, rP_prefix, check_BLAST=check_BLAST, **rP_setting)
        while len(rP) == 0:
            print('Fail to generate rP, update setting: %s' % plex_id)
            # detect failure mode with the most cases
            fail_mode = max(fail, key=fail.get)
            # update setting
            if fail_mode == 'SNP_fail' and rP_setting['check_SNP'] == True:
                rP_setting['check_SNP'] = False
            elif fail_mode == 'GC_fail':
                rP_setting['min_GC'] -= 0.05
                rP_setting['max_GC'] += 0.05
            elif fail_mode == 'complexity_fail':
                rP_setting['min_complexity'] -= 0.05
            # generate primers again
            rP, fail = rP_generator.get(rP_design, rP_prefix, check_BLAST=check_BLAST, **rP_setting)
        plex_info['rP_candidate'] = rP
        plex_info['rP_setting'] = rP_setting
        
    print('Finished in %.3fs' % (time()-tik))
    return all_plex_info


def optimize(all_plex_info, config):
    '''
    Optimize fP and rP with the SADDLE algorithm.
    Input:
        all_plex_info: output of get_primer()
        config: design configurations
    Output:
        all_plex_info: add optimized fP and rP to each plex of input all_plex_info
    '''
    print('Optimizing primer dimers...')
    # set random seed
    np.random.seed(config['seed'])
    random.seed(config['seed'])

    # automatically determine # SADDLE iterations
    pool_size = max(len(all_plex_info)/N_POOLS - 100, 0)
    InitSATemp = 1000 + 10*pool_size
    NUMSTEPS = 10 + int(pool_size/10)
    ZEROSTEPS = NUMSTEPS
    TimePerStep = 1000
    # InitSATemp = config['InitSATemp']
    # NUMSTEPS = config['NumSteps']
    # ZEROSTEPS = config['ZeroSteps']
    # TimePerStep = config['TimePerStep']

    # generate primer pair candidates (fP-rP)
    n_pair = []
    for plex_id, plex_info in all_plex_info.items():
        optimize = []
        plex_amp_len = []
        for i, fp in plex_info['fP_candidate'].iterrows():
            for j, rp in plex_info['rP_candidate'].iterrows():
                insert_len = plex_info['insert_coords'][1] - plex_info['insert_coords'][0] + 1
                # check amplicon length
                amp_len = fp['primer_len'] + fp['dist'] + insert_len + rp['dist'] + rp['primer_len']
                plex_amp_len.append(amp_len)
                optimize.append([fp, rp])
        plex_info['optimize'] = optimize
        n_pair.append(len(optimize))
    print('total primer pairs %d' % sum(n_pair))
    print('average pairs per plex %.2f' % (sum(n_pair)/len(n_pair)))
    
    # optimize each tube
    all_lc = []
    tik = time()
    for i_tube in range(1, N_POOLS+1):
        print('\npool %d, simulated annealing...' % i_tube)
        # plex_id of current tube
        curr_tube = [plex_id for plex_id, plex_info in all_plex_info.items() if plex_info['tube'] == i_tube]

        # load existing primers
        #existing_primer = config['existing_primer']
        existing_primer = []

        # randomly select one primer set (one pair each plex)
        curr_index = {}
        curr_fp = {}
        curr_rp = {}
        curr_badness = 0
        for plex_id in curr_tube:
            # randomly select one primer pair
            i = floor(rand() * len(all_plex_info[plex_id]['optimize']))
            curr_index[plex_id] = i
            fp = all_plex_info[plex_id]['optimize'][i][0]
            rp = all_plex_info[plex_id]['optimize'][i][1]
            curr_fp[plex_id] = fp['seq']
            curr_rp[plex_id] = rp['seq']
            curr_badness += fp['badness'] + rp['badness']
        inter_badness, comp_badness = PrimerSetBadnessFast(list(curr_fp.values()), list(curr_rp.values()), existing_primer)
        curr_badness += inter_badness
        print('initial loss = %.3f' % curr_badness)

        # optimization
        learning_curve = []
        SATemp = InitSATemp
        for step in range(NUMSTEPS + ZEROSTEPS):
            print('SA temperature = %.3f, SADDLE Loss = %.3f' % (SATemp, curr_badness))
            for t in range(TimePerStep):
                # pick a plex to change primer pair
                mutplex = random.choice(curr_tube) # randomly select a plex_id in current tube
                n_pair = len(all_plex_info[mutplex]['optimize'])
                if n_pair == 1:
                    continue
                mutindex = floor(rand() * n_pair) # index of new primer pair

                # copy current design
                new_index = deepcopy(curr_index)
                new_fp = deepcopy(curr_fp)
                new_rp = deepcopy(curr_rp)
                
                # calculate new badness
                new_index[mutplex] = mutindex
                new_fp[mutplex] = all_plex_info[mutplex]['optimize'][mutindex][0]['seq']
                new_rp[mutplex] = all_plex_info[mutplex]['optimize'][mutindex][1]['seq']
                new_badness = 0
                for plex_id, i in new_index.items():
                    fp = all_plex_info[plex_id]['optimize'][i][0]
                    rp = all_plex_info[plex_id]['optimize'][i][1]
                    new_badness += fp['badness'] + rp['badness']
                inter_badness, new_comp_badness = PrimerSetBadnessFast(list(new_fp.values()), list(new_rp.values()), existing_primer)
                new_badness += inter_badness
                
                # calculate probability of accepting mutation
                if new_badness < curr_badness:
                    accept = True
                else:
                    if SATemp == 0:
                        accept = False
                    elif rand() < 2**((curr_badness - new_badness)/SATemp):
                        accept = True
                    else:
                        accept = False
                        
                # implement acceptance
                if accept:
                    curr_index[mutplex] = new_index[mutplex]
                    curr_fp[mutplex] = new_fp[mutplex]
                    curr_rp[mutplex] = new_rp[mutplex]
                    curr_badness = new_badness
                    comp_badness = new_comp_badness
                    
                learning_curve.append(curr_badness)
            # reduce SATemp
            SATemp = max(0, (SATemp - InitSATemp/NUMSTEPS))
            
        # save optimization result
        for n, (plex_id, i) in enumerate(curr_index.items()):
            all_plex_info[plex_id]['fP'] = all_plex_info[plex_id]['optimize'][i][0]
            all_plex_info[plex_id]['rP'] = all_plex_info[plex_id]['optimize'][i][1]
            all_plex_info[plex_id]['fP_badness'] = comp_badness[0][n]
            all_plex_info[plex_id]['rP_badness'] = comp_badness[1][n]

        all_lc.append(learning_curve)
    
    print('Primer dimer optimization finished in %.3fs' % (time()-tik))
    print()
    return all_plex_info, all_lc


def to_df(all_plex_info, config):
    '''
    Format sequences and coordinates as a DataFrame
    Input:
        all_plex_info: output of optimize()
        config: design configurations.
    Output:
        DataFrame of all amplicons
    '''
    print('Formatting design output...')

    l_fP_adp = len(config['fP_prefix']) # length of adapter/prefix
    l_rP_adp = len(config['rP_prefix'])

    ref_name = []
    plex_name = []
    pool = []

    # sequence specific part of primers
    fp = []
    rp = []

    # coordinates (1-based)
    start = []
    insert_start = []
    insert_end = []
    end = []

    # sequences
    amp = []
    insert = []

    # full primer (with adapters/prefixes)
    fp_full = []
    rp_full = []

    # ARTIC columns
    art_genome = []
    art_start = []
    art_stop = []
    art_primer_name = []
    art_pool = []
    art_strand = []
    art_primer_seq = []

    #print('amplicons containing SNPs or iSNVs in primers:')
    for plex_id, plex_info in all_plex_info.items():
        ref_name.append(plex_info['reference'])
        plex_name.append(plex_id)
        pool.append(plex_info['tube'])
        
        curr_fp = plex_info['fP'] # fP object
        curr_rp = plex_info['rP']

        curr_fp_seq = curr_fp['seq'][l_fP_adp:] # fP sequence
        curr_rp_seq = curr_rp['seq'][l_rP_adp:]
        # if curr_fp_seq != curr_fp_seq.lower() or curr_rp_seq != curr_rp_seq.lower():
        #     print(plex_id)

        fp.append(curr_fp_seq)
        rp.append(curr_rp_seq)

        # calculate coordinates
        pc = plex_info['primers_coords'] # primer design coords ([fp_start, fp_end, rp_start, rp_end])
        curr_start = pc[1] - curr_fp['dist'] - curr_fp['primer_len'] + 1
        curr_insert_start = curr_start + curr_fp['primer_len']
        curr_insert_end = pc[2] + curr_rp['dist'] - 1
        curr_end = curr_insert_end + curr_rp['primer_len']
        start.append(curr_start)
        insert_start.append(curr_insert_start)
        insert_end.append(curr_insert_end)
        end.append(curr_end)

        curr_context_seq = plex_info['context_seq']
        amp.append(curr_context_seq[curr_start-pc[0]:curr_end-pc[0]+1])
        insert.append(curr_context_seq[curr_insert_start-pc[0]:curr_insert_end-pc[0]+1])

        fp_full.append(curr_fp['seq']) # fP with adapter/prefix
        rp_full.append(curr_rp['seq'])

        # ARTIC columns (use 0-based coordinates for BED format)
        # first fP, then rP
        art_genome.extend([plex_info['reference'], plex_info['reference']])
        art_start.extend([curr_start-1, curr_insert_end])
        art_stop.extend([curr_insert_start-1, curr_end])
        art_primer_name.extend([plex_id+'_LEFT', plex_id+'_RIGHT'])
        art_pool.extend([plex_info['tube'], plex_info['tube']])
        art_strand.extend(['+', '-'])
        art_primer_seq.extend([curr_fp_seq.upper(), curr_rp_seq.upper()])

    # format as DataFrame
    df = pd.DataFrame({
        'reference': ref_name, 
        'amplicon_id': plex_name, 
        'pool': pool, 
        'fP': fp, 
        'rP': rp, 
        'start': start, 
        'insert_start': insert_start, 
        'insert_end': insert_end, 
        'end': end, 
        'amplicon': amp, 
        'insert': insert, 
        'fP_full': fp_full, 
        'rP_full': rp_full, 
    })

    # ARTIC DataFrame
    art_df = pd.DataFrame({
        'genome': art_genome, 
        'start': art_start, 
        'stop': art_stop, 
        'primer_name': art_primer_name, 
        'pool': art_pool, 
        'strand': art_strand, 
        'primer_seq': art_primer_seq
    })
    return df, art_df


def tiling(ref_path: str, out_path: str, title: str, max_amp_len: int, min_amp_len: int, 
    w_egc: float, w_lc: float, w_ns: float, w_var: float, temperature: float, salinity: float, 
    dG_max: float, min_GC: float, max_GC: float, min_complexity: float, max_len: int, 
    check_var: bool, fP_prefix: str, rP_prefix: str, seed: int, threads: int):
    '''
    Design tiled amplicons. 
    Input:
        ref_path: Path to the Olivar reference file (.olvr), or the directory of reference files for multiple targets.
        out_path: Output directory [./].
        title: Name of this design [olivar-design].
        max_amp_len: Maximum amplicon length [420].
        min_amp_len: Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
        w_egc: Weight for extreme GC content [1.0].
        w_lc: Weight for low sequence complexity [1.0].
        w_ns: Weight for non-specificity [1.0].
        w_var: Weight for variations [1.0].
        temperature: PCR annealing temperature [60.0].
        salinity: Concentration of monovalent ions in units of molar [0.18].
        dG_max: Maximum free energy change of a primer in kcal/mol [-11.8].
        min_GC: Minimum GC content of a primer [0.2].
        max_GC: Maximum GC content of a primer [0.75].
        min_complexity: Minimum sequence complexity of a primer [0.4].
        max_len: Maximum length of a primer [36].
        check_var: Filter out primer candidates with variations within 5nt of 3' end [False]. 
            Setting check_var=True is not recommended when a lot of variations are provided, 
            since this would significantly reduce the number of primer candidates. 
        fP_prefix: Prefix of forward primer. Empty string '' by default.
        rP_prefix: Prefix of reverse primer. Empty string '' by default.
        seed: Random seed for optimizing primer design regions and primer dimer [10].
        threads: Number of threads [1].
    '''
    config = {
        'ref_path': ref_path, 
        'out_path': out_path, 
        'title': title, 
        'max_amp_len': max_amp_len, 
        'min_amp_len': min_amp_len, 
        'w_egc': w_egc, 
        'w_lc': w_lc, 
        'w_ns': w_ns, 
        'w_var': w_var, 
        'temperature': temperature, 
        'salinity': salinity, 
        'dG_max': dG_max, 
        'min_GC': min_GC, 
        'max_GC': max_GC, 
        'min_complexity': min_complexity, 
        'max_len': max_len, 
        'check_SNP': check_var, 
        'fP_prefix': fP_prefix, 
        'rP_prefix': rP_prefix, 
        'seed': seed, 
        'threads': threads
    }

    if not os.path.exists(config['out_path']):
        os.makedirs(config['out_path'])

    # store config values
    design_title = config['title']
    design_ref_path = config['ref_path']

    # validate input .olvr file(s)
    ref_path_dict = dict() # {ref_name: ref_path}
    # check input is a file or a directory
    if os.path.isfile(design_ref_path) and design_ref_path.endswith(REFEXT):
        # use the name of the .olvr file as reference name
        file = os.path.basename(design_ref_path)
        ref_path_dict[file[:-len(REFEXT)]] = design_ref_path
    elif os.path.isdir(design_ref_path):
        for file in sorted(os.listdir(design_ref_path)):
            file_path = os.path.join(design_ref_path, file)
            if os.path.isfile(file_path) and file_path.endswith(REFEXT):
                # use file name as reference name
                ref_path_dict[file[:-len(REFEXT)]] = file_path
        if ref_path_dict:
            print(f'Found {len(ref_path_dict)} {REFEXT} file(s) under "{design_ref_path}"')
            for file_path in ref_path_dict.values():
                print(os.path.basename(file_path))
        else:
            raise FileNotFoundError(f'No {REFEXT} file found in the directory "{design_ref_path}".')
    else:
        raise FileNotFoundError(f'Input is neither a {REFEXT} file nor a directory.')
    print()
    
    # design PDRs for each reference
    all_plex_info = {}
    all_ref_info = {}
    for ref_name, ref_path in ref_path_dict.items():
        config['title'] = ref_name
        config['ref_path'] = ref_path
        print(f'Designing PDRs for {ref_name} using {ref_path}...')
        temp_dict, risk_arr, gc_arr, comp_arr, hits_arr, var_arr, all_loss, seq_record = design_context_seq(config)
        all_plex_info.update(temp_dict)
        all_ref_info[ref_name] = {
            'risk_arr': risk_arr, 
            'gc_arr': gc_arr, 
            'comp_arr': comp_arr, 
            'hits_arr': hits_arr, 
            'var_arr': var_arr, 
            'all_loss': all_loss, 
            'seq_record': seq_record, 
        }
        print()
    
    # revert modified config values
    config['title'] = design_title
    config['ref_path'] = design_ref_path

    all_plex_info_primer = get_primer(all_plex_info, config)
    all_plex_info_optimize, learning_curve = optimize(all_plex_info_primer, config)
    df, art_df = to_df(all_plex_info_optimize, config)

    # format and save design output
    design_out = {
        'config': config, 
        'df': df, 
        'art_df': art_df, # ARTIC output format
        'all_plex_info': all_plex_info_optimize, 
        'learning_curve': learning_curve, 
        'all_ref_info': all_ref_info
    }
    design_path = os.path.join(config['out_path'], '%s.olvd' % config['title'])
    with open(design_path, 'wb') as f:
        pickle.dump(design_out, f, protocol=5) # protocol 5 needs python>=3.8
        print('Design file saved as %s' % design_path)

    # save human readable files
    save(design_out, config['out_path'])


def save(design_out, out_path: str):
    '''
    Load from a previous Olivar design object (in tiling()) or file (.olvd) and save output files.
    Input:
        design_out: Output of tiling(), or load the olvd file with pickle.
        out_path: Output directory.
    '''
    # load data
    if type(design_out) is str:
        print(f'Loading Olivar design from {design_out}...')
        try:
            with open(design_out,  'rb') as f:
                design_out = pickle.load(f)
                print(f'Successfully loaded Olivar design.')
        except FileNotFoundError:
            raise FileNotFoundError('Olivar design file not found.')
    
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        
    config = design_out['config']
    lc = design_out['learning_curve']
    df = design_out['df']
    all_ref_info = design_out['all_ref_info']

    # column name mapper for backward compatibility
    df.rename(columns={
        'amp_id': 'amplicon_id', 
        'amp': 'amplicon'
    }, inplace=True)

    # save configurations
    save_path = os.path.join(out_path, '%s.json' % config['title'])
    with open(save_path, 'w') as f:
        json.dump(config, f, indent=4)
    print(f'Configurations saved as {save_path}')

    # save sequences and coordinates as csv
    save_path = os.path.join(out_path, '%s.csv' % config['title'])
    df.to_csv(save_path, index=False)
    print(f'Primer pools saved as {save_path}')

    # save ARTIC output format
    try:
        art_df = design_out['art_df']
        save_path = os.path.join(out_path, '%s.scheme.bed' % config['title'])
        art_df.sort_values(by='start', axis=0, ascending=True, inplace=True)
        art_df.to_csv(save_path, sep='\t', header=False, index=False)
        print(f'Primer pools (ARTIC format) saved as {save_path}')
    except KeyError:
        print('ARTIC/PrimalScheme format is not included in older versions of Olivar, skipped.')

    #------------- SADDLE Loss -------------#
    fig = make_subplots(
        rows=1, cols=2, 
        subplot_titles=('Primer dimer optimization (pool-1)', 
            'Primer dimer optimization (pool-2)')
    )

    fig.add_trace(
        go.Scatter(
            y=lc[0], 
            line=dict(color='#1f77b4'), 
            showlegend=False
        ), 
        row=1, col=1
    )
    fig.add_trace(
        go.Scatter(
            y=lc[1], 
            line=dict(color='#1f77b4'), 
            showlegend=False
        ), 
        row=1, col=2
    )

    fig.update_xaxes(title_text='iterations', row=1, col=1)
    fig.update_xaxes(title_text='iterations', row=1, col=2)

    fig.update_yaxes(title_text='SADDLE Loss', row=1, col=1)
    fig.update_yaxes(title_text='SADDLE Loss', row=1, col=2)

    # save html figure
    save_path = os.path.join(out_path, f'{config['title']}_SADDLE_Loss.html')
    with open(save_path, 'w') as f:
        f.write(plotly.io.to_html(fig))
    print(f'SADDLE optimization plot saved as {save_path}')
    #------------- SADDLE Loss -------------#

    for ref_name, ref_info in all_ref_info.items():
        print(f'Saving output files and figures for {ref_name}...')

        risk_arr = ref_info['risk_arr']
        gc_arr = ref_info['gc_arr']
        comp_arr = ref_info['comp_arr']
        hits_arr = ref_info['hits_arr']
        var_arr = ref_info['var_arr']
        seq_record = ref_info['seq_record']
        all_loss = ref_info['all_loss']

        # save reference sequence
        save_path = os.path.join(out_path, '%s.fasta' % ref_name)
        with open(save_path, 'w') as f:
            SeqIO.write([seq_record], f, 'fasta')
        print(f'Reference sequence saved as {save_path}')

        # save risk array
        risk = pd.DataFrame({
            'position': range(1, len(risk_arr)+1), 
            'base': list(seq_record.seq), 
            'extreme GC': gc_arr, 
            'low complexity': comp_arr, 
            'non-specificity': hits_arr, 
            'variations': var_arr, 
            'risk': risk_arr
        })
        save_path = os.path.join(out_path, '%s_risk.csv' % ref_name)
        risk.to_csv(save_path, index=False)
        print(f'Risk scores saved as {save_path}')

        #------------- PDR Loss -------------#
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                y=all_loss, 
                line=dict(color='#1f77b4'), 
                showlegend=False
            )
        )
        fig.update_layout(
            title=f'Optimization of primer design regions (PDRs) for {ref_name}', 
            xaxis_title='iterations (sorted)', 
        )
        fig.update_yaxes(title_text='Loss', type='log')
        # save html figure
        save_path = os.path.join(out_path, f'{ref_name}_PDR_Loss.html')
        with open(save_path, 'w') as f:
            f.write(plotly.io.to_html(fig))
        print(f'PDR optimization plot saved as {save_path}')
        #------------- PDR Loss -------------#

        #------------- risk array and primers -------------#
        # create figure
        fig = go.Figure()

        # set figure scale
        base_offset = (0.5*config['w_egc'] + 0.5*config['w_lc'] + config['w_ns'])/6
        base_offset *= 3

        # plot risk array
        r = np.arange(len(risk_arr))
        fig.add_trace(
            go.Scatter(
                x=r, y=gc_arr,
                hoverinfo='skip',
                mode='lines',
                line=dict(width=0, color='#1f77b4'),
                name='extreme GC',
                stackgroup='one' # define stack group
            )
        )
        fig.add_trace(
            go.Scatter(
                x=r, y=comp_arr,
                hoverinfo='skip',
                mode='lines',
                line=dict(width=0, color='#ff7f0e'),
                name='low complexity',
                stackgroup='one'
            )
        )
        fig.add_trace(
            go.Scatter(
                x=r, y=hits_arr,
                hoverinfo='skip',
                mode='lines',
                line=dict(width=0, color='#2ca02c'),
                name='non-specificity',
                stackgroup='one'
            )
        )
        fig.add_trace(
            go.Scatter(
                x=r, y=var_arr,
                hoverinfo='skip',
                mode='lines',
                line=dict(width=0, color='#d62728'),
                name='variations',
                stackgroup='one'
            )
        )

        # plot primers
        fp_rp_diff = 0.1 # distance between fP and rP
        head_dx = 7 # arrow head length
        head_dy = base_offset*7/90 # arrow head height
        primer_x = [] # x coords for primer plot
        primer_y = [] # y coords for primer plot
        hover_text = []

        # pool 1
        fp_offset = 2 * base_offset
        rp_offset = (2-fp_rp_diff) * base_offset
        for i, row in df[(df['reference']==ref_name) & (df['pool']==1)].iterrows():
            fp_start = row['start']-1
            fp_stop = row['insert_start']-1
            rp_start = row['insert_end']
            rp_stop = row['end']
            # [fp_start, fp_stop] and [fp_offset, fp_offset] is the body of fP
            # [fp_stop-head_dx, fp_stop] and [fp_offset+head_dy, fp_offset] is the head of fP
            primer_x.extend([fp_start, fp_stop, None, fp_stop-head_dx, fp_stop, None, 
                rp_start, rp_stop, None, rp_start, rp_start+head_dx, None])
            primer_y.extend([fp_offset, fp_offset, None, fp_offset+head_dy, fp_offset, None, 
                rp_offset, rp_offset, None, rp_offset, rp_offset-head_dy, None])
            hover_text.extend(['pool-1 %s fP: %d - %d' % (row['amplicon_id'], row['start'], fp_stop)]*6 + \
                ['pool-1 %s rP: %d - %d' % (row['amplicon_id'], rp_start+1, rp_stop)]*6)
        
        # pool 2
        fp_offset = base_offset
        rp_offset = (1-fp_rp_diff) * base_offset
        for i, row in df[(df['reference']==ref_name) & (df['pool']==2)].iterrows():
            fp_start = row['start']-1
            fp_stop = row['insert_start']-1
            rp_start = row['insert_end']
            rp_stop = row['end']
            # [fp_start, fp_stop] and [fp_offset, fp_offset] is the body of fP
            # [fp_stop-head_dx, fp_stop] and [fp_offset+head_dy, fp_offset] is the head of fP
            primer_x.extend([fp_start, fp_stop, None, fp_stop-head_dx, fp_stop, None, 
                rp_start, rp_stop, None, rp_start, rp_start+head_dx, None])
            primer_y.extend([fp_offset, fp_offset, None, fp_offset+head_dy, fp_offset, None, 
                rp_offset, rp_offset, None, rp_offset, rp_offset-head_dy, None])
            hover_text.extend(['pool-2 %s fP: %d - %d' % (row['amplicon_id'], row['start'], fp_stop)]*6 + \
                ['pool-2 %s rP: %d - %d' % (row['amplicon_id'], rp_start+1, rp_stop)]*6)
        
        # plot all primers
        fig.add_trace(
            go.Scatter(
                x=primer_x, y=primer_y, 
                mode='lines', # connect points with lines
                connectgaps=False, # gaps (np.nan or None) in the provided data arrays are not connected
                line=dict(color='rgb(0,0,0)'), # black
                showlegend=False, 
                hovertemplate='%{text}<extra></extra>', # <extra></extra> hides the "trace 4"
                text=hover_text
            )
        )

        # title and axis
        fig.update_layout(
            hoverlabel=dict(
                bgcolor='white'
            ), 
            # title
            title_text="risk components are stacked together", 
            title_x=0.98, 
            titlefont=dict(
                size=18,
            ), 

            # axis
            xaxis=dict(
                tickformat='%d', 
                tickfont=dict(
                    size=14, 
                ), 
                showgrid=False, 
                rangeslider=dict(
                    visible=True
                )
            ), 
            yaxis=dict(
                range=[0, base_offset*6], 
                tickfont=dict(
                    size=14, 
                ), 
                showgrid=False
            ), 

            # axis title
            xaxis_title='position', 
            yaxis_title='risk', 

            # global font
            font=dict(
                size=16,
            )
        )

        # save html figure
        save_path = os.path.join(out_path, '%s.html' % ref_name)
        with open(save_path, 'w') as f:
            f.write(plotly.io.to_html(fig))
        print(f'Risk and primer viewer saved as {save_path}')
        #------------- risk array and primers -------------#
        print()


def validate(primer_pool: str, pool: int, BLAST_db: str, out_path: str, 
    title: str, max_amp_len: int, temperature: float, threads: int):
    '''
    Run analysis on a primer pool for multiplex PCR
    Input:
        primer_pool: Path to the csv file of a primer pool. 
            Required columns: "amplicon_id" (amplicon name), "fP" (sequence of forward primer), "rP" (sequence of reverse primer), "pool" (pool number, e.g., 1).
        pool: Primer pool number [1]. 
        BLAST_db: Optional, path to the BLAST database. 
            Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        out_path: Output directory [./]. 
        title: Name of validation [olivar-val].
        max_amp_len: Maximum length of predicted non-specific amplicon [1500]. 
        temperature: PCR annealing temperature [60.0]. 
        threads: Number of threads [1]. 
    '''
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    seq_list = []
    seq_names = []
    df = pd.read_csv(primer_pool, sep=',', index_col=False)

    # column name mapper for backward compatibility
    df.rename(columns={
        'amp_id': 'amplicon_id', 
        'amp': 'amplicon'
    }, inplace=True)
    
    for i, row in df[df['pool']==pool].iterrows():
        seq_list.extend([row['fP'], row['rP']])
        seq_names.extend(['%s_fP' % row['amplicon_id'], '%s_rP' % row['amplicon_id']])
    print(f'Successfully loaded {primer_pool}, pool-{pool}.')

    # non-specific simulation
    if BLAST_db:
        df_ns, df_count, all_hits = ns_simulation(BLAST_db, seq_list, seq_names, max_amp_len, threads)
        save_path = os.path.join(out_path, f'{title}_pool-{pool}_ns-amp.csv')
        df_ns.to_csv(save_path, index=False)
        print('Non-specific amplicons saved as %s' % save_path)
        save_path = os.path.join(out_path, f'{title}_pool-{pool}_ns-pair.csv')
        df_count.to_csv(save_path, index=False)
        print('Non-specific primer pairs saved as %s' % save_path)
    else:
        all_hits = np.zeros(len(seq_list)) - 1
        print('No BLAST database provided, skipped.')

    # calculate Badness
    #all_rep = [seq2arr(p) for p in seq_list]
    _, all_bad = PrimerSetBadnessFast(seq_list)
    all_bad = all_bad[0]

    # calculate dG
    gen = design.primer_generator(temperature=temperature, salinity=0.18)
    all_dG = [gen.dG_init+gen.StacksDG(p) for p in seq_list]

    df_val = pd.DataFrame({
        'name': seq_names, 
        'seq': seq_list, 
        'length': [len(p) for p in seq_list], 
        '%GC': [basic.get_GC(p) for p in seq_list], 
        'complexity': [basic.get_complexity(p) for p in seq_list], 
        'dG (%dC)' % temperature: all_dG, 
        'self_align': [design.WAlignScore(p) for p in seq_list], 
        'dimer_score': all_bad, 
        'BLAST_hits': all_hits
    })
    save_path = os.path.join(out_path, f'{title}_pool-{pool}.csv')
    df_val.to_csv(save_path, index=False)
    print('Validation file saved as %s' % save_path)
