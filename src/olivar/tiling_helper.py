import logging
logger = logging.getLogger('main')

import numpy as np
import pandas as pd
import pickle
import random
import multiprocessing

from time import time
from numpy.random import rand, default_rng
from copy import deepcopy
from tqdm import tqdm
from math import floor

# imternal modules
import basic
import design

# algorithmic parameters
PRIMER_DESIGN_LEN = 40 # length of primer design region [40]
UPPER_GC = 0.75 # extreme GC upper bond [0.75]
LOWER_GC = 0.25 # extreme GC lower bond [0.75]
LOW_COMPLEXITY = 0.4 # low complexity lower bond [0.4]
CHOICE_RATE = 0.3 # bottom CHOICE_RATE lowest risk primer design regions are choosen from candidate region
RISK_TH = 0.1 # top RISK_TH highest risk primer design regions are considered as loss

N_POOLS = 2 # number of primer pool (not to be changed)

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
    iterMul = config['iterMul']
    
    # set random number generator
    rng_parent = default_rng(seed)
    
    with open(ref_path,  'rb') as f:
        olv_ref = pickle.load(f)

    seq_raw = olv_ref['seq'] # SNPs are capitalized in build()
    seq_id = olv_ref['seq_id'] # ID of the fasta record (everything before the space character in the FASTA header)
    start = olv_ref['start'] # start of design region, 1-based, closed
    stop = olv_ref['stop'] # stop of design region, 1-based, closed
    gc_arr = olv_ref['gc_arr']
    comp_arr = olv_ref['comp_arr']
    var_arr = olv_ref['var_arr']
    hits_arr = olv_ref['hits_arr']
    logger.info('Successfully loaded reference file %s' % ref_path)
    
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

    N = iterMul * 500*len(risk_arr)//max_amp_len # number of primer sets to generate
    #N = 100
    rand_int = rng_parent.integers(2**32, size=N) # random seeds for each iteration

    # single thread
    # design = [generate_context((risk_arr, start, stop, max_amp_len, min_amp_len, rand_int[i])) for i in tqdm(range(N))]

    # multi threads
    logger.info('reference sequence length: %d' % len(seq_raw))
    logger.info('design region: %d:%d' % (start, stop))
    logger.info('region length: %d' % (stop-start+1))
    logger.info('number of PDR sets to be tested: %d' % N)
    logger.info('Designing PDRs...')
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
    logger.info('Finished in %.3fs' % (time()-tik))

    # Loss of all iterations
    all_loss = [d[2] for d in design]
    all_loss.sort(reverse=True)

    # find the best arangement of primer design regions
    best_design = sorted(design, key=lambda x:x[2])
    all_context_seq, all_risk, loss = best_design[0]
    logger.info('Loss of the best design: %.3f' % loss)
    
    # prepare output
    all_plex_info = {}
    for i, (context_seq, risk) in enumerate(zip(all_context_seq, all_risk)):
        cs = seq_raw[context_seq[0]-1:context_seq[3]] # actual context sequence
        plex_info = {
            'reference': seq_id, 
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
        all_plex_info['%s_%d' % (seq_id, i+1)] = plex_info
    logger.info('total amplicons: %d' % (i+1))
    cover_start = all_context_seq[0][1]+1
    cover_stop = all_context_seq[-1][2]-1
    logger.info('covered region: %d:%d' % (cover_start, cover_stop))
    logger.info('coverage of reference sequence: %.3f%%' % (100*(cover_stop-cover_start+1)/len(seq_raw)))
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
    logger.info('Generating primer candidates...')
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
            logger.info('Fail to generate fP, update setting: %s' % plex_id)
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
            logger.info('Fail to generate rP, update setting: %s' % plex_id)
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
        
    logger.info('Finished in %.3fs' % (time()-tik))
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
    logger.info('Optimizing primer dimers...')
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
    logger.info('total primer pairs %d' % sum(n_pair))
    logger.info('average pairs per plex %.2f' % (sum(n_pair)/len(n_pair)))
    
    # optimize each tube
    all_lc = []
    tik = time()
    for i_tube in range(1, N_POOLS+1):
        logger.info('\npool %d, simulated annealing...' % i_tube)
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
        inter_badness, comp_badness = design.PrimerSetBadnessFast(list(curr_fp.values()), list(curr_rp.values()), existing_primer)
        curr_badness += inter_badness
        logger.info('initial loss = %.3f' % curr_badness)

        # optimization
        learning_curve = []
        SATemp = InitSATemp
        for step in range(NUMSTEPS + ZEROSTEPS):
            logger.info('SA temperature = %.3f, SADDLE Loss = %.3f' % (SATemp, curr_badness))
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
                inter_badness, new_comp_badness = design.PrimerSetBadnessFast(list(new_fp.values()), list(new_rp.values()), existing_primer)
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
    
    logger.info('Primer dimer optimization finished in %.3fs' % (time()-tik))
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
    logger.info('Formatting design output...')

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
        art_primer_name.extend([plex_id+'_LEFT_1', plex_id+'_RIGHT_1'])
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
