import logging
logger = logging.getLogger('main')

import os
import pickle
from Bio import SeqIO
import pandas as pd
import numpy as np
import multiprocessing
from time import time

# internal modules
import basic
from ncbi_tools import BLAST_batch_short

REFEXT = '.olvr' # extension for Olivar reference file

def run_build(fasta_path: str, var_path: str, BLAST_db: str, out_path: str, title: str, threads: int):
    '''
    Build the Olivar reference file for tiled amplicon design
    Input:
        fasta_path: Path to the FASTA reference sequence.
        var_path: Optional, path to the csv file of SNP coordinates and frequencies. 
            Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.
        BLAST_db: Optional, path to the BLAST database. 
            Note that this path should end with the name of the BLAST database (e.g., "example_input/Human/GRCh38_primary").
        out_path: Output directory [./]. 
        title: Name of the Olivar reference file [FASTA record ID]. 
        threads: Number of threads [1]. 
    '''
    n_cpu = threads

    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    # load the first record in fasta_path
    logger.info(f'Loading the first record of "{fasta_path}"')
    for record in SeqIO.parse(fasta_path, 'fasta'):
        seq_raw = str(record.seq).lower() # SNPs in upper case
        break # read first record
    seq_id = str(record.id) # ID of the fasta record (everything before the space character in the FASTA header)
    if title is None:
        title = seq_id
    save_path = os.path.join(out_path, f'{title}{REFEXT}')
    if os.path.exists(save_path):
        raise FileExistsError(f'Olivar reference file already exists "{save_path}"')

    word_size = 28
    offset = 14 # word_size should be divisible by offset
    logger.info('Building Olivar reference with word_size=%d, offset=%d' % (word_size, offset))

    n_cycle = word_size/offset
    if n_cycle != int(n_cycle):
        raise ValueError('word_size should be divisible by offset.')
    else:
        n_cycle = int(n_cycle)
    
    # check ambiguous bases
    ambiguous_bases = set(seq_raw) - {'a', 't', 'c', 'g'}
    if ambiguous_bases:
        raise ValueError(f'Ambiguous bases are not supported. {ambiguous_bases} found in reference fasta file.')

    # load SNP/iSNV coordinates
    # all coordinates are 1-based and inclusive
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
                logger.warning("'FREQ' is not provided")
                df[cname] = np.nan
                continue
            if received_dtype != expected_dtype:
                try:
                    # Attempt to cast the column to the expected dtype
                    df[cname] = df[cname].astype(expected_dtype)
                    logger.info(f"Column '{cname}' was converted to {expected_dtype}.")
                except ValueError:
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
        logger.info(f'Variation coordinates and frequencies loaded from "{var_path}"')
    else:
        logger.info('No sequence variation coordinates provided, skipped.')

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
        logger.info('Calculating GC content and sequence complexity...')
        tik = time()
        all_gc = np.array(pool.map(basic.get_GC, all_word))
        all_complexity = np.array(pool.map(basic.get_complexity, all_word))
        logger.info(f'Finished in {time()-tik:.3f}s')
    if BLAST_db:
        logger.info(f'Calculating non-specificity with BLAST database "{BLAST_db}"')
        all_hits, _ = BLAST_batch_short(all_word, db=BLAST_db, n_cpu=n_cpu, mode='rough') # tabular output
    else:
        logger.info('No BLAST database provided, skipped.')
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
        'seq_id': seq_id, 
        'seq_record': record, 
        'start': start, 
        'stop': stop, 
        'gc_arr': gc_arr, 
        'comp_arr': comp_arr, 
        'var_arr': var_arr[start-1:stop], 
        'hits_arr': hits_arr, 
        'all_hits': all_hits
    }
    with open(save_path, 'wb') as f:
        pickle.dump(olv_ref, f, protocol=5) # protocol 5 needs python>=3.8
        logger.info('Reference file saved as %s' % save_path)