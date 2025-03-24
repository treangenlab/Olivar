#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Main workflow of Olivar tiling.
Architecture:
main.py
    build()
        run_build()
    tiling()
        design_context_seq()
            generate_context()
                find_min_loc()
        get_primer()
        optimize()
        to_df()
    save()
    specificity()
    sensitivity()
'''


__author__ = 'Michael X. Wang'


import os
CURRENT_DIR = os.path.dirname(__file__)

import sys
import io
import json
import pickle # needs python>=3.8

import logging

# logging settings
logFormatter = logging.Formatter(
    fmt='%(asctime)s | %(levelname)-8s | %(message)s', 
    style='%', 
    datefmt='%Y-%m-%d %H:%M:%S'
)
# add handler for printing to stdout
consoleHandler = logging.StreamHandler(sys.stdout)
consoleHandler.setFormatter(logFormatter)
# create logger
logger = logging.getLogger('main')
logger.setLevel(logging.DEBUG) # show all messages
logger.addHandler(consoleHandler)

import pandas as pd
import numpy as np

# internal modules
import basic
import design # SADDLE
from ncbi_tools import ns_simulation
from design import PrimerSetBadnessFast
from msa_tools import run_validate, run_preprocess, run_cmd

# helper modules
from build_helper import run_build
from tiling_helper import design_context_seq, get_primer, optimize, to_df
from save_helper import save

REFEXT = '.olvr' # extension for Olivar reference file
DESIGNEXT = '.olvd' # extension for Olivar design file

def build(fasta_path: str, msa_path: str, var_path: str, BLAST_db: str, out_path: str, title: str, threads: int, align: bool, min_var: float):
    '''
    Build the Olivar reference file for tiled amplicon design, handling gaps in the MSA.
    Input:
        fasta_path: Optional, Path to the FASTA reference sequence.
        msa_path: Optional, Path to the MSA file (Multiple Sequence Alignment in FASTA format).
        var_path: Optional, path to the CSV file of SNP coordinates and frequencies. 
            Required columns: "START", "STOP", "FREQ". "FREQ" is considered as 1.0 if empty. Coordinates are 1-based.
        BLAST_db: Optional, path to the BLAST database. 
        out_path: Output directory [./]. 
        title: Name of the Olivar reference file [MSA record ID]. 
        threads: Number of threads [1]. 
        align: Conrol whether do alignment for MSA file or not [False]. 
        min_var: Minimum threshold of frequencies of SNP [0.01].
    '''
    if not msa_path and not fasta_path:
        raise ValueError("Either 'msa_path' or 'fasta_path' must be provided.")
    if msa_path:
        if fasta_path:
            logger.warning("Both 'msa_path' and 'fasta_path' provided. Ignoring 'fasta_path' and processing with MSA + BLAST.")
        if not os.path.exists(msa_path):
            raise FileNotFoundError(f"MSA file '{msa_path}' not found.")
        fasta_path = None
        if align:
            logger.info("Running alignment with MAFFT...")
            msa_filename = os.path.splitext(os.path.basename(msa_path))[0]
            aligned_msa_path = os.path.join(out_path, os.path.basename(msa_path).replace('.fasta', '_aligned.fasta'))
            msa_str = run_cmd('mafft', '--auto', '--thread', str(threads), msa_path)

            if not os.path.exists(out_path):
                os.makedirs(out_path)
            with open(aligned_msa_path, 'w') as f:
                f.write(msa_str)
            msa_path = aligned_msa_path
        #fasta_path, var_path = run_preprocess(io.StringIO(msa_str), msa_filename, out_path, threads)
        fasta_path, var_path = run_preprocess(msa_path, msa_filename, out_path, threads, min_var)
    else:
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file '{fasta_path}' not found.")
    run_build(fasta_path, var_path, BLAST_db, out_path, title, threads)


def tiling(ref_path: str, out_path: str, title: str, max_amp_len: int, min_amp_len: int, 
    w_egc: float, w_lc: float, w_ns: float, w_var: float, temperature: float, salinity: float, 
    dG_max: float, min_GC: float, max_GC: float, min_complexity: float, max_len: int, 
    check_var: bool, fP_prefix: str, rP_prefix: str, seed: int, threads: int, num_iterations: int):
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
        num_iterations: Number of iterations [1].
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
        'threads': threads,
        'num_iterations': num_iterations
    }

    if not os.path.exists(config['out_path']):
        os.makedirs(config['out_path'])

    # store config values
    design_ref_path = config['ref_path']

    # validate input .olvr file(s)
    ref_path_dict = dict() # {ref_name: ref_path}


    # log to file, must happen after working_dir is created
    # if working_dir exist, it should be a directory as well
    fileHandler = logging.FileHandler(os.path.join(out_path, f'{config['title']}.log'), mode='a')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

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
            logger.info(f'Found {len(ref_path_dict)} {REFEXT} file(s) under "{design_ref_path}"')
            for file_path in ref_path_dict.values():
                logger.info(os.path.basename(file_path))
        else:
            raise FileNotFoundError(f'No {REFEXT} file found in the directory "{design_ref_path}".')
    else:
        raise FileNotFoundError(f'Input is neither a {REFEXT} file nor a directory.')
    
    # design PDRs for each reference
    all_plex_info = {}
    all_ref_info = {}
    for ref_name, ref_path in ref_path_dict.items():
        config['ref_path'] = ref_path
        logger.info(f'Designing PDRs for {ref_name} using {ref_path}...')
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
    
    # revert modified config values
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
        logger.info('Design file saved as %s' % design_path)

    # save human readable files
    save(design_out, config['out_path'])

def specificity(primer_pool: str, pool: int, BLAST_db: str, out_path: str, 
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
    logger.info(f'Successfully loaded {primer_pool}, pool-{pool}.')

    # non-specific simulation
    if BLAST_db:
        df_ns, df_count, all_hits = ns_simulation(BLAST_db, seq_list, seq_names, max_amp_len, threads)
        save_path = os.path.join(out_path, f'{title}_pool-{pool}_ns-amp.csv')
        df_ns.to_csv(save_path, index=False)
        logger.info('Non-specific amplicons saved as %s' % save_path)
        save_path = os.path.join(out_path, f'{title}_pool-{pool}_ns-pair.csv')
        df_count.to_csv(save_path, index=False)
        logger.info('Non-specific primer pairs saved as %s' % save_path)
    else:
        all_hits = np.zeros(len(seq_list)) - 1
        logger.info('No BLAST database provided, skipped.')

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
    logger.info('Validation file saved as %s' % save_path)

def sensitivity(primer_pool: str, msa_path: str, pool: int, out_path: str, 
    title: str, temperature: float=60.0, sodium: float=0.18, threads: int=1, align: bool=False):
    '''
    Visualize an MSA, along with its primers/probes (if provided), and validate their alignment and sensitivity.
    Input:
        primer_pool: Path to the csv file of a primer pool. 
            Required columns: "amplicon_id" (amplicon name), "fP" (sequence of forward primer), "rP" (sequence of reverse primer), "pool" (pool number, e.g., 1).
        msa_path: Path to the MSA file.
        pool: Primer pool number [1]. 
        out_path: Output directory [./]. 
        title: Name of validation [olivar-val].
        temperature: PCR annealing temperature [60.0]. 
        sodium: The sum of the concentrations of monovalent ions (Na+, K+, NH4+), in molar [0.18].
        threads: Number of threads [1]. 
        align: Conrol whether do alignment for MSA file or not. [False]. 
    '''
    if not msa_path:
        raise ValueError("'msa_path' must be provided for sensitivity calculation.")
    
    if align:
        logger.info("Running alignment with MAFFT...")
        aligned_msa_path = os.path.join(out_path, os.path.basename(msa_path).replace('.fasta', '_aligned.fasta'))
        msa_str = run_cmd('mafft', '--auto', '--thread', str(threads), msa_path)
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        with open(aligned_msa_path, 'w') as f:
            f.write(msa_str)
        msa_path = aligned_msa_path
    
    run_validate(msa_path = msa_path, primer_pool = primer_pool, pool = pool, temperature = temperature, sodium = sodium, out_path = out_path, title = title, n_cpu = threads)
