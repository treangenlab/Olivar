#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Michael X. Wang'
__version__ = '0.0.1'

import os
CURRENT_DIR = os.path.dirname(__file__)

import sys
import argparse
import logging
import datetime
import subprocess
import multiprocessing
from time import time

logger = logging.getLogger('main')

import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from plotly import graph_objects as go
import pandas as pd

# seed random generator
rng = np.random.default_rng(seed=0)

# Define DNA bases and their mapping to numbers
DNA = ('A', 'T', 'C', 'G', '-')
DNANUM = {'A':0, 'T':1, 'C':2, 'G':3}

# Define parameters for dG calculation
paraH = np.array([
    [-7.6, -7.2, -8.4, -7.8], 
    [-7.2, -7.6, -8.2, -8.5], 
    [-8.5, -7.8, -8.0, -10.6], 
    [-8.2, -8.4, -9.8, -8.0], 
])
paraS = np.array([
    [-21.3, -20.4, -22.4, -21.0], 
    [-21.3, -21.3, -22.2, -22.7], 
    [-22.7, -21.0, -19.9, -27.2], 
    [-22.2, -22.4, -24.4, -19.9], 
])

# Define ambiguous bases
AMBFREQ = {        # A,   T,   C,   G,   -
    'A': np.array([1.0, 0.0, 0.0, 0.0, 0.0]), 
    'T': np.array([0.0, 1.0, 0.0, 0.0, 0.0]), 
    'C': np.array([0.0, 0.0, 1.0, 0.0, 0.0]), 
    'G': np.array([0.0, 0.0, 0.0, 1.0, 0.0]), 
    '-': np.array([0.0, 0.0, 0.0, 0.0, 1.0]), 
    'U': np.array([0.0, 1.0, 0.0, 0.0, 0.0]), 
    'N': np.array([1/4, 1/4, 1/4, 1/4, 0.0]), 
    'R': np.array([1/2, 0.0, 0.0, 1/2, 0.0]), 
    'Y': np.array([0.0, 1/2, 1/2, 0.0, 0.0]), 
    'K': np.array([0.0, 1/2, 0.0, 1/2, 0.0]), 
    'M': np.array([1/2, 0.0, 1/2, 0.0, 0.0]), 
    'S': np.array([0.0, 0.0, 1/2, 1/2, 0.0]), 
    'W': np.array([1/2, 1/2, 0.0, 0.0, 0.0]), 
    'B': np.array([0.0, 1/3, 1/3, 1/3, 0.0]), 
    'D': np.array([1/3, 1/3, 0.0, 1/3, 0.0]), 
    'H': np.array([1/3, 1/3, 1/3, 0.0, 0.0]), 
    'V': np.array([1/3, 0.0, 1/3, 1/3, 0.0]), 
}
# ambiguous base mapping
AMBMAPPING = (
    ('A', 'A'), ('A', 'N'), ('A', 'R'), ('A', 'M'), ('A', 'W'), ('A', 'D'), ('A', 'H'), ('A', 'V'), 
    ('T', 'T'), ('T', 'N'), ('T', 'Y'), ('T', 'K'), ('T', 'W'), ('T', 'B'), ('T', 'D'), ('T', 'H'), 
    ('C', 'C'), ('C', 'N'), ('C', 'Y'), ('C', 'M'), ('C', 'S'), ('C', 'B'), ('C', 'H'), ('C', 'V'), 
    ('G', 'G'), ('G', 'N'), ('G', 'R'), ('G', 'K'), ('G', 'S'), ('G', 'B'), ('G', 'D'), ('G', 'V'), 
    ('U', 'U'), ('U', 'T'), ('U', 'N'), ('U', 'Y'), ('U', 'K'), ('U', 'W'), ('U', 'B'), ('U', 'D'), ('U', 'H'), 
)
AMBREP = {
    'N': ('A', 'T', 'C', 'G'), 
    'R': ('A', 'G'), 
    'Y': ('T', 'C'), 
    'K': ('T', 'G'), 
    'M': ('A', 'C'), 
    'S': ('C', 'G'), 
    'W': ('A', 'T'), 
    'B': ('T', 'C', 'G'), 
    'D': ('A', 'T', 'G'), 
    'H': ('A', 'T', 'C'), 
    'V': ('A', 'C', 'G')
}
ALPHABET = set(AMBFREQ)

# create substitution matrix for local pairwise alignment
from Bio.Align.substitution_matrices import Array
SUBMTX = Array(tuple(AMBFREQ), dims=2) - 1 # set mismatch_score = -1
for b1, b2 in AMBMAPPING:
    SUBMTX[b1, b2] = 1  # set match_score = -1
    SUBMTX[b2, b1] = 1
# create pairwise aligner
from Bio.Align import PairwiseAligner
aligner_local = PairwiseAligner(
    # match_score and mismatch_score are defined in substitution_matrix
    mode='local', 
    open_gap_score=-2, 
    extend_gap_score=-1, 
    substitution_matrix=SUBMTX
)
aligner_global = PairwiseAligner(
    mode='global', 
    open_gap_score=-2, 
    extend_gap_score=-1, 
)


def print_time_delta(seconds):
    logger.info(f'Finished in {datetime.timedelta(seconds=seconds)}')


def log_and_raise(exception: Exception, msg: str):
    logger.error(msg)
    raise exception(msg)


def run_cmd(*args) -> str:
    for a in args:
        if type(a) is not str:
            log_and_raise(TypeError, 'Command line arguments should be converted to str')
    cmd_out = subprocess.run(args, capture_output=True, text=True)
    # capture error message
    if cmd_out.returncode != 0:
        log_and_raise(Exception, cmd_out.stderr)
    return cmd_out.stdout


def mp_wrapper(func: callable, all_args: iter, n_cpu: int=1, text: str|None=None) -> list:
    '''Wrapper for multiprocessing.Pool()
    Args:
        func (callable): function for multiprocessing
        all_args (iter): k-mer length for minimizer sketch
        n_cpu (int): number of processes to run in parallele [1]
        text (str): message to be printed when multiprocessing starts [None]

    Returns:
        func_out (list): func outputs in a list in the same order as all_args
    '''
    if text:
        logger.info(f'{text} (threads={n_cpu})')
    if n_cpu == 1:
        func_out = [func(*args) for args in tqdm(all_args, ascii=' >')]
    elif n_cpu > 1:
        with multiprocessing.Pool(processes=n_cpu) as pool:
            tik = time()
            func_out = pool.starmap(func, all_args)
    else:
        log_and_raise(ValueError, 'n_cpu should be an integer')
    if text:
        print_time_delta(time()-tik)
    return func_out


def n_chunk(ls: list, n: int):
    '''Yield n rounghly same size chunks of from ls.
    '''
    l = len(ls)
    remainder = l%n
    size = l//n
    stop = 0
    for i in range(n):
        start = stop
        if i < remainder:
            stop = start + size + 1
            yield ls[start: stop]
        else:
            stop = start + size
            yield ls[start: stop]


def count_bases(seq) -> np.ndarray:
    '''Calculate frequencies of 'A', 'T', 'C', 'G', '-'
    '''
    try:
        c = sum([AMBFREQ[b] for b in seq])
    except KeyError:
        # from None: to suppress "During handling of the above exception, another exception occurred"
        # so that only the ValueError is shown to the user
        # otherwise, the KeyError is also shown, which is redundant
        raise ValueError(
            f'{set(seq)-ALPHABET} does not belong to the nucleotide alphabet, but found in input sequence.'
        ) from None
    #return {'A': c[0], 'T': c[1], 'C': c[2], 'G': c[3], '-': c[4]}, c
    return c


def msa_to_matrix(msa: iter) -> np.ndarray:
   return np.array([list(str(s).upper()) for s in msa])


def seq2arr(seq: str) -> list:
    '''DNA sequence (string) to an array (list of numbers)
    Args: 
        seq (str): DNA sequence
    
    Returns:
        seq_arr (list): list of numbers ('A':0, 'T':1, 'C':2, 'G':3)
    '''
    try:
        seq_arr = [DNANUM[b] for b in seq]
    except KeyError:
        raise ValueError(
            f'{set(seq)-set(DNANUM)} does not belong to the nucleotide alphabet, but found in input sequence (ambiguous bases are not allowed here).'
        ) from None
    return seq_arr


class Primer(Seq):
    '''The primer object, inherited from Bio.Seq.Seq
    Args: 
        seq (str): primer sequence (ambiguous bases are supported)
        temperature (float): annealing temperature in Degree Celsius [60]
        sodium (float): the sum of the concentrations of monovalent ions (Na+, K+, NH4+), in molar [0.18]
    '''
    def __init__(self, seq: str, temperature: float=60.0, sodium: float=0.18) -> None:
        # sanity check
        seq = seq.upper()
        if not set(seq).issubset(ALPHABET):
            raise ValueError(f'{set(seq)-ALPHABET} does not belong to the nucleotide alphabet, but found in input sequence.')
        elif len(seq) < 2:
            raise ValueError(f'Primer sequence is shorter than 2nt.')
        
        # initialize with Bio.Seq.Seq
        super().__init__(seq)
        
        # calculate params for dG calculation
        self._paraS = paraS + 0.368 * np.log(sodium)
        self._paraG = paraH - (temperature + 273.15) * self._paraS / 1000
        self._dG_init = 0.2 - (temperature + 273.15) * (-5.7)/1000
        
        # get all non-ambigous sequences and calculate dG and GC content
        self.expanded = self._expand()
        self.dG_expanded = {s: self._get_dG(s) for s in self.expanded}
        self.dG = float(np.mean(list(self.dG_expanded.values())))
        self.GC_expanded = {s: gc_fraction(s) for s in self.expanded}
        self.GC = float(np.mean(list(self.GC_expanded.values())))

        self.temperature = temperature
        self.sodium = sodium

    def _expand(self) -> tuple:
        '''return a tuple of all non-ambigous sequences
        '''
        n_perm = 1 # number of non-ambiguous sequences
        amb_pos = list() # position of ambiguous bases
        for pos, b in enumerate(self.__str__()):
            try:
                amb = AMBREP[b]
                amb_pos.append((amb, pos))
                n_perm *= len(amb)
            except KeyError:
                continue
        # stack all sequence permutations into a matrix
        mtx = np.array([list(self.__str__())]*n_perm)
        # change columns of ambiguous bases
        n_repeat = n_perm
        n_tile = 1
        for amb, pos in amb_pos:
            n_repeat /= len(amb)
            mtx[:, pos] = np.tile(np.repeat(amb, n_repeat), n_tile)
            n_tile *= len(amb)
        return tuple([Seq(''.join(row)) for row in mtx])
    
    def _get_dG(self, seq: str) -> float:
        '''return dG in kcal/mol
        '''
        rep = seq2arr(seq)
        dG = self._dG_init
        for i, _ in enumerate(rep[:-1]):
            dG += self._paraG[rep[i], rep[i+1]]
        return float(dG)


class AttachedPrimer(Primer):
    '''A primer attached to a template, inherited from Primer
    Args: 
        seq (str): primer sequence (ambiguous bases are supported)
        temperature (float): annealing temperature in Degree Celsius [60]
        sodium (float): the sum of the concentrations of monovalent ions (Na+, K+, NH4+), in molar [0.18]
        template (str): template sequence to be attached (should only consist of A, T, C, G)
    '''
    def __init__(self, seq: str, temperature: float, sodium: float, template: str) -> None:
        # initialize with Primer
        super().__init__(seq, temperature, sodium)
        seq = self.__str__()

        # find where the primer is aligned on the template (local alignment)
        # primer must be the target (first arg), so that its sequence is unchanged
        template = template.upper()
        aln_fw = aligner_local.align(seq, template)[0] # align to forward strand
        aln_rv = aligner_local.align(seq, template, strand='-')[0] # align to reverse strand
        alignment = max(aln_fw, aln_rv, key=lambda x: x.score) # choose the best strand
        # primer (target) might not aligned head to tail
        left_overhang = int(alignment.coordinates[0][0])
        right_overhang = len(seq) - int(alignment.coordinates[0][-1])
        # get aligned coords of template (query)
        c1 = int(alignment.coordinates[1][0])
        c2 = int(alignment.coordinates[1][-1])
        if c1 <= c2:
            # forward strand
            self.strand = '+'
            self.start = c1 - left_overhang
            self.stop = c2 + right_overhang
            self.template = Seq(alignment.query[self.start:self.stop])
        else:
            # reverse strand
            self.strand = '-'
            self.start = c2 - right_overhang
            self.stop = c1 + left_overhang
            self.template = Seq(alignment.query[self.start:self.stop]).reverse_complement()
        
        # since local alignment might have overhangs, align again with global alignment
        aln_global = aligner_global.align(seq, self.template)[0]
        self.alignment = aln_global
        # get the formated string for printing (self.printed)
        self._format()

    def _format(self) -> str:
        # get the gapped alignment strings (with '-'s)
        primer, template = self.alignment
        middle_str = ''
        # find mismatches and gaps
        for b1, b2 in zip(primer, template):
            if b1 == '-' or b2 == '-':
                middle_str += '-'
            elif b1 == b2:
                middle_str += '|'
            else:
                # check if b1 is an ambiguous bases (we only expect the primer to have ambiguous bases)
                try:
                    amb = AMBREP[b1]
                    if b2 in amb:
                        middle_str += '|'
                    else:
                        middle_str += '.'
                except KeyError:
                    middle_str += '.'
        self._middle_str = middle_str

        # generate the formated string
        if self.strand == '+':
            start = self.start
            stop = self.stop
        else:
            start = self.stop
            stop = self.start
        padding_spaces = ' '*(len(str(start))-1) # align the lines for printing
        # print_str = f'GC: {100*self.GC:.2f}%\n'
        # print_str += f'dG ({self.temperature}{chr(176)}C): {self.dG:.2f} kcal/mol\n'
        # print_str += 'alignment:\n'
        print_str  = f'{padding_spaces}0 {primer} {len(primer)}\n'
        print_str += f'{padding_spaces}  {middle_str}\n'
        print_str +=           f'{start} {template} {stop}\n'
        self.printed = print_str


class MSA(object):
    '''Load, process and plot multiple sequence alignment (MSA) in fasta format
    Args: 
        fasta_path (str): path to the MSA file in FASTA format
        n_cpu (int): use multiple threads
    '''
    def __init__(self, fasta_path: str, n_cpu: int=1) -> None:
        self.primers = dict()

        logger.info(f'Loading MSA from {fasta_path} (threads={n_cpu})...')
        tik = time()
        msa = [record.seq for record in SeqIO.parse(fasta_path, 'fasta')]
        # convert to numpy matrix (all characters are converted to uppercase)
        if len(msa) <= n_cpu or n_cpu == 1:
            self.matrix = msa_to_matrix(msa)
        elif n_cpu > 1:
            with multiprocessing.Pool(processes=n_cpu) as pool:
                self.matrix = np.concatenate(
                    pool.map(msa_to_matrix, n_chunk(msa, n_cpu)), 
                    axis=0
                )
        else:
            raise ValueError(f'n_cpu should be a positive integer, but received {n_cpu}')
        self.row, self.col = self.matrix.shape
        
        # Bio.Align class can also read an MSA, and convert to matrix, but slower
        # self.msa = Align.read(fasta_path, 'fasta')
        # self.matrix = np.array(self.msa, dtype='U')

        self._get_consensus(n_cpu)
        print_time_delta(time()-tik)
    
    def _get_consensus(self, n_cpu: int=1) -> None:
        '''Calculate the concensus sequence
        '''
        # count the number of each base for each MSA column
        # returns a matrix with 5 rows (bases) and self.col columns
        logger.info(f' - Generating consensus ({self.row} rows, {self.col} columns)...')
        count_bases_args = [
            (self.matrix[:, i],) # count_bases() only has one argument
            for i in range(self.col)
        ]
        self.vote = np.array(mp_wrapper(
            count_bases, count_bases_args, n_cpu, 
        )).T # transpose the matrix, so that it has 5 rows (bases) and self.col columns
        # find the most frequent base for each column
        consensus = np.array([DNA[i] for i in np.argmax(self.vote, axis=0)])
        
        # get the gapless consensus string, and a mapping of its positions to MSA columns, 
        # so that when a sequence is aligned to the gapless consensus, we can also know where it aligns to the MSA
        consensus_gapless = ''
        consensus2msa = dict() # position in consensus_gapless -> MSA column index
        msa2consensus = dict() # MSA column index -> position in consensus_gapless
        loc = 0 # position on consensus
        for i_col, b in enumerate(consensus):
            msa2consensus[i_col] = loc
            if b != '-':
                consensus_gapless += b
                consensus2msa[loc] = i_col
                loc += 1
        consensus2msa[loc] = i_col
        self.consensus = consensus_gapless
        self._consensus2msa = consensus2msa
        self._msa2consensus = msa2consensus
        self.consensus_array = consensus

    def attach_primer(self, seq: str, temperature: float, sodium: float, name: str|None=None) -> None:
        if name is None:
            name = f'Primer-{len(self.primers)}'
        primer = AttachedPrimer(seq, temperature, sodium, self.consensus)

        # fetch MSA columns of the primer
        msa_slice = self.matrix[:, 
            self._consensus2msa[primer.start]: self._consensus2msa[primer.stop]
        ]
        msa_slice = [''.join(s).replace('-', '') for s in msa_slice]
        if primer.strand == '-':
            msa_slice = [Seq(s).reverse_complement() for s in msa_slice]
        n_perfect_match = len([True for s in msa_slice if s in primer.expanded])
        
        # format print string
        print_str  = f'{name}\n'
        print_str += f'GC: {100*primer.GC:.2f}%\n'
        print_str += f'dG ({primer.temperature}{chr(176)}C): {primer.dG:.2f} kcal/mol\n'
        print_str += f'primer/consensus:\n{primer.printed}'
        print_str += f'sensitivity: {100*n_perfect_match/self.row:.2f}% ({n_perfect_match}/{self.row})\n'
        primer.printed = print_str

        self.primers[name] = primer
        logger.info(f'Successfully attached primer {name}')
        return primer
    
    def variant_call(self) -> dict:
        '''Output location and frequencies of variations using the consensus as reference. 
        '''
        logger.info('Calling variations using the consensus sequence as reference...')
        tik = time()
        msa = self.matrix
        ref = self.consensus_array
        num_seq = self.row

        # because of insertions, locations on the reference is not the same as column number of MSA
        loc = [] # create an array of locations of ref (gaps/insertions do not increase location)
        curr_loc = 1
        last_loc = 1
        is_ins = [] # label insertions (list of T/F)
        ins_loc_len = {} # loc of ins: (start col of ins, length of ins)
        for i_col, base in enumerate(ref):
            if base == '-':
                is_ins.append(True)
                loc.append(last_loc)
                try:
                    ins_loc_len[last_loc][1] += 1
                except KeyError:
                    ins_loc_len[last_loc] = [i_col, 1]
            else:
                is_ins.append(False)
                loc.append(curr_loc)
                last_loc = curr_loc
                curr_loc += 1
        # count the frequency of insertions at each ins_loc
        # freq of ins cannot be calculated the same way as substitutions and deletions
        # sub and del can be calculated column by column, while ins cannot
        ins_loc_freq = {} # loc of ins: freq of ins
        for ins_loc, (ins_col, ins_len) in ins_loc_len.items():
            temp_mtx = msa[:, ins_col:ins_col+ins_len]
            # count the number of rows containing '-'
            ins_num = sum(np.any(temp_mtx!='-', axis=1))
            #if ins_num != 0:
            ins_loc_freq[ins_loc] = ins_num/num_seq

        # variant calling of sub and del
        var_dict = {} # loc: freq
        for i_col, (curr_loc, ins_flag) in enumerate(zip(loc, is_ins)):
            if ins_flag:
                continue
            curr_col = msa[:, i_col]
            curr_freq = sum(curr_col!=ref[i_col])/num_seq
            if curr_freq != 0:
                var_dict[curr_loc] = curr_freq

        # add frequencies of ins
        for ins_loc, ins_freq in ins_loc_freq.items():
            try:
                var_dict[ins_loc] += ins_freq
            except KeyError:
                var_dict[ins_loc] = ins_freq
        print_time_delta(time()-tik)
        return var_dict
    
    def plot(self, save_path: str|None=None) -> None:
        fig = go.Figure()
        # plot MSA
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 'lightgray']
        x = np.array(range(self.col + 1))
        # x_plot: [0, 1, 1, 2, 2, 3, 3, 4, 4, ..., self.col-1, self.col-1, self.col]
        x_plot = np.repeat(x, 2)[1:-1]
        for i, base in enumerate(DNA):
            freq_plot = np.repeat(self.vote[i, :], 2)
            fig.add_scatter(
                x=x_plot, 
                y=freq_plot, 
                hoverinfo='skip', 
                name=base, 
                mode='lines', 
                line=dict(width=0, color=colors[i]), 
                stackgroup='one'
            )
        
        # plot primers
        all_primer_x = []
        all_primer_y = []
        all_hover_text = []
        all_mismatch_x = []
        all_mismatch_y = []
        primer_y_prev = self.row/2
        for primer in self.primers.values():
            if primer.strand == '+':
                primer_x = [
                    self._consensus2msa[primer.start], 
                    self._consensus2msa[primer.stop], 
                    None # add None so that primers will not be connected
                ]
            else:
                primer_x = [
                    self._consensus2msa[primer.stop], 
                    self._consensus2msa[primer.start], 
                    None
                ]
            all_primer_x.extend(primer_x)
            # place randomly on y-axis
            primer_y = rng.uniform(self.row*0.25, self.row*0.75)
            while abs(primer_y-primer_y_prev) < self.row*0.1:
                primer_y = rng.uniform(self.row*0.25, self.row*0.75)
            primer_y_prev = primer_y
            all_primer_y.extend([primer_y, primer_y, None])
            # hover text should be html style
            #hover_text = f'{primer_name}\n{primer.printed}'
            hover_text = primer.printed.replace('\n', '<br>')
            all_hover_text.extend([hover_text]*3)

            # lable mismatches and gaps (gaps pending)
            if primer.strand == '+':
                middle_str = primer._middle_str
                template = primer.alignment[1]
            else:
                middle_str = primer._middle_str[::-1]
                template = primer.alignment[1][::-1]
            curr_consensus_pos = primer.start
            mismatch_x = []
            for i, s in enumerate(middle_str):
                if s == '|':
                    curr_consensus_pos += 1
                    continue
                elif s == '.':
                    mismatch_x.append(self._consensus2msa[curr_consensus_pos]+0.5)
                    curr_consensus_pos += 1
                else: # gap
                    if template[i] == '-':
                        pass # pending
            mismatch_y = [primer_y]*len(mismatch_x)
            all_mismatch_x.extend(mismatch_x)
            all_mismatch_y.extend(mismatch_y)

        fig.add_scatter(
            x=all_primer_x, 
            y=all_primer_y, 
            showlegend=False, 
            name='', 
            marker=dict(size=12, symbol='arrow-bar-up', angleref='previous'), 
            line=dict(width=5, color='black'), 
            connectgaps=False, 
            hovertemplate='%{text}', 
            text=all_hover_text
        )
        fig.add_scatter(
            x=all_mismatch_x, 
            y=all_mismatch_y, 
            name='mismatch to consensus', 
            mode='markers', 
            marker=dict(
                size=12, 
                symbol='x-thin', 
                line=dict(width=3, color='blueviolet')
            ), 
            hoverinfo='skip'
        )

        # title and axis
        fig.update_layout(
            # hover label
            hoverlabel=dict(
                bgcolor='white', 
                font_family='Overpass, monospace' # use mono font to show alignment correctly
            ), 

            # title
            title=dict(
                text=f'MSA of {self.row} sequences ({self.col} columns); length of consensus: {len(self.consensus)}', 
                #text=f'MSA of {self.row} H1N1 HA sequences (GISAID 04/25/24-10/10/24); length of consensus: {len(self.consensus)}', 
                xanchor='left', 
                x=0.06, 
                # yanchor='bottom', 
                # y=0.8, 
                font=dict(size=18)
            ), 

            # legend
            legend=dict(
                orientation='h', 
                xanchor='right', 
                x=1, 
                yanchor='bottom', 
                y=1, 
                traceorder='normal' # A, T, C, G are shown left to right
            ), 

            # axis
            xaxis=dict(
                # setting customized tick lables will show all lables at all times, zoom-in or not
                # tickmode='array', 
                # tickvals=list(self._msa2consensus), 
                # ticktext=list(self._msa2consensus.values()), 
                tickfont=dict(
                    size=16, 
                ), 
                showgrid=False, 
                rangeslider=dict(
                    visible=True
                )
            ), 
            yaxis=dict(
                range=[0, self.row], 
                tickfont=dict(
                    size=16, 
                ), 
                showgrid=False
            ), 

            # axis title
            xaxis_title='MSA position (primer coordinates are 0-start, half-open)', 
            yaxis_title='# bases', 

            # global font
            font=dict(
                size=18, 
                family='Arial'
            )
        )
        if save_path is None:
            fig.show()
        else:
            if save_path[-5:] != '.html':
                save_path += '.html'
            fig.write_html(f'{save_path}')


def run_variant_call(msa_path: str, msa_filename: str, prefix: str|None=None, n_cpu: int=1, min_var: float=0.01) -> None:
    msa = MSA(msa_path, n_cpu)
    var_dict = msa.variant_call()

    if prefix is None:
        prefix = msa_path

    if not os.path.exists(prefix):
        os.makedirs(prefix)
    
    consens_path = os.path.join(prefix, f'{msa_filename}_consensus.fasta')
    SeqIO.write(
        SeqRecord(Seq(msa.consensus), id=msa_filename, description='consensus sequence'), 
        consens_path, format='fasta'
    )

    var_path = os.path.join(prefix, f'{msa_filename}_var.csv')
    with open(var_path, 'w') as f:
        f.write('START,STOP,FREQ\n')
        for pos, freq in var_dict.items():
            if freq >= min_var:
                f.write(f'{pos},{pos},{freq}\n')
    return consens_path, var_path

def run_validate(msa_path: str, primer_pool: str, pool: int, out_path: str, title: str, 
    temperature: float=60.0, sodium: float=0.18, n_cpu: int=1) -> None:
    
    msa = MSA(msa_path, n_cpu)
    consensus_record = SeqIO.SeqRecord(Seq(msa.consensus), id='consensus', description=f'of {msa.row} sequences')

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

    out_text = ''
    for i in range(len(seq_list)):
        primer = msa.attach_primer(seq_list[i], temperature, sodium, seq_names[i])
        out_text += primer.printed + '\n'
    with open(os.path.join(out_path, f'{title}_pool-{pool}.out'), 'w') as f:
        f.write(out_text)
        SeqIO.write(consensus_record, f, 'fasta')
    msa.plot(os.path.join(out_path, title + f'_pool-{pool}.html'))

def run_preprocess(msa_path: str, msa_filename: str, prefix: str|None=None, n_cpu: int=1, min_var: float=0.01):
    if prefix is None:
        prefix = msa_path

    logger.info("Running variant calling to extract consensus and SNPs...")
    consens_path, var_path = run_variant_call(msa_path, msa_filename, prefix, n_cpu, min_var)

    return consens_path, var_path


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--version', '-v', action='version', version='%(prog)s v' + __version__)
    subparsers = parser.add_subparsers(dest='subparser_name', help='sub-commands')

    # variant call
    snps_help = 'Input an MSA, output the consensus sequence, as well as a list of the location and frequencies of variations/SNPs, using the consensus as reference'
    snps_parser = subparsers.add_parser('snps', help=snps_help, description=snps_help)
    snps_parser.add_argument(
        'msa_path', type=str, metavar='msa-fasta', 
        help='Path to the MSA file in FASTA format.')
    snps_parser.add_argument(
        '--prefix', '-o', type=str, default=None, metavar='<string>', 
        help='Prefix for output files (.consensus.fasta and .csv). If not provided, use the MSA file path.')
    snps_parser.add_argument(
        '--threads', '-p', type=int, default=1, metavar='<int>', 
        help='Number of threads [1].')

    # validate
    validate_help = 'Visualize an MSA, along with its primers/probes (if provided), and validate their alignment and sensitivity.'
    validate_parser = subparsers.add_parser(
        'validate', help=validate_help, description=validate_help, formatter_class=argparse.RawTextHelpFormatter
    )
    validate_parser.add_argument(
        'msa_path', type=str, metavar='msa-fasta', 
        help='Path to the MSA file in FASTA format.')
    validate_parser.add_argument(
        'primer_pool', type=str, metavar='csv-file', 
        help='Path to the csv file of a primer pool. Required columns: "amp_id" (amplicon name), "fP" (sequence of forward primer), "rP" (sequence of reverse primer).')
    validate_parser.add_argument(
        '--pool', type=int, default=1, metavar='<int>', 
        help='Primer pool number [1].')
    validate_parser.add_argument(
        '--temperature', '-t', type=float, default=60, metavar='<float>', 
        help='Annealing temperature in Degree Celsius [60.0].')
    validate_parser.add_argument(
        '--sodium', '-s', type=float, default=0.18, metavar='<float>', 
        help='The sum of the concentrations of monovalent ions (Na+, K+, NH4+), in molar [0.18].')
    validate_parser.add_argument(
        '--output', '-o', type=str, default='./', metavar='<string>', 
        help='Output path (output to current directory by default).')
    validate_parser.add_argument(
        '--title', '-t', type=str, default='olivar-val', metavar='<string>', 
        help='Name of validation [olivar-val].')
    validate_parser.add_argument(
        '--threads', '-p', type=int, default=1, metavar='<int>', 
        help='Number of threads [1].')

    args = parser.parse_args()

    if args.subparser_name == 'snps':
        run_variant_call(
            msa_path=args.msa_path, 
            msa_filename=args.msa_filename, 
            prefix=args.prefix, 
            n_cpu=args.threads
        )
    elif args.subparser_name == 'validate':
        run_validate(
            msa_path=args.msa_path, 
            primers_pool=args.primer_pool, 
            pool=args.pool,
            temperature=args.temperature, 
            sodium=args.sodium, 
            out_path=args.output, 
            title=args.title,
            n_cpu=args.threads
        )
    else:
        print('A sub-command needs to be specified (snps or validate).')
        print("Use '--help' to print detailed descriptions of command line arguments")
