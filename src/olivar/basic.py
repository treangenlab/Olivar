#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Some basic tools for processing DNA sequences and coordinates
'''


__author__ = 'Michael X. Wang'


import os
CURRENT_DIR = os.path.dirname(__file__)


from collections import Counter

import numpy as np

import logging
logger = logging.getLogger('main')

IUPAC_CODES = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
}

complement = {
            'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N',
            'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n',
            'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K',
            'B':'V', 'V':'B', 'D':'H', 'H':'D',
            'r':'y', 'y':'r', 's':'s', 'w':'w', 'k':'m', 'm':'k',
            'b':'v', 'v':'b', 'd':'h', 'h':'d',
            '0':'1', '1':'0', '2':'3', '3':'2', '4':'4',
            0:1, 1:0, 2:3, 3:2, 4:4
}

degenerate_bases = {
        frozenset(['A', 'C']): 'M',
        frozenset(['A', 'G']): 'R',
        frozenset(['A', 'T']): 'W',
        frozenset(['C', 'G']): 'S',
        frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'T']): 'K',
        frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'T']): 'H',
        frozenset(['A', 'G', 'T']): 'D',
        frozenset(['C', 'G', 'T']): 'B',
        frozenset(['A', 'C', 'G', 'T']): 'N',

        frozenset(['a', 'c']): 'm',
        frozenset(['a', 'g']): 'r',
        frozenset(['a', 't']): 'w',
        frozenset(['c', 'g']): 's',
        frozenset(['c', 't']): 'y',
        frozenset(['g', 't']): 'k',
        frozenset(['a', 'c', 'g']): 'v',
        frozenset(['a', 'c', 't']): 'h',
        frozenset(['a', 'g', 't']): 'd',
        frozenset(['c', 'g', 't']): 'b',
        frozenset(['a', 'c', 'g', 't']): 'n'
}

def revcomp(seq):
    '''
    give the reverse complement of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        reverse complement of input
    '''
    if isinstance(seq, str) or isinstance(seq, list):
        try:
            bases = [complement[base] for base in seq]
        except KeyError:
            raise ValueError(f"Unrecoginized bases are found in sequence.")
        return ''.join(reversed(bases)) if isinstance(seq, str) else list(reversed(bases))
    else:
        raise ValueError('Only string or list input is accepted for reverse complement.')

def get_GC(seq):
    '''
    get the GC ratio of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        GC ratio of input sequence
    '''
    return len([base for base in list(seq) if base in ['G','g',3,'C','c',2]])/len(seq)


def seq2arr(seq_str):
    '''
    convert string sequence to array (list of numbers)
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input: 
        string sequence
    output:
        array sequence
    '''
    dic = {'A':0, 'T':1, 'C':2, 'G':3, 'a':0, 't':1, 'c':2, 'g':3}
    try:
        seq_list = [dic[s] for s in list(seq_str)]
    except KeyError:
        raise ValueError(f"Base(s) other than 'A', 'T', 'C', 'G' is found in '{seq_str[:10]}...', ambiguous bases are not accepted.")
    return seq_list


def arr2seq(seq_array):
    '''
    convert array to string sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input: 
        array sequence
    output:
        string sequence
    '''
    dic = {0:'A', 1:'T', 2:'C', 3:'G'}
    seq_list = [dic[s] for s in seq_array]
    return ''.join(seq_list)


def seq2num(seq):
    '''
    use quaternary to convert array or string sequence to a number
    CAUTION: diffent sequence may output the same number (e.g. 'AACT' and 'AAACT')
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        an integer
    '''
    bases = list(seq)
    dic = {'A':0, 'T':1, 'C':2, 'G':3, 
    'a':0, 't':1, 'c':2, 'g':3, 
    0:0, 1:1, 2:2, 3:3}
    bases = [dic[base] for base in bases]
    len_bases = len(seq)
    out_num = 0
    for i in range(len_bases):
        out_num = out_num + 4**(len_bases-i-1)*bases[i]
    return out_num


def num2seq(num, l=None):
    '''
    The reverse function of seq2num
    CAUTION: a single number represents multiple sequences if l=None
    input:
        num: int
        l: length of output sequence; if none, 'A' at sequence beginning is ignored
    output:
        sequence
    '''
    base = ['A', 'T', 'C', 'G']
    seq = []
    if l is None:
        while True:
            new_num = num // 4
            b = num % 4
            seq.append(base[b])
            if new_num == 0:
                break
            num = new_num
        return ''.join(seq[::-1])
    elif num < 4**l:
        for _ in range(l):
            new_num = num // 4
            b = num % 4
            seq.append(base[b])
            num = new_num
        return ''.join(seq[::-1])
    else:
        logger.info('Input number out of range according to input sequence length.')
        return None


def seq2complex(seq_str):
    '''
    Convert string sequence to a list of complex
    '''
    dic = {'A':1, 'T':-1, 'C':1j, 'G':-1j, 'N': 0, 'a':1, 't':-1, 'c':1j, 'g':-1j, 'n': 0}
    seq_list = [dic[s] for s in list(seq_str)]
    return seq_list


def complex2seq(seq_complex):
    '''
    Convert a list of complex back to string
    '''
    seq_str = []
    for c in seq_complex:
        if c.real > c.imag:
            if -c.real < c.imag:
                seq_str.append('A')
            elif -c.real > c.imag:
                seq_str.append('G')
            else:
                seq_str.append('N')
        elif c.real < c.imag:
            if -c.real < c.imag:
                seq_str.append('C')
            elif -c.real > c.imag:
                seq_str.append('T')
            else:
                seq_str.append('N')
        else:
            seq_str.append('N')
    return ''.join(seq_str)


def randseq(l, seed=None):
    '''
    Generate random sequence of length l
    input:
        l: sequence length
        seed: random seed for numpy
    output:
        random sequence of length l
    '''
    if seed is not None:
        np.random.seed(seed)
    return arr2seq(np.random.randint(4, size=l).tolist())


def randarr(l, seed=None):
    '''
    Generate random array of length l
    input:
        l: array length
        seed: random seed for numpy
    output:
        random array of length l
    '''
    if seed is not None:
        np.random.seed(seed)
    return np.random.randint(4, size=l).tolist()


def get_complexity(seq):
    # calculate Shannon entropy for word size of 1, 2 and 3
    hash1 = Counter()
    hash2 = Counter()
    hash3 = Counter()

    # count occurance
    for j in range(len(seq)-2):
        hash1[seq[j]] += 1
        hash2[seq[j:j+2]] += 1
        hash3[seq[j:j+3]] += 1
    hash1[seq[-2]] += 1
    hash1[seq[-1]] += 1
    hash2[seq[-2:]] += 1

    # calculate probability mass function
    p1 = np.array(list(hash1.values()))/(len(seq))
    p2 = np.array(list(hash2.values()))/(len(seq)-1)
    p3 = np.array(list(hash3.values()))/(len(seq)-2)

    SE1 = -np.sum(p1 * np.log2(p1))
    SE2 = -np.sum(p2 * np.log2(p2))
    SE3 = -np.sum(p3 * np.log2(p3))
    return min([SE1/2, SE2/4, SE3/6])

def get_combinations(seq):
    """Calculate the total number of possible sequence combinations from a degenerate sequence."""
    score = 1
    for char in seq:
        score *= len(IUPAC_CODES.get(char, [char]))  # Default to 1 if char is standard
    return score

def get_degenerate_base(bases):
    """Return the IUPAC degenerate base code for a set of nucleotides"""
    if isinstance(bases, str):
        unique_bases = set(bases)
    else:
        unique_bases = {b for b in bases}
    
    if len(unique_bases) == 1:
        return unique_bases.pop()
    
    return degenerate_bases.get(frozenset(unique_bases), 'N')
