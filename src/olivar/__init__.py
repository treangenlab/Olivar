#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Olivar: multiplex PCR tiling design (https://www.nature.com/articles/s41467-024-49957-9)
Olivar implements a novel algorithm to reduce non-specific amplifications in PCR, 
while avoiding primer design at SNPs and other undesired regions at the same time. 
Olivar also optimize for primer dimers with the SADDLE algorithm (https://www.nature.com/articles/s41467-022-29500-4). 
'''


__author__ = 'Michael X. Wang'
__version__ = '1.3.0'


import os
import sys
CURRENT_DIR = os.path.dirname(__file__)
sys.path.append(CURRENT_DIR)

from main import build, tiling, save, specificity, sensitivity
