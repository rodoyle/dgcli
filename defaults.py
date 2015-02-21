#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Contains the default values for various commands
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import functools
import logging
import itertools
import operator

DEFAULT_GENE = {"name": "PSMA4"}

DEFAULT_GENOME = {'version': 'GRCh38.p2'}

DEFAULT_NUCLEASE = {'name': 'wtCas9'}

DEFAULT_FILTERS = [
     'filter_gc_content',
     'filter_consecutive_bases',
     'filter_pam_sites',
     'filter_restriction_site',
]

DEFAULT_SCORING_FUNCTION = {
    # 'mitv1_db': {
    #     'maxmismatch': 3,
    #     'minscore': 0.70,
    #     'offset': 0
    # },
    'doench': {'minscore': 0.20},
    'filters': DEFAULT_FILTERS,
}

DEFAULT_NTERM_PAIRS = {
    'gene': DEFAULT_GENE,
    "genome": DEFAULT_GENOME,
    "nuclease": {"name": 'Cas9D10A'},
    'scoring_function': DEFAULT_SCORING_FUNCTION,
    'focal_point': 'n-terminal',
    "async": False,
}

DEFAULT_GENE_WALKER = {
    'gene': DEFAULT_GENE,
    'genome': DEFAULT_GENOME,
    'nuclease': DEFAULT_NUCLEASE,
    'scoring_function': DEFAULT_SCORING_FUNCTION,
    'async': False,
}
