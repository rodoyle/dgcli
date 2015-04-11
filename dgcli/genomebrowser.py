#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Library module of requests to genomebrowser service endpoints.

This module is separate from the CLI interface so that the underlying functions
may be imported and used in scripts or iPython sessions without the click
decorators.

The aim of each function is to fetch the designed object, handle any arising
errors, and return a parsed JSON version that may be passed to any constructor
or "rehydrator" function.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os

import json
import requests

from logging import log

GB_MODELS = ['gene', 'transcript', 'exon', 'cds', 'cdsregion', 'breakpoint',
             'translocation', 'fragilesite', 'guide', 'trackedguide']
GENE_EXPANSION = [
    ['coding_sequences', 'regions', 'coordinates'],
    ['transcripts', 'exons', 'coordinates'],
    ['coordinates']
]
TRANSCRIPT_EXPANSION = [['exons', 'coordinates'], ['coordinates']]
CDS_EXPANSION = [['regions', 'coordinates'], ['coordinates']]
BREAKPOINT_EXPANSION = [['gene', 'coordinates'], ['coordinates']]
GB_CRUD_URL = 'api/genomebrowser/crud'


def make_fetch_instruction(object_type, filter_attrs, count=False):
    """
    Make an instruction for fetching the details of the target object.
    :param annotation_dict:
    :return:
    """

    if object_type == 'gene':
        expansion = GENE_EXPANSION
    elif object_type == "transcript":
        expansion = TRANSCRIPT_EXPANSION
    elif object_type == any(('fragilesite', 'breakpoint')):
        expansion = BREAKPOINT_EXPANSION
    elif object_type == 'genome':
        expansion = [['chromosomes']]
    elif object_type == 'chromosomes':
        expansion = [['genome']]
    else:
        expansion = [['coordinates']]

    return {
        'object': object_type,
        'read': {
            'filter': filter_attrs,
            'expand': expansion
        }
    }


def make_slice_instruction(genome, chromosome, start_end, tracks, strand, sequence):
    """
    Make a valid slice instruction to get features from a genome
    :return:
    """
    genome_instruction = {
        "chromosome": {
            "name": chromosome,
        },
        "genome": {
            "version": genome,
        },
        "target_region": {
            "start": start_end[0],
            "end": start_end[1],
            "strand": strand if strand else 0,  # optional, default 0
        },
        "region_type": tracks,
        "sequence": sequence
    }


def fetch_annotation(target_server, annotation_dict):
    """Load the detail and child objects of an annotation from an annotated
     genome.
    The annotation dict should contain enough information to specify your gene.
    Example: {'name': 'BRCA2', 'accession': 'ENSG....'}
    """
    url = os.path.join(target_server, GB_CRUD_URL)
    resp = requests.post(
        url,
        json.dumps(make_fetch_instruction(annotation_dict)),
        headers={b'content-type': b'application/json'}
    )
    try:
        return resp.json()['read'][0]

    except KeyError:
        msg = "Improper Query to GenomeBrowser"
        log.error(msg)
        log.error(resp.text)
        raise
    except IndexError:
        msg = "No annotations found matching filter {0}".format(annotation_dict)
        log.warn(msg)
        log.error(resp.text)
        raise
