#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module of utility functions for CRISPR related requests to BioRPC services
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import json
import os
import logging
import pprint

import requests

log = logging.getLogger(__name__)

BIORPC_URL = 'rpc'

# FUTURE - make on and offtarget scores filters
GUIDE_FILTERS = ['filter_gc_content', 'filter_consecutive_bases',
                  'filter_pam_sites', 'filter_restriction_site']


def make_slice_instruction(chr_name, start_end, nuclease_name):
    """
    Make json RPC instruction for loading nuclease site by genomic coordinates
    :param chr_name: chromosome name
    :param start_end:
    :param nuclease_name:
    :return:
    """
    params = {
        "chromosome": {
            'name': chr_name,
        },
        "start_end": (start_end),
        "nuclease": {'name': nuclease_name},
        'async': False,
    }

    return params


def slice_nuclease_track(target_server, chr_name, start_end, nuclease_name):
    """
    Load the nuclease sites for a specific interval of genomic coordinates
    :return:
    """
    target_url = os.path.join(target_server, BIORPC_URL)
    #Future - update server-side code to take a genome argument
    #Future - allow a list of several nucleases
    body = make_slice_instruction(chr_name, start_end, nuclease_name)

    resp = requests.post(target_url,
                         json.dumps(body),
                         headers={b'content-type': b'application/json'})
    try:
        return resp.json()['result']
    except KeyError:
        log.warn("No RESULT in response")
        pprint.pprint(resp.json())
        raise
    except ValueError:
        pprint.pprint(resp.text)
        raise


def make_score_guide_instruction(guides, genome_dict, activity_params, offtarget_params):
    """
    Utility to build the body of a batch guide analysis request.
    This is a simple version of the method for small batches. make_scoring_job
    is for larger tasks like libraries.

    :param target_server:
    :param guide:
    :param genome_dict:
    :param sf_dict:
    :return:
    """
    for id_, g in enumerate(guides):
        yield {
            'guide': g,  #FUTURE - validate guide
            'genome': genome_dict,  # constant
            'scoring_function': {
                'doench': activity_params,  # rename to activity
                'mitv1_db': offtarget_params  # rename to offtarget
            },
        }


def walk_gene(gene, genome, async):
    """Test you can precompute guides for a gene"""
    target_url = RPC_URL
    initial_request = {
        'jsonrpc': '2.0',
        'method': 'walk_gene',
        'params': {
            "genome": {'version': genome},
            'gene': {'name': gene},
            "nuclease": {'name': 'wtCas9'},
            'scoring_function': {
                # 'mitv1_db': {
                #     'maxmismatch': 3,
                #     'minscore': 0.70,
                #     'offset': 0
                # },
                'doench': {'minscore': 0.15},
                #'filters': None,
                'filters': [
                    'filter_gc_content',
                    'filter_consecutive_bases',
                    'filter_pam_sites',
                    'filter_restriction_site',
                ],
            },
        "async": async,
        },
        'id': 1,
    }
    resp = requests.post(target_url,
                         json.dumps(initial_request),
                         headers={b'content-type': b'application/json'},
                         )
    try:
        guides = resp.json()['result']
        click.echo(pprint.pprint(guides))
        return guides
    except KeyError:
        log.warn("No RESULT in response")
        pprint.pprint(resp.json())
    except ValueError:
        pprint.pprint(resp.text)



def run_pair_tornado(gene, genome, nuclease, async, output):
    nuclease_dict = {'name': nuclease}
    target_url = RPC_URL
    initial_request = {
        'jsonrpc': '2.0',
        'method': 'pair_tornado',
        'params': {
            'gene': {'name': gene},
            "genome": {'version': genome},
            "nuclease": nuclease_dict,
            'scoring_function': {
                # 'mitv1_db': {
                #     'maxmismatch': 3,
                #     'minscore': 0.50,
                #     'offset': 0
                # },
                # 'filters': [
                #     'filter_gc_content',
                #     'filter_consecutive_bases',
                #     'filter_pam_sites',
                #     'filter_restriction_site',
                # ],
                'doench': {'minscore': 0.15},
                'filters': None,
            'focal_point': 'n-terminal',
            },
        "async": async,
        },
        'id': 1,
    }
    resp = requests.post(target_url,
                         json.dumps(initial_request),
                         headers={b'content-type': b'application/json'},
                         )
    try:
        if async:
            task_id = resp.json()['result']['task_id']
            return task_id
        else:
            yaml.dump(resp.json()['result'], output)
    except KeyError:
        click.echo(resp.text)
        raise
    except ValueError:
        click.echo(resp.text)
        raise