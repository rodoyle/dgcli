#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utility and comms functions.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import json
import logging
import multiprocessing
from collections import deque

import os
import requests
import yaml


from dgcli import async_job

log = logging.getLogger(__name__)

def make_post(endpoint_url, credentials, body_dict, output, on_error):
    """
    Utility function to make post requests
    :param url:
    :param body_dict:
    :return:
    """

    resp = requests.post(endpoint_url, json.dumps(body_dict), auth=credentials)
    result_key = 'read' if 'read' in body_dict else None

    if resp.ok:
        try:
            read_list = resp.json()[result_key] if result_key else resp.json()
            yaml.dump_all(read_list, output, indent=2)
        except KeyError:
            on_error(resp)
        except ValueError:
            on_error(resp)
    else:
        on_error(resp)


def make_file_upload_request(endpoint_url, credentials, file_path, output, on_error):

    resp = requests.post(endpoint_url,
                         files={'inventory_file': open(file_path, 'rb')},
                         auth=credentials)
    if resp.ok:  # expect and HTTPCreated 201 Response
        data = resp.json()
        yaml.dump(data, output, indent=2)
    else:
        on_error(resp)




def make_json_rpc(endpoint_url, credentials, id_, method, params, output, on_error, async=False):
    """
    Utility function to make json RPC requests
    Local to CLI module as it is designed to work with standard
    :param id_:
    :param method:
    :param params:
    :return:
    """
    params.update({'async': async})  # add the async setting
    body = {'jsonrpc': 2.0, 'method': method, 'id': id_, 'params': params}
    resp = requests.post(endpoint_url, json.dumps(body), auth=credentials)

    if resp.ok:
        try:
            result = resp.json()['result']  # JSON RPC has strict grammar to help  # NOQA
            yaml.dump_all(result, output)
        except KeyError:
            #likely an error, return in a parable way
            on_error(resp)
        except ValueError:
            on_error(resp)
    else:
        on_error(resp)


def make_batch_rpc(endpoint_url, credentials, jobs, method, output,
                   on_error, async=False):
    """
    Make a batch request. Unless there is a very good reason to do a sync
    batch (like bulk delete), async is strongly recommended to present timeout.

    :param endpoint_url:
    :param jobs:
    :param method:
    :param params:
    :param output:
    :param on_error:
    :param async:
    :return:
    """
    queue = deque()
    output_path = output.name  # should be a file name

    for id_, j in jobs:  # list or iterable
        queue.append({
            'jsonrpc': 2.0,
            'id': id_,
            'method': method,
            'params': j,
            'async': async
        })

    if async:
        lock = multiprocessing.Lock()
        multiprocessing.Process(target=async_job.start_batch,
                                args=(endpoint_url, credentials, jobs,
                                      output_path, on_error, lock)).start()


def iter_repository(repo_root, extension_mapping):
    """
    Return an iterator over the files in a repository in a depth first manner.
    At each terminal node, check the file extensions. If the extension, matches
    something in the extension mapping, invoke the parser in the extension
    mapping.
    """
    for root, dirs, files in os.walk(repo_root):
        #print root, dirs, files
        for folder in dirs:  # filter out the hidden files
            if (root == "."):
                # print "Ignoring", folder
                dirs.remove(folder)
        for name in files:
            extension = os.path.splitext(name)
            if extension in extension_mapping:
                parser = extension_mapping[extension]
                with open(name, 'r') as open_file:
                    parsed_data = parser(open_file)
                    #not really a url, but useful non the less to get at related
                    # files
                    parsed_data['url'] = os.path.join(root, name)
                yield parsed_data


def write_to_xls(workbook, molecules):
    """"
    Given an array of JSON  objects, write them out to an upload ready XLSX
    File. Column headers will have the full attribute path of the data,
    on record will be stored per line.
    **mandatory child objects will be denormalized into the same row as their
    parent**
    """
    feat_sheet = workbook.add_worksheet('dnafeatures')
    seen_feat = {}
    ft_row = 0
    ft_cols = ['category.name', 'accession', 'name', 'pattern.bases',
               'pattern.sha1', 'description', 'properties']
    col = 0
    # Header
    for col, heading in enumerate(ft_cols):
        feat_sheet.write_string(ft_row, col, heading)
        col += 1
        feat_sheet.write_string(ft_row, col, "Extracted From")

    ft_row += 1

    for mol in molecules:
        annotations = mol.get('dnafeatures')
        for ann in annotations:
            ft = ann.get('dnafeature')
            ft_accession = ft.get('accession')
            definition = "{mol_accession}:{start}..{end}".format(
                mol_accession=mol['accession'],
                start=ann['start'],
                end=ann['end']
            )
            if ft_accession in seen_feat:
                # now check the sequence
                seen_bases = seen_feat[ft_accession]['pattern']['sha1']
                new_bases = ft['pattern']['sha1']
                if new_bases == seen_bases:
                    log.debug("Skipping already-seen {0}".format(ft_accession))
                    continue
                else:
                    msg = "{0} re-defined with new nucleotide pattern in {1}".format(
                        ft_accession, definition)
                    log.warn(msg)

            for col, attr in enumerate(ft_cols):
                parent = ft
                if '.' in attr:
                    attr_path = attr.split('.')
                    for child_name in attr_path:
                        parent = parent[child_name]
                    value = parent
                else:
                    value = ft.get(attr)
                if isinstance(value, basestring):
                    feat_sheet.write_string(ft_row, col, value)
                elif isinstance(value, dict):
                    if value == {}:
                        continue
                    flat_value = json.dumps(value)
                    feat_sheet.write_string(ft_row, col, flat_value)
                else:
                    feat_sheet.write(ft_row, col, value)

            # Extraction Source
            feat_sheet.write_string(ft_row, (col + 1), definition)
            ft_row += 1
            seen_feat[ft_accession] = ft  # stash the feature

