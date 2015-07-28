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

import dgparse


from dgcli import async_job

log = logging.getLogger(__name__)


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


def iterate_records_from_files(record_files):
    """
    Iterate over all the records in a collection of files
    :return:
    """
    records = []
    for record_path in record_files:
        path, format = os.path.splitext(record_path)
        parser = dgparse.PARSERS[format]
        with open(record_path, 'r') as record_file:
            records.extend(parser(record_file))
    return records

