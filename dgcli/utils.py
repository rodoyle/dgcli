#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utility and comms functions.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import json
import multiprocessing
from collections import deque

import requests
import yaml

from dgcli import async_job


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




