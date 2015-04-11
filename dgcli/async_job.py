#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for running big jobs asynchronously.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

from collections import deque
import json
import requests
import yaml
import sys
import time


def write_to_output_stream(result, output_path, lock=None):
    """"Util to write results to output stream"""
    if output_path == 'stdout':
        try:
            if lock:
                lock.acquire()
            yaml.dump(result, sys.stdout)
        finally:
            if lock:
                lock.release()
    else:
        with open(output_path, 'w') as output:
            yaml.dump(result, output)


def process_response(resp, output, on_error, lock=None):
    """Loop through task queue, prune out errors, write results out, and
    return any remaining items.
    """
    try:
        result = resp.json()['result']
        if isinstance(result, dict):  # we're probably good?
            if 'task_id' in result: # we have a pending task
                return result
            else:
                write_to_output_stream(result, output, lock)
                return None
        elif result == 'FAILURE':
            on_error(resp)
            return None
        else:
            return result
    except KeyError:
        on_error(resp)
        return None
    except AssertionError:
        on_error(resp)
        return None


def check_if_done(target_url, credentials, task, output, on_error, lock):
    """
    Check the status of an asynchronous Job
    :param task_id:
    :return:
    """
    body = {
        'jsonrpc': '2.0',
        'method': 'check_status',
        'id': task['task_id'],
    }
    # ping the server to see if job is done
    resp = requests.post(target_url,
                         json.dumps(body),
                         headers={b'content-type': b'application/json'},
                         auth=credentials)
    return process_response(resp, output, on_error, lock)


def poll_tasks(endpoint_url, credentials, task_queue, output, on_error, lock):
    """" Master processes that will poll tasks in the background"""
    # set up initial request
    while True:
        if len(task_queue) < 1:
            break
        task = task_queue.popleft()
        pending = check_if_done(endpoint_url, credentials, task, output, on_error,
                                lock)
        if pending:
            task_queue.append(task)
        time.sleep(1)


def start_batch(endpoint_url, credentials, jobs, output_path, on_error, lock):
    """Process a batch"""
    # Start the batch
    resp = requests.post(endpoint_url, json.dumps(jobs), auth=credentials)
    # response is a list this time
    task_queue = deque()
    try:
        assert resp.ok
        group_result = resp.json()
        # some jobs may have failed, while some passed
        for r in group_result:
            pending = process_response(r, output_path, on_error)
            if pending:
                task_queue.append(pending)
        if len(task_queue) > 0:
            poll_tasks(endpoint_url, credentials, task_queue, output_path,
                       on_error, lock)
    except AssertionError:
        on_error(resp)