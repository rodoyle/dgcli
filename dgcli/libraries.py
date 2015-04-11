#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for triggering library generation.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import xlrd
import yaml


def parse_excel_file(filehandle):
    """
    Parse an excel file and return a list of targets
    :param filehandle:
    :return:
    """
    targets = []

    target_file = xlrd.open_workbook(filehandle)
    sheet = target_file.sheet_by_index(0)  # get the first sheet
    for row_idx in range(1, sheet.nrows):
        row = sheet.row_values(row_idx)
        targets.append({
            "gene": {
                "name": row[0],  # the name
                "accession": row[1]  # the accession
            }
        })
    return targets


def append_targets(spec_file_path, target_file, output, update):
    """Parse an xlsx or csv formatted target list
    Append a new target JSON object to the library spec JSON object
    """
    with open(spec_file_path, b'r') as spec_file:
        specs = yaml.load(spec_file)
    if target_file.name.endswith('xlsx'):  # if excel file
        specs['targets'] = parse_excel_file(target_file.name)
    if update:
        with open(spec_file_path, b'w') as spec_file:
            yaml.dump(specs, spec_file, indent=2)
    else:
        yaml.dump(specs, output, indent=2)


def design_library(spec_file, async, dryrun, output, dump):
    """
    Make a request to the server to start a library design job
    :param spec_file:
    :param async:
    :param dryrun:
    :param output:
    :return:
    """
    target_url = RPC_URL
    specs = json.load(spec_file)
    request = {
        "jsonrpc": '2.0',
        "method": 'design_library',
        "id": 1,
        "params": {
            "genome": specs.get('genome', defaults.DEFAULT_GENOME),
            "nuclease": specs.get('nuclease', defaults.DEFAULT_NUCLEASE),
            "defaults": specs.get('defaults', defaults.DEFAULT_GENE_WALKER),
            "targets": specs.get('targets'),
            #"callbacks": specs.get('callbacks', None),
            "name": specs.get('name', "Custom_Library"),
            "description": specs.get('description'),
            "dry_run": dryrun,
            "async": async,
        },
    }
    if dump:  # dump the request body
        json.dump(request, sys.stdout, indent=2)
        return
    resp = requests.post(target_url,
                         json.dumps(request),
                         headers={b'content-type': b'application/json'},
                         )
    try:
        if async:
            task_id = resp.json()['result']['task_id']
            return task_id
        else:
            json.dump(resp.json()['result'], output, indent=2)
    except KeyError:
        click.echo(resp.text)
        raise
    except ValueError:
        click.echo(resp.text)
        raise












