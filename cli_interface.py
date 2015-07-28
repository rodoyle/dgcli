#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Command line interface to the DeskGen Genome Editing Platform.
"""

from __future__ import unicode_literals

import json
import logging
import os
import pprint
import sys

import click
import requests
import xlsxwriter

import dgparse

import dgcli.genomebrowser as gb
from dgcli.genomebrowser import GB_MODELS
from dgcli import genome_editing as ge
from dgcli import utils, libraries, enevolv_writer
from dgcli import inventory as inv
from dgcli.configuration import load


logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

RC_PATH = os.path.expanduser('~/.dgrc')  # Get normalized location of RC files


@click.group()
@click.pass_context
def cli(ctx):
    config = load(RC_PATH)
    ctx.obj = config


@cli.command('login')
@click.pass_context
@click.option('--email')
@click.option('--password')
def login_cmd(ctx, email, password):
    try:
        ctx.obj.target_server.inventory.login(email, password)
    except inv.AuthenticationError:
        click.echo("Error Authenticating", err=True)


@cli.command('read')
@click.argument('type')
@click.option('--filters', default=None)
@click.option('--output', type=click.File(), default=sys.stdout)
@click.pass_context
def read_cmd(ctx, type_, filters, output):
    """Load all of DNA Molecule, Design, Feature, or Annotation from a
    remote source"""
    filters = filters if filters else {}
    if object in GB_MODELS:
        service = 'genomebrowser'
        body = gb.make_fetch_instruction(type_, filters)
    else:
        service = 'inventory'
        body = {'object': type_, 'filters': filters}

    endpoint_url = os.path.join(ctx.target_server.url, 'api/{0}/crud'.format(service))
    credentials = ctx.user.credentials
    resp = requests.post(endpoint_url, json.dumps(body), auth=credentials)
    try:
        read_list = resp.json()['read']
        json.dump(read_list, output, indent=2)
    except KeyError:
        pprint.pprint(resp.text)
    except ValueError:
        pprint.pprint(resp.text)


@cli.command('create')
@click.argument('record_files', type=click.STRING, nargs=-1)
@click.option('--record_type')
@click.option('--format', type=click.Choice(['json', 'xlsx', 'csv', 'yaml']))
@click.pass_context
def create_cmd(ctx, record_files, record_type, format):
    """Create a record on a remote server from local data or fail"""
    # Trigger parsing and generate list of json records
    records = utils.iterate_records_from_files(record_files)
    for record in records:
        data, errors = ctx.obj.target_server.inventory.create(record)
        if errors:
            for k, error_message in errors.iteritems():
                click.echo(error_message, err=True)


@cli.command('slice_genome')
@click.argument('chromosome')
@click.argument('start_end', nargs=2, type=click.IntRange())
@click.argument('tracks', type=click.Choice(GB_MODELS))
@click.option('--genome')
@click.option('--strand', default=0)
@click.option('--output', type=click.File('wb'), default=sys.stdout)
@click.option('--sequence', default=False)
@click.option('--nuclease', default='wtCas9')
@click.pass_context
def slice_genome_cmd(ctx, genome, chromosome, start_end, tracks, strand, output,
                     sequence, nuclease):
    """"
    Load all the track data in the genome browser in an slice interval from
    a remote source.
    """
    # FUTURE - genomebrowser should just support dynamic tracks for guides,
    # talens, and other features. For now we use this work around.
    # FUTURE - allow other intersections besides overlaps
    guides = False
    if 'guides' in tracks:
        guides = True
        del tracks[tracks.index('guides')]
    ctx.target_server.genomebrowser.slice(genome, chromosome, start_end, strand, sequence)  # NOQA
    # this is a work around as RPC does the guides while GB does the rest
    if guides:
        nuc_slice = ge.make_slice_instruction(chromosome, start_end, nuclease)
        utils.make_json_rpc(0, 'load_guides', nuc_slice, output)


@cli.command('score_guides')
@click.argument('guides', type=click.File(), default=sys.stdin)
@click.option('--genome')
@click.option('--output', '-o', type=click.File(), default=sys.stdout)
@click.option('--async', default=False, type=click.BOOL)
@click.option('--activity')
@click.option('--offtarget')
@click.pass_context
def score_guides_cmd(ctx, guides, genome, async, activity, offtarget, output):
    """
    Analyze a batch of guide RNAs
    :return:
    """
    ctx.target_server.workers.score_guide(guides, genome, activity, offtarget)


@cli.command()
@click.argument('spec_file_path', type=click.Path())
@click.option('--target_file', default=sys.stdin, type=click.File())
@click.option('--output', default=sys.stdout, type=click.File())
@click.option('--update', default=False, type=click.BOOL)
def append_targets(spec_file_path, target_file, output, update):
    """Parse an xlsx or csv formatted target list
    Append a new target JSON object to the library spec JSON object
    """
    with open(spec_file_path, b'r') as spec_file:
        click.echo(spec_file)
        specs = json.load(spec_file)
    if target_file.name.endswith('xlsx'):  # if excel file
        specs['targets'] = libraries.parse_excel_file(target_file.name)
    if update:
        with open(spec_file_path, b'w') as spec_file:
            json.dump(specs, spec_file, indent=2)
    else:
        json.dump(specs, output, indent=2)


@cli.command('extract')
@click.argument('source', type=click.STRING, nargs=-1)
@click.option('--output', default='extracted_dnafeatures.xlsx')
@click.option('--debug', type=click.BOOL, default=False)
def extract_cmd(source, output, debug):
    """
    Extract all instances of a model from a set of files matching a pattern.
    :param object_type:
    :returnj:
    """
    workbook = xlsxwriter.Workbook(output)
    records = []
    for fname in source:
        msg = "Parsing file {0}".format(fname)
        click.echo(msg)
        with open(fname, 'r') as snapfile:
            try:
                # parse the snapgene file and adapt to dtg's schema
                adapted_data = dgparse.snapgene.parse(snapfile)
                if adapted_data['dnafeatures']:
                    # update the record with data from the filename
                    adapted_data.update(inv.parse_file_name(fname))
                    records.append(adapted_data)
                else:
                    click.echo("{0} didn't contain features".format(os.path.basename(fname)))
                if output == "stdio":
                    click.echo(adapted_data)
            except dgparse.exc.ParserException:
                click.echo("Error Parsing {0}".format(fname), err=True)
                if debug:
                    raise
    if output != "stdio":
        enevolv_writer.write_to_xls(workbook, records)

if __name__ == '__main__':
    cli(obj={})
