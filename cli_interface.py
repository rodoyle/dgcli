"""System level test of the guide services
"""

from __future__ import unicode_literals

import functools
import json
import logging
import os
import pprint
import sys
import yaml
import xlrd

import click
import requests

import defaults
import dgcli.genomebrowser as gb
from dgcli.genomebrowser import GB_MODELS
from dgcli import genome_editing as ge
from dgcli import utils, libraries

logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

RC_PATH = os.path.expanduser('~/.dgrc')  # Get normalized location of RC files
with open(RC_PATH, 'r') as rc_file:
    CONFIG = yaml.load(rc_file) # Load Run Ctrl config into RAM

# Check the ENV for any local overrides
CONFIG.update(os.environ)
BIODATA = CONFIG.get('BIODATA', '/opt/biodata')

# Derived Constants
RPC_URL = os.path.join(CONFIG['target_server'], 'rpc')
GB_SLICE_URL = os.path.join(CONFIG['target_server'], 'api/genomebrowser')
DESIGNS = os.path.join(CONFIG['biodata_root'], 'dnadesigns')
CREDENTIALS = (CONFIG['email'], CONFIG['password'])


@click.group()
def cli():
    pass


def filter_extensions(filepath):
    if filepath.endswith(('gb', 'gbk', 'genbank')):
        return True
    else:
        False


@cli.command()
def set_default(key, value):
    CONFIG[key] = value


@cli.command()
def view_default(key):
    click.echo(CONFIG[key])


def load_design_file(design_file_name):
    """Parse and prepare the DnaDesign Object"""
    design_filepath = os.path.join(BIODATA, 'dnadesigns', design_file_name)
    with open(design_filepath) as dnadesign_file:
        data = parsers.genbank.parse(dnadesign_file)[0]
        data['sha1'] = data['sequence']['sha1']
        data['sequence'] = data['sequence']['seq']
        data.pop('file_contents')
        data['name'] = design_file_name
    return data


def annotate_design(dna_design):
    """Annotate the DNA Design"""
    pass


def validate_solution(solution_json):
    """Validate a solution object"""
    is_circular = solution_json['is_circular']
    if is_circular:
        assert solution_json['backbone']


@cli.command('fetch')
@click.argument('object')
@click.option('--filters', default=None)
@click.option('--output', type=click.File(), default=sys.stdout)
def fetch_cmd(object, filters, output):
    """Load all of DNA Molecule, Design, Feature, or Annotation from a
    remote source"""
    filters = filters if filters else {}
    if object in GB_MODELS:
        service = 'genomebrowser'
        body = gb.make_fetch_instruction(object, filters)
    else:
        service = 'inventory'
        body = {'object': object, 'filters': filters}

    endpoint_url = os.path.join(CONFIG['target_server'], 'api/{0}/crud'.format(service))
    credentials = (CONFIG['email'], CONFIG['password'])
    resp = requests.post(endpoint_url, json.dumps(body), auth=credentials)
    try:
        read_list = resp.json()['read']
        json.dump(read_list, output, indent=2)
    except KeyError:
        pprint.pprint(resp.text)
    except ValueError:
        pprint.pprint(resp.text)


@cli.command('slice_genome')
@click.argument('chromosome', help="Target Chromosome, ie. chr1")
@click.argument('start_end', nargs=2, type=click.IntRange(), help="Interval")
@click.argument('tracks', type=click.Choice(GB_MODELS))
@click.option('--genome', default=CONFIG['genome'], help="Target Genome")
@click.option('--strand', default=0)
@click.option('--output', type=click.File('wb'), default=sys.stdout)
@click.option('--sequence', default=False)
@click.option('--nuclease', default='wtCas9')
def slice_genome_cmd(genome, chromosome, start_end, tracks, strand, output,
                     sequence, nuclease):
    """"
    Load all the track data in the genome browser in an slice interval from
    a remote source.
    """
    target_url = GB_SLICE_URL
    # FUTURE - genomebrowser should just support dynamic tracks for guides,
    # talens, and other features. For now we use this work around.
    # FUTURE - allow other intersections besides overlaps
    guides = False
    if 'guides' in tracks:
        guides = True
        del tracks[tracks.index('guides')]

    instruction = gb.make_slice_instruction(genome, chromosome, start_end, strand, sequence)  # NOQA
    utils.make_post(target_url, instruction, output)
    # this is a work around as RPC does the guides while GB does the rest
    if guides:
        nuc_slice = ge.make_slice_instruction(chromosome, start_end, nuclease)
        utils.make_json_rpc(0, 'load_guides', nuc_slice, output)


@cli.command('score_guides')
@click.argument('guides', type=click.File(), default=sys.stdin, help='stream of mapped guides')
@click.option('--genome', default=CONFIG['genome'])
@click.option('--output', '-o', type=click.File(), default=sys.stdout)
@click.option('--async', default=False, type=click.BOOL)
@click.option('--activity', default=CONFIG['activity_params'])
@click.option('--offtarget', default=CONFIG['offtarget_params'])
def score_guides_cmd(guides, genome, async, activity, offtarget, output):
    """
    Analyze a batch of guide RNAs
    :return:
    """
    endpoint_url = RPC_URL
    credentials = CREDENTIALS
    method = 'score_guide'
    jobs = ge.make_score_guide_instruction(guides, genome, activity, offtarget)

    def on_error(response):
        click.echo(response.text)

    utils.make_batch_rpc(endpoint_url, credentials, jobs, method, output,
                         on_error, async=False)


@cli.command('run_design')
@click.argument('spec_file', type=click.File(), default=sys.stdin)
@click.option('--async', default=False, type=click.BOOL)
@click.option('--output', '-o', type=click.File(), default=sys.stdout)
def start_design(spec_file, async, output):
    """Start DNA Design Job to generate guides, pairs, oligos, vectors, donors,
    etc. Libraries are simply batches/chains of the above.

    The specfile defines a canvas of more atomic methods which are composed
    together and run either synchronously or asynchronously.

    """
    utils.make_json_rpc()

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

if __name__ == '__main__':
    cli()







