"""System level test of the guide services
"""

from __future__ import unicode_literals

import csv
import json
import logging
import os
import pprint
import sys
import uuid
import yaml
import xlrd

import click
import requests

import defaults

logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

BIODATA = os.environ.get('BIODATA', '/opt/biodata')
TARGET_SERVER = 'https://api.deskgen.com'  #default

CONFIG = {
    'BIODATA': BIODATA,
    'email': None,
    'password': None
}


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


@cli.command()
@click.argument('design_path', type=click.Path())
@click.argument('target_server')
@click.argument('email')
@click.argument('password')
def autoclone(design_path, target_server, email, password):
    """
    AutoClone a DNA Design
    """
    url = os.path.join(target_server, 'rpc')
    with open(design_path, 'r') as design_file:
        # first try opening a valid json file from a previous command
        dna_design = design_file.read()
        # next try parsing the file


    body = {
        'jsonrpc': '2.0',
        'method': 'find_cloning_solutions',
        'id': 1,
        'params': {
            'async': False,
            'dnadesign': dna_design,
        }
    }
    resp = requests.post(url,
                         json.dumps(body),
                         auth=(email, password))  # required

    if resp.ok:
        msg = "{0} cloned".format(dna_design['name'])
        log.info(msg)
        click.echo(pprint.pprint(resp.json()['result']))
    else:
        msg = "Cloning Solution Generation failed for {0}".format(dna_design['name'])  # NOQA
        log.warn(msg)
        log.warn(pprint.pprint(resp.text))


@cli.command()
@click.argument('target_server')
@click.argument('object')
@click.option('--name', default=None)
def read(target_server, object, name=None):
    """Load all of an object"""
    if object in ('genome', 'chromosome', 'gene', 'transcript',
                  'exon', 'cds', 'cdsregion', 'trackedguide'):
        service = 'genomebrowser'
    else:
        service = 'inventory'
    endpoint_url = os.path.join(target_server, 'api/{0}/crud'.format(service))
    filters = {}
    if name:
        filters.update({'name': name})
    body = {
        'object': object,
        'read': {
            'filter': filters,
        }
    }
    credintials = ('edwardp@deskgen.com', 'VALIS')
    resp = requests.post(endpoint_url, json.dumps(body), auth=credintials)
    try:
        read_list = resp.json()['read']
        click.echo(pprint.pprint(read_list))
        return read_list
    except KeyError:
        pprint.pprint(resp.text)
    except ValueError:
        pprint.pprint(resp.text)

@cli.command()
@click.argument('target_server')
@click.argument('gene_name')
@click.argument('output', type=click.File('wb'))
def load_gene(target_server, gene_name, output):
    """Load a gene and all the target region data from the genome browser"""
    gene_instruction = {
        'object': 'gene',
        'read': {
            'filter': {'name': gene_name},
            'expand': [['coding_sequences', 'regions', 'coordinates'],
                       ['transcripts', 'exons', 'coordinates'],
                       ['coordinates']
                       ],
        }
    }
    url = os.path.join(target_server, 'api/genomebrowser/crud')
    resp = requests.post(
        url,
        json.dumps(gene_instruction),
        headers={b'content-type': b'application/json'}
    )
    try:
        gene_data = resp.json()['read'][0]
        yaml.dump(gene_data, output)
    except KeyError:
        msg = "Improper Query to GenomeBrowser"
        log.error(msg)
        log.error(resp.text)
        raise
    except IndexError:
        msg = "No Genes Found for filter {0}".format(gene_name)
        log.warn(msg)
        log.error(resp.text)
        raise


@cli.command()
@click.argument('target_server')
@click.argument('chr_name')
@click.argument('start_end', nargs=2, type=click.IntRange())
@click.option('--nuclease', default='wtCas9')
def load_guides(target_server, chr_name, start_end, nuclease):
    """
    Test the load guide services. This is a synchronous method.
    This really need to be improved to handle multiple genomes better.
    :return:
    """
    endpoint_url = 'rpc'
    target_url = os.path.join(target_server, endpoint_url)
    params = {
        "chromosome": {
            'name': chr_name,
        },
        "start_end": (start_end),
        "nuclease": {'name': nuclease},
        'async': False,
    }

    body = {
        'jsonrpc': '2.0',
        'method': 'load_guides',
        'params': params,
        'id': 1,
    }

    resp = requests.post(target_url,
                         json.dumps(body),
                         headers={b'content-type': b'application/json'})
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
@click.argument('target_server')
@click.argument('guide')
@click.argument('genome')
@click.option('--async', default=False, type=click.BOOL)
def score_guide(target_server, guide, genome, async):
    """
    Test the score_guide method
    :return:
    """
    endpoint_url = 'rpc'
    target_url = os.path.join(target_server, endpoint_url)
    body = {
        'jsonrpc': '2.0',
        'method': 'score_guide',
        'params': {
            "genome": {'version': genome},
            "guide": {
                'type_': 'GuideRna',
                'bases': guide,
                'pam': 'AGG',
                'coordinates': {u'chromosome_id': 1,
                   u'start_end': [32890642, 32890665],
                   u'strand': 1},
                'hits': None,
                'score': None,
                'mismatch': None,
                'gc_content': None,

            },
            "scoring_function": {
                'doench': {
                    'minscore': 0.0  # Make higher to prevent slow scoring
                },
                'mitv1_db': {
                    "maxmismatch": 3,
                    'offset': 0,
                }
            },
        "async": async,
        },
        'id': 1,
    }
    resp = requests.post(target_url,
                         json.dumps(body),
                         headers={b'content-type': b'application/json'},
                         )
    if async:
        task_id = resp.json()['result']['task_id']
        click.echo(task_id)
    else:
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
@click.argument('target_server', default='https://api-staging.deskgen.com')
@click.argument('gene')
@click.argument('genome')
@click.option('--async', default=False, type=click.BOOL)
def walk_gene(target_server, gene, genome, async):
    """Test you can precompute guides for a gene"""
    endpoint_url = 'rpc'
    target_url = os.path.join(target_server, endpoint_url)
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


@cli.command()
@click.argument('target_server', default='https://api-staging.deskgen.com')
@click.argument('gene')
@click.argument('genome', default='GRCh37.p13')
@click.argument('nuclease', default='Cas9D10A')
@click.option('--async', default=False, type=click.BOOL)
@click.option('--output', default=sys.stdout, type=click.File())
def run_pair_tornado(target_server, gene, genome, nuclease, async, output):
    nuclease_dict = {'name': nuclease}
    endpoint_url = 'rpc'
    target_url = os.path.join(target_server, endpoint_url)
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

@cli.command()
@click.argument('spec_file', default=sys.stdin, type=click.File())
@click.option('--async', default=False, type=click.BOOL)
@click.option('--dryrun', default=False, type=click.BOOL)
@click.option('--output', default=sys.stdout, type=click.File())
def design_library(spec_file, async, dryrun, output):
    """
    Make a request to the server to start a library design job
    :param spec_file:
    :param async:
    :param dryrun:
    :param output:
    :return:
    """
    endpoint_url = 'rpc'
    target_url = os.path.join(TARGET_SERVER, endpoint_url)
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
            "callbacks": specs.get('callbacks', None),
            "name": specs.get('name', "Custom_Library"),
            "description": specs.get('description'),
            "dry_run": dryrun,
            "async": async,
        },
    }
    resp = requests.post(target_url,
                         json.dumps(request),
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

@cli.command()
@click.argument('target_server')
@click.argument('task_id')
def check_status(target_server, task_id):

    endpoint_url = 'rpc'
    target_url = os.path.join(target_server, endpoint_url)
    body = {
        'jsonrpc': '2.0',
        'method': 'check_status',
        'id': task_id,
    }
    resp = requests.post(target_url,
                         json.dumps(body),
                         headers={b'content-type': b'application/json'},
                         )
    click.echo(resp.text)


def run_polling(target_server, scoring_queue, email, password):

    # Check our tasks
    endpoint_url = 'rpc'
    target_url = os.path.join(target_server, endpoint_url)
    failed_tasks = list()
    passed_tasks = list()

    while True:
        if len(scoring_queue) < 0:
            break
        task_id = scoring_queue.popleft()
        body = {
        'jsonrpc': '2.0',
        'method': 'check_status',
        'id': task_id,
        }
        resp = requests.post(target_url,
                         json.dumps(body),
                         headers={b'content-type': b'application/json'},
                         auth=(email, password))
        if resp.json()['result'] == 'PENDING':
            scoring_queue.append(task_id)
        log.info(resp.json()['result'])


@cli.command()
@click.argument('target_server')
@click.argument('genome')
def load_genome(target_server, genome):
    url = os.path.join(target_server, 'api/genomebrowser/crud')
    body = {
        'object': 'genome',
        'read': {
            'filter': {'version': genome},
            'expand': [['chromosomes']]
        }
    }
    resp = requests.post(
        url,
        json.dumps(body),
        headers={b'content-type': b'application/json'}
    )
    click.echo(resp.text)


def parse_excel_file(filehandle):
    """
    Parse an excel file and return a list of targets
    :param filehandle:
    :return:
    """
    targets = []

    target_file = xlrd.open_workbook(filehandle)
    click.echo("Opended XLS workbook of targets")
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
        specs['targets'] = parse_excel_file(target_file.name)
    if update:
        with open(spec_file_path, b'w') as spec_file:
            json.dump(specs, spec_file, indent=2)
    else:
        json.dump(specs, output, indent=2)

if __name__ == '__main__':
    cli()







