#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for working with the Inventory service
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import json
import logging
import os
import re
import yaml
from collections import namedtuple

import requests
from click import echo
import dgparse

log = logging.getLogger(__name__)

INVENTORY_CRUD = 'api/inventory/crud'
UPLOAD_URL = 'api/inventory/upload'
DESIGN_UPLOAD_URL = 'uploadDesignFile'


UNIQUE_CONSTRAINTS = {
    'dnamoleculesequence': ['sha1'],
    'dnafeaturepattern': ['sha1'],
    'dnafeaturecategory': ['name'],
    'organisation': ['domain'],
    'user': ['email'],
    'dnamolecule': ['organisation_id', 'accession'],
    'dnadesign': ['user_id', 'sha1'],
    'dnamolecule_dnafeature': ['dnamolecule_id', 'dnafeature_id', 'start', 'end', 'strand'],
    'dnadesign_dnafeature': ['dnadesign_id', 'dnafeature_id', 'start', 'end', 'strand'],
}

RELATIONS = {
    'user': ['organisation'],
    'dnamolecule': ['user', 'organisation', 'sequence'],
    'dnamoleculefile': ['dnamolecule'],
    'dnadesign': ['user', 'organisation'],
    'dnafeature': ['user', 'organisation', 'pattern', 'category'],
    'dnamolecule_dnafeature': ['dnamolecule', 'dnafeature'],
    'dnadesign_dnafeature': ['dnadesign', 'dnafeature'],
}


def make_create_request(object_type, objects):
    """Utility to make a CRUD request to create several objects.
    The objects are always added for the user making the request.
    """
    # this function is sort of dumb but will almost certainly be extended
    # with parsing, validation, and error handling code.
    if not isinstance(objects, list):
        objects = [objects]
    return {
        'object': object_type,
        'create': objects,
    }


def make_read_instruction(object_type, filters=None):
    """
    Return the JSON body for making a read request
    :param object_type:
    :param filters:
    :return:
    """
    if not filters:
        filters = {}

    return {
        'object': object_type,
        'read': filters
    }


def get_or_create_object(endpoint, credentials, object_type, object_data):
    """Check if an object exists by querying for it. If it exists, return the
    ID. If it doesn't exist, create it, and return the ID of the new object"""
    # We need to know what filters are unique for this object
    filters = {}
    for constraint in UNIQUE_CONSTRAINTS.get(object_type):
        filters[constraint] = object_data[constraint]

    read_body = make_read_instruction(object_type, filters)

    try:
        resp = requests.post(endpoint, json.dumps(read_body), auth=credentials)
        return resp.json()['read']['id']
    except KeyError:
        create_inst = make_create_request(object_type, object_data)
        resp = requests.post(endpoint, json.dumps(create_inst),
                             auth=credentials)
        return resp.json()['read']['id']


def _generate_key(tablename, object):
    """Generate a deduplication key for the object"""
    key = [tablename]
    key.extend([object[constraint] for constraint in UNIQUE_CONSTRAINTS[tablename]])
    return tuple(key)


def _get_relation_tablename(relation):

    if relation == 'sequence':
        return 'dnamoleculesequence'
    if relation == 'pattern':
        return 'dnafeaturepattern'
    if relation == 'category':
        return 'dnafeaturecategory'
    return relation


def fetch_remote(requestor, models):
    """Fetch the data from the remote repository
    We only fetch the models that are locally defined.
    """
    remote = {}
    for query in map(make_read_instruction, models):
        tablename = query['object']
        resp = requestor(query)
        data = resp.json()['read']
        for obj in data:
            key = _generate_key(tablename, obj)
            remote[key] = data['id']  # only keep the ID.
    return remote


def _make_push(requestor, current_model, new_objects, remote):
    req = make_create_request(current_model, new_objects)
    resp = requestor(req)
    data = resp.json()['read']
    for item in data:
        item_key = _generate_key(current_model, item)
        remote[item_key] = item['id']  # store the newly mad ID


def push_remote(requestor, local, remote):
    """Bulk push the new objects and register the returned ids"""
    current_model = None
    new_objects = None
    to_push = sorted([(k, d) for k, d in local.iteritems])

    for key, data in to_push:

        if not current_model:
            current_model = key[0]
        if not new_objects:
            new_objects = [data]

        if key[0] != current_model:
            # we're switching models
            _make_push(requestor, current_model, new_objects, remote)
            # restart on the new slice
            current_model = key[0]
            new_objects = [data]
        else:
            new_objects.append(data)

    # clean up
    if new_objects:
        _make_push(requestor, current_model, new_objects, remote)


def find_new(local, remote):
    """Find the new elements"""
    for key in local.iterkeys():
        path = local.pop(key)  # remove the item already present
        tablename = key[0]
        if key not in remote:
            extension = os.path.splitext(local[key])
            parser = dgparse.formats.get_parser(extension)

            with open(path, 'rb') as file_handle:
                data = parser(file_handle)  # this needs some more logic
                # handle the children of this new item
                for relation in RELATIONS[tablename]:
                    relation_table = _get_relation_tablename(relation)
                    relation_key = _generate_key(relation_table, data[relation])

                    if relation_key in remote:  # the parent already exists
                        data.pop(relation)
                        data[relation + '_id'] = remote[relation_key]
                    elif (relation_key in local) and isinstance(local[relation_key], int):
                        data.pop(relation)  # remove the relation
                        data[relation + '_id'] = local[relation_key] # replace with xref
                    else:
                        local[relation_key] = data[relation]
                local[key] = data

    # local should only contain new items now.


def get_file_conventions(conf_path):
    """Return a function that will consider paths, check they conform to a
    schema, and return a list of the valid files"""

    with open(conf_path) as conf_handle:
        conf = yaml.load(conf_handle)

    validators = list()
    for identifier, params in conf.iteritems():
        name_tokens = params['file_pattern']['attributes']
        pattern = params['file_pattern']['regex']

        def validator(filename):
            match_obj = re.match(filename, pattern)
            if match_obj:
                data = {attr: match_obj.group(idx) for idx, attr in enumerate(name_tokens)}
                params.pop('file_pattern')
                data.update(params)
                return namedtuple('validator', **data)(**data)
            else:
                return None

        validators.append(validator)
    return validators


def upload_files(local_root, schema_path, file_conv_path):
    """"Loop over a directory and upload all valid files encountered"""
    with open(schema_path) as sp:
        schema = json.load(sp)

    validators = get_file_conventions(file_conv_path)
    for root, dirs, files in os.walk(local_root):
        for name in files:
            echo("inspecting {0}".format(name))
            for v in validators:
                abspath = os.path.join(root, name)
                valid_data = v(name)
                if valid_data:
                    url = schema['definitions'][valid_data.type]['upload_url']

                    yield (url, abspath)

def parse_file_name(path):
    """
    Parse file names in accession, name, and format
    :param path:
    :return dictionary of matched items:
    """
    matched = re.search("([A-Z]{2,}\d+)[\s\_]([^/]+)\.(\w+)$", path)
    if matched:
        return {
            'accession': matched.group(1),
            'name': matched.group(2),
            'format': matched.group(3),
        }
