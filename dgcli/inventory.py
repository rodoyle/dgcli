#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for working with the Inventory service
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import json
import os

import requests

import dgparse

INVENTORY_CRUD = 'api/inventory/crud'

EXTENSION_MAPPING = {
    '.dna': 'dnafeature',
    # '.seq': 'dnamoleculefile',
}

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

PARSE_ORDER = ['organisation', 'user', 'dnafeaturepattern', 'dnafeaturecategory',
               'dnafeature', 'dnadesign', 'dnamoleculesequence', 'dnamolecule',
               'dnamoleculefile',
               'dnadesigndnafeature', 'dnamoleculednafeature']


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


def build_local(local_root):
    """Initialize the local registry"""
    local_registry = {}
    for root, dirs, files in os.walk(local_root):
        for name in files:
            extension = os.path.splitext(name)
            if extension in EXTENSION_MAPPING:
                model = EXTENSION_MAPPING[extension]
                # map unique constraints to file path
                key = tuple([model].extend(UNIQUE_CONSTRAINTS[model]))
                local_registry[key] = os.path.join(root, name)
    return local_registry


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
        for object in data:
            key = _generate_key(tablename, object)
            remote[key] = data['id']  #only keep the ID.
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
                        data[relation + '_id'] = local[relation_key] #replace with xref
                    else:
                        local[relation_key] = data[relation]
                local[key] = data

    # local should only contain new items now.


def sync(requestor, local_root):
    """
    Given a repistory root and a remote server, traverse the repository
    and compare the contents to the remote server to get a list of objects to
    be pushed, grouped by table.

    :param requestor: An object capable of making post requests.
    :param local_root:
    :param remote:
    :return:
    """
    local = build_local(local_root)
    remote = fetch_remote(requestor, PARSE_ORDER)
    find_new(local, remote)  #modify in place
    push_remote(requestor, local, remote)





