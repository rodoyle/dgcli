#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for working with the Inventory service
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import json

import requests

INVENTORY_CRUD = 'api/inventory/crud'

EXTENSION_MAPPING = {
    '.dna': 'dnamolecule',
    '.seq': 'dnamoleculefile'
}

UNIQUE_CONSTRAINTS = {
    'dnamoleculesequence': ['sha1'],
    'dnafeaturepattern': ['sha1'],
    'dnafeaturecategory': ['name']
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


def make_read_instruction(object_type, filters):
    """
    Return the JSON body for making a read request
    :param object_type:
    :param filters:
    :return:
    """
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

