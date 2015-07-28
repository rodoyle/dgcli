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


class InventoryService(object):
    """
    Represents the Inventory Service and its public methods.
    """

    def __init__(self, target_server, credentials, on_success, on_error):
        self.target_server = target_server
        self.credentials = credentials
        self.on_success = on_success
        self.on_error = on_error
        self.crud_url = os.path.join(self.target_server, 'api/inventory/crud')
        self.upload_url = os.path.join(self.target_server, 'api/inventory/upload')

    def create(record):
        """Creates a record"""
        resp = requests.post(endpoint, dgparse.schema, auth=credentials)



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
    else:
        return {'accession': os.path.basename(path)}

