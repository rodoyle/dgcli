#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for working with the Inventory service
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import logging
import os
import re

import requests
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


class AuthenticationError(Exception):
    pass


class InventoryService(object):
    """
    Represents the Inventory Service and its public methods.
    """

    def __init__(self, target_server):
        self.target_server = target_server
        self.crud_url = os.path.join(self.target_server, 'api/inventory/crud')
        self.upload_url = os.path.join(self.target_server, 'api/inventory/upload')
        self.login_url = os.path.join(self.target_server, 'api/inventory/authenticate')
        self.session = requests.Session()

    def set_credentials(self, username, password):
        self.session.auth = (username, password)

    def _get_schema(self, record):
        type_ = record.get('type_')
        return dgparse.VALIDATORS[type_]

    def login(self, username, password):
        self.set_credentials(username, password)
        resp = self.session.post(self.login_url, {})
        # Copy Session Token
        if resp.ok:
            self.session.cookies['session'] = resp.cookies['session']
        else:
            raise AuthenticationError

    def create(self, record):
        """Creates a record"""
        schema = self._get_schema(record)
        body, errors = schema.dumps(record)
        if errors:
            return record, errors
        resp = self.session.post(self.crud_url, json=body)
        if resp.ok:
            return resp.json()['read'], {}
        else:
            try:
                return record, resp.json()
            except ValueError:
                return record, {'create': resp.text}


def construct_service(config):
    inv_serve = InventoryService(config.target_server)
    return inv_serve



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

