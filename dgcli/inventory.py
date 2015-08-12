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

from marshmallow import Schema, fields

import requests
import dgparse

MODELS = {
    'oligo': 'oligo',
    'primer': 'oligo',
    'plasmid': 'dnamolecule',
    'dnafeature': 'dnafeature',
    'dnamolecule': 'dnamolecule',
    'dnadesign': 'dnamolecule',
    'construct': 'dnamolecule'
}

log = logging.getLogger(__name__)


class CreateInstructionSchema(Schema):
    object = fields.String(required=True, attribute="object")
    create = fields.List(fields.Raw, required=True)


class AuthenticationError(Exception):
    pass


class InventoryClient(object):
    """
    Represents the Inventory Service and its public methods.
    """
    crudreq_schema = CreateInstructionSchema()

    def __init__(self, target_server):
        self.target_server = target_server
        self.crud_url = os.path.join(self.target_server, 'api/inventory/crud')
        self.upload_url = os.path.join(self.target_server, 'api/inventory/upload')
        self.login_url = os.path.join(self.target_server, 'api/inventory/authenticate')
        self.session = requests.Session()

    def set_credentials(self, username, password):
        self.session.auth = (username, password)

    def _get_schema(self, record_type):
        return dgparse.VALIDATORS[record_type]

    def login(self, username, password):
        self.set_credentials(username, password)
        resp = self.session.post(self.login_url, {})
        # Copy Session Token
        if resp.ok:
            self.session.cookies['session'] = resp.cookies['session']
        else:
            raise AuthenticationError

    def _handle_response(self, record, resp, key='create'):
        if resp.ok:
            return resp.json()[key], {}
        else:
            try:
                return record, resp.json()
            except ValueError:
                return record, {key: resp.text}

    def create(self, object_, record):
        """Creates a record"""
        schema = self._get_schema(object_)
        data, errors = schema.dump(record)
        errors.update(schema.validate(record))
        if errors:
            return record, errors
        instruction = {
            'object': MODELS[object_],
            'create': data if isinstance(data, list) else [data],
        }
        # validate instruction
        body, errors = self.crudreq_schema.dump(instruction)
        resp = self.session.post(self.crud_url, json=body)
        return self._handle_response(record, resp, key='create')

    def upload(self, openfile):
        """Upload a file"""
        instruction = {"inventoryfile": openfile}
        resp = self.session.post(self.upload_url, files=instruction)
        return self._handle_response(openfile, resp, key='id_')


def construct_service(config):
    inv_serve = InventoryClient(config.target_server)
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

