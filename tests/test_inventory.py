#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests of the Inventory Module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os

import mock
import pytest

from dgcli import inventory as inv
from dgcli.configuration import load
from collections import namedtuple

MockResponse = namedtuple('MockResponse', ('ok json text'))


@pytest.fixture
def invsrv():
    config = load(os.path.join(os.path.dirname(__file__), '../config/sample.dgrc'))
    invsrv = inv.InventoryClient(config.target_server.url)
    invsrv.set_credentials(config.user.email, config.user.password)
    return invsrv


@pytest.fixture
def construct_response():
    def response(expected_body):
        return MockResponse(True, lambda: expected_body, str(expected_body))
    return response


@pytest.fixture
def plasmid_record():
    return {
        'type_': 'dnamolecule',
        'accession': 'EB123',
        'name': "TestPlasmid",
        'sequence': {
            'bases': b"ACGT"*3,
        }
    }


def test_inventory_create_service(invsrv, construct_response, plasmid_record):
    expected_response = {'create': [plasmid_record]}

    # because this in a class you need to pass the (unused) self argument.
    def fake_endpoint(self, _, json):  # Requests wants this json kwarg
        assert json['create']
        return construct_response(expected_response)
    # Important: need to patch the Requests Session method
    with mock.patch('dgcli.inventory.requests.sessions.Session.post', fake_endpoint):
        data, errors = invsrv.create(plasmid_record)
    assert data == [plasmid_record]
    assert not errors

