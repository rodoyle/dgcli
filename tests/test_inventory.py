#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os

import mock
import pytest

from dgcli import inventory as inv
from dgcli.configuration import load


@pytest.fixture
def invsrv():
    config = load(os.path.join(os.path.dirname(__file__), '../config/sample.dgrc'))
    invsrv = inv.InventoryService(config.target_server.url)
    invsrv.set_credentials(config.user.email, config.user.password)
    return invsrv


def test_inventory_create_service(invsrv):
    record = {'type_': 'dnamolecule', 'accession': 'EB123'}
    expected_response = {'create': [record]}
    with mock.patch('dgcli.inventory.requests.post', mock.Mock(return_value=expected_response))
        data, errors = invsrv.create(record)
    assert data == record
    assert not errors

