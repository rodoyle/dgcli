#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import pytest
from dgcli import inventory as inv


@pytest.mark.parametrize("filename", [
    "EB0009 Leu2 integration vector 1x TetO CFP HisG without tTA.gb"
])
def test_get_file_conventions(filename):
    conf_path = os.path.normpath(
        os.path.join(__file__, "../../config/file_convention.yaml")
    )
    hits = 0
    for validator in inv.get_file_conventions(conf_path):
        result = validator(filename)
        if result:
            hits += 1
    assert hits == 1