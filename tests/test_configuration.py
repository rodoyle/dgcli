"""
Tests of configuration (.dgrc) loading and validation
"""
import os

from dgcli import configuration


def test_load_config():
    test_file_path = os.path.join(os.path.dirname(__file__), '../config/sample.dgrc')
    config = configuration.load(test_file_path)
    assert config
