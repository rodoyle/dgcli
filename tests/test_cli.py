"""
Tests of the command line interface itself.
"""
import os
import cli_interface
from click.testing import CliRunner
import pytest

from dgcli.configuration import load


@pytest.fixture
def cli_runner():
    return CliRunner()


def test_create_cmd(cli_runner):
    file_path = os.path.join(os.path.dirname(__file__), '../data/plasmid.csv')
    ctx = {'obj': load(os.path.join(os.path.dirname(__file__), '../config/sample.dgrc'))}
    setattr(cli_interface.create_cmd, 'context_settings', ctx)
    result = cli_runner.invoke(cli_interface.create_cmd, [file_path])
    assert result.exit_code == 0
    assert result.output