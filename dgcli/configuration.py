"""
Configuration and Context Object
Contains parameters available to all commands
"""

import os

from marshmallow import Schema, fields
import yaml


class ServerConfig(object):
    """
    Represents the Target server
    """
    def __init__(self, url):
        self.url = url


class UserConfig(object):
    def __init__(self, email, password, biodata_root):
        self.email = email
        self.password = password
        self.credentials = (self.email, self.password)
        self.biodata_root = biodata_root


class Config(object):

    def __init__(self, config_file_path, target_server, user):
        self.config_file_path = config_file_path
        self.target_server = target_server
        self.user = user


class ServerSchema(Schema):
    url = fields.Url(required=True)

    def make_object(self, data):
        return ServerConfig(**data)


class UserSchema(Schema):
    name = fields.String()
    email = fields.Email(required=True)
    password = fields.String(required=True)
    biodata_root = fields.String()

    def make_object(self, data):
        return UserConfig(**data)


class ConfigSchema(Schema):
    target_server = fields.Nested(ServerSchema)
    user = fields.Nested(UserSchema)


def load(config_file_path):
    config_schema = ConfigSchema()
    with open(config_file_path, 'r') as config_file:
        params = yaml.load(config_file)
        params.update(os.environ)
        data, errors = config_schema.load(params)
    data['config_file_path'] = config_file_path
    config = Config(**data)
    return config
