#! /usr/bin/env python


from setuptools import setup, find_packages

setup(name='dgcli',
      version='0.0.1',
      url='https://api.deskgen.com',
      author='Desktop Genetics Ltd',
      maintainer='Desktop Genetics Ltd',
      maintainer_email='devs@desktopgenetics.com',
      description="Command Line Interface for the DeskGen Platform",
      packages=find_packages(),
      entry_points={
          'console_scripts': [
              'dg = cli_interface:cli',
          ]
      },
      install_requires=[
          'click',
          'requests',
          'pyyaml',
          'xlrd',
      ],
      tests_require=[
          'tox',
          'pytest',
          'mock',
      ],
      data_files=[
          'config/file_convention.yaml',
          'config/remote_schema.yaml',
          'config/sample.dgrc'
      ]
)
