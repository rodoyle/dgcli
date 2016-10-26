#! /usr/bin/env python

from setuptools import setup, find_packages

setup(name='dgcli',
      version='0.0.2',
      url='https://www.deskgen.com',
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
          'dgparse',
          'click',
          'click_plugins',
          'requests',
          'pyyaml',
          'openpyxl',
          'xlsxwriter',
          'xlrd',
          'marshmallow',
          'responses'
      ],
      tests_require=[
          'pytest',
          'mock',
      ],
      data_files=[
          'config/sample.dgrc'
      ]
)
