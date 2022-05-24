"""
Install the parsevasp python package.

usage: pip install .
"""

import json
import os

from setuptools import find_packages, setup

SETUP_JSON_PATH = os.path.join(os.path.dirname(__file__), 'setup.json')
README_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'README.rst')

with open(README_PATH, 'r', encoding='utf8') as readme:
    LONG_DESCRIPTION = readme.read()

if __name__ == '__main__':
    with open(SETUP_JSON_PATH, 'r', encoding='utf8') as info:
        SETUP_KWARGS = json.load(info)
    setup(
        packages=find_packages(exclude=['test']),
        keywords=['VASP', 'parser', 'python', 'xml'],
        long_description=LONG_DESCRIPTION,
        **SETUP_KWARGS
    )
