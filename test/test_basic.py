import os
import pytest
import numpy
from parsevasp import xml

testdir = os.path.dirname(__file__)


def test_vasprun_parser():
    vaspxml = xml.Xml(testdir + "/vasprunbasic.xml")
    for status in ['initial', 'final', 'all']:
        assert vaspxml.get_energies(status)
