#!/usr/bin/python
import sys
import io
import os
import logging
testdir = os.path.dirname(__file__)
srcdir = '../parsevasp'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))
import xml

logging.basicConfig()
logger = logging.getLogger('Testing')

vaspxml = xml.Xml(logger, testdir + "/vasprunpartial.xml")
