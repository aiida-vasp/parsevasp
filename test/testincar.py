#!/usr/bin/python
import sys
import io
import os
import logging
testdir = os.path.dirname(__file__)
srcdir = '../parsevasp'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))
import incar

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('Testing')

incar = incar.Incar(file_path = testdir + "/INCAR")

incar.modify("PREC", "N")
incar.modify("TEST", [2, 3, 4])
incar.modify("TEST", "4 5 6")

incar.write(file_path = testdir + "/INCARMOD", comments = True)
print incar.entries

print incar.get_dict()
