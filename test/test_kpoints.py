import os
import pytest
import numpy as np
from parsevasp.kpoints import Kpoints

@pytest.fixture(scope = 'module', params=[0])
def kpoints_parser(request, tmpdir_factory):
    """Load KPOINTS file.

    """
    testdir = os.path.dirname(__file__)
    kpointsfile = testdir + "/KPOINTS"
    tmpfile = str(tmpdir_factory.mktemp('data').join('KPOINTS'))
    kpoints_truncate(request.param, kpointsfile, tmpfile)
    kpoints = Kpoints(file_path = tmpfile)
    
    return kpoints

def kpoints_truncate(index, original, tmp):
    """Truncate the KPOINTS file.

    """
    
    with open(original, 'r') as kpointsfile:
        content = kpointsfile.read().splitlines()
    truncated_content = '\n'.join(content[:-index or None])
    with open(tmp, 'w') as kpointsfile:
        kpointsfile.write(str(truncated_content))

    return

def test_kpoints_exist(kpoints_parser):
    """Check if kpoints_parser exists.

    """

    assert kpoints_parser.get_dict()
