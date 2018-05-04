import os
import pytest
import numpy as np
from parsevasp.poscar import Poscar

@pytest.fixture(scope = 'module', params=[0])
def poscar_parser(request, tmpdir_factory):
    """Load POSCAR file.

    """
    testdir = os.path.dirname(__file__)
    poscarfile = testdir + "/POSCAR"
    tmpfile = str(tmpdir_factory.mktemp('data').join('POSCAR'))
    poscar_truncate(request.param, poscarfile, tmpfile)
    poscar = Poscar(file_path = tmpfile)
    
    return poscar

def poscar_truncate(index, original, tmp):
    """Truncate the POSCAR file.

    """
    
    with open(original, 'r') as poscarfile:
        content = poscarfile.read().splitlines()
    truncated_content = '\n'.join(content[:-index or None])
    with open(tmp, 'w') as poscarfile:
        poscarfile.write(str(truncated_content))

    return

def test_poscar_exist(poscar_parser):
    """Check if poscar_parser exists.

    """

    assert poscar_parser.get_dict()
