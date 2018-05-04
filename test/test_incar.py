import os
import pytest
import numpy as np
from parsevasp.incar import Incar

@pytest.fixture(scope = 'module', params=[0])
def incar_parser(request, tmpdir_factory):
    """Load INCAR file.

    """
    testdir = os.path.dirname(__file__)
    incarfile = testdir + "/INCAR"
    tmpfile = str(tmpdir_factory.mktemp('data').join('INCAR'))
    incar_truncate(request.param, incarfile, tmpfile)
    incar = Incar(file_path = tmpfile)
    
    return incar

def incar_truncate(index, original, tmp):
    """Truncate the INCAR file.

    """
    
    with open(original, 'r') as incarfile:
        content = incarfile.read().splitlines()
    truncated_content = '\n'.join(content[:-index or None])
    with open(tmp, 'w') as incarfile:
        incarfile.write(str(truncated_content))

    return

def test_incar_exist(incar_parser):
    """Check if incar_parser exists.

    """

    assert incar_parser.get_dict()
