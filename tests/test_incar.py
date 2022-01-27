"""Test incar."""
import os

import numpy as np
import pytest

from parsevasp.incar import Incar, IncarItem


@pytest.fixture()
def incar_dict():
    """Create a dictionary of valid INCAR items.

    """

    incar_dict = {'encut': 350, 'Sigma': '.5e-1 #comment', 'lreal': False, 'PREC': 'Accurate'}
    return incar_dict


@pytest.fixture(scope='module', params=[0])
def incar_parser(request, tmpdir_factory):
    """Load INCAR file.

    """
    testdir = os.path.dirname(__file__)
    incarfile = testdir + '/INCAR'
    tmpfile = str(tmpdir_factory.mktemp('data').join('INCAR'))
    incar_truncate(request.param, incarfile, tmpfile)
    incar = Incar(file_path=tmpfile)

    return incar


@pytest.fixture(scope='module', params=[0])
def incar_parser_file_object(request, tmpdir_factory):
    """Load INCAR file using a file object.

    """
    testdir = os.path.dirname(__file__)
    incarfile = testdir + '/INCAR'
    tmpfile = str(tmpdir_factory.mktemp('data').join('INCAR'))
    incar_truncate(request.param, incarfile, tmpfile)
    incar = None
    with open(tmpfile) as file_handler:
        incar = Incar(file_handler=file_handler)

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


def test_incar_parser_read(incar_parser):
    """Check if incar_parser exists.

    """

    assert incar_parser.get_dict()


def test_incar_parser_parameters(incar_parser):
    """Check parameters of the INCAR.

    """

    dictionary = incar_parser.get_dict()
    assert dictionary['emin'] == 5.5
    assert dictionary['emax'] == 7.5
    assert dictionary['nedos'] == 100000
    assert dictionary['prec'] == 'A'
    assert dictionary['loptics'] == True
    assert dictionary['encut'] == 350
    assert dictionary['dipol'] == [1, 2, 2]
    assert dictionary['ismear'] == -5
    assert dictionary['algo'] == 'V'


def test_incar_parser_write(incar_parser, tmp_path):
    """Check the write functions for both file paths and file objects.

    """
    incar = incar_parser.get_dict()
    # Write the content
    incar_write_path = tmp_path / 'INCAR'
    incar_parser.write(file_path=incar_write_path)
    # Then reload and compare
    incar_reloaded = Incar(file_path=incar_write_path).get_dict()
    assert incar == incar_reloaded
    # Write again with file object
    with open(incar_write_path, 'w') as handler:
        incar_parser.write(file_handler=handler)
    # Then reload again and compare
    with open(incar_write_path, 'r') as handler:
        incar_reloaded = Incar(file_handler=handler).get_dict()
    assert incar == incar_reloaded


def test_incar_parser_parameters_file_object(incar_parser_file_object):
    """Check parameters of the INCAR using a file object

    """
    dictionary = incar_parser_file_object.get_dict()
    assert dictionary['emin'] == 5.5
    assert dictionary['emax'] == 7.5
    assert dictionary['nedos'] == 100000
    assert dictionary['prec'] == 'A'
    assert dictionary['loptics'] == True
    assert dictionary['encut'] == 350
    assert dictionary['dipol'] == [1, 2, 2]
    assert dictionary['ismear'] == -5
    assert dictionary['algo'] == 'V'


def test_incar_from_dict(incar_dict):
    """Test passing a dictionary.

    """
    incar_io = Incar(incar_dict=incar_dict)
    comp_dict = {'encut': 350, 'sigma': 0.05, 'lreal': False, 'prec': 'Accurate'}
    assert str(sorted(incar_io.get_dict())) == str(sorted(comp_dict))


def test_incar_parser_from_string():
    """Test passing a string.

    """

    test_str = 'LOPTICS = .True.\nAddgrid=.false.'
    incar_io = Incar(incar_string=test_str)
    incar_dict = incar_io.get_dict()
    assert incar_dict.pop('loptics') is True
    assert incar_dict.pop('addgrid') is False
    assert not incar_dict


def test_incar_parser_from_string_complexr():
    """Test passing a more complex string.

    """

    test_string = '''LOPTICS = .True.
    EVENONLY = .False. # this is a comment; FLOAT\t=\t1.45e-03
    ISMEAR = THIS ; SIGMA = THAT
    NBANDS = 45  # endline comment; may contain '#' and ';' NOPARAM = this is not a parameter
    DIPOL = 1 2 -33 5
    '''
    parsed = Incar(incar_string=test_string)
    incar_dict = parsed.get_dict()
    assert incar_dict['loptics'] is True
    assert incar_dict['evenonly'] is False
    assert incar_dict['ismear'] == 'THIS'
    assert incar_dict['sigma'] == 'THAT'
    assert incar_dict['dipol'] == [1, 2, -33, 5]
    assert incar_dict['nbands'] == 45
    assert 'noparam' not in incar_dict
    assert 'float' not in incar_dict


def test_incar_parser_invalid_tag():
    """Test passing a tag that is not recognized.

    """

    test_string = '''SOMEINVALIDTAG = .TRUE.'''
    with pytest.raises(SystemExit):
        parsed = Incar(incar_string=test_string)


def test_incar_parser_invalid_tag():
    """Test passing a tag that is not recognized and its override.

    """

    test_string = '''SOMEINVALIDTAG = .TRUE.'''
    parsed = Incar(incar_string=test_string, validate_tags=False)
    incar_dict = parsed.get_dict()
    assert list(incar_dict.keys())[0] == 'someinvalidtag'


def test_incar_item():
    """Test the incar item class.

    """

    item = IncarItem(tag='encut', value=350, comment='    test comment ')
    assert item.get_tag() == 'encut'
    assert item.get_value() == 350
    assert item.get_comment() == 'test comment'
