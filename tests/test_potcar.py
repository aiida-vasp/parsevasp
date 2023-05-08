"""Test potcar."""
import os

import numpy as np
import pytest

from parsevasp.potcar import Potcar


@pytest.fixture(scope='module')
def potcar_reference_metadata():
    """Reference metadata for a dummy POTCAR"""
    metadata = {
        'VRHFIN': 'In: s1p1',
        'LEXCH': 'CA',
        'EATOM': 1111.1111,
        'TITEL': 'PAW In_d 11Feb1111',
        'LULTRA': False,
        'IUNSCR': 1,
        'RPACOR': 1.111,
        'POMASS': 111.111,
        'ZVAL': 11.111,
        'RCORE': 1.111,
        'RWIGS': 1.111,
        'ENMAX': 111.111,
        'ENMIN': 111.111,
        'ICORE': 1,
        'LCOR': True,
        'LPAW': True,
        'EAUG': 111.111,
        'DEXC': 1.111,
        'RMAX': 1.111,
        'RAUG': 1.111,
        'RDEP': 1.111,
        'RDEPT': 1.111,
    }

    return metadata


@pytest.fixture(scope='module')
def potcar_parser(request):
    """Load POTCAR file."""

    testdir = os.path.dirname(__file__)
    potcarfile = os.path.join(testdir, 'POTCAR')
    potcar = Potcar(file_path=potcarfile)

    return potcar


@pytest.fixture(scope='module')
def potcar_parser_file_object():
    """Load POTCAR file using a file object."""

    testdir = os.path.dirname(__file__)
    potcarfile = os.path.join(testdir, 'POTCAR')
    potcar = None
    with open(potcarfile, 'r', encoding='utf8') as file_handler:
        potcar = Potcar(file_handler=file_handler)

    return potcar


@pytest.mark.parametrize(
    'potcar_object,reference_values', [
        ('potcar_parser', 'potcar_reference_metadata'),
        ('potcar_parser_file_object', 'potcar_reference_metadata'),
    ]
)
def test_potcar_metadata(potcar_object, reference_values, request):
    """Test if the metadata produced matches the expected one"""
    potcar_object = request.getfixturevalue(potcar_object)
    reference_values = request.getfixturevalue(reference_values)

    for key, value in reference_values.items():
        assert key in potcar_object.metadata, f'key "{key}" not in the metadata'
        assert value == potcar_object.metadata[
            key
        ], f'referance value "{value}" does not match to found value {potcar_object.metadata[key]} for key "{key}"'


@pytest.mark.parametrize(
    'potcar_object,reference_values', [
        ('potcar_parser', 'potcar_reference_metadata'),
        ('potcar_parser_file_object', 'potcar_reference_metadata'),
    ]
)
def test_potcar_attributes(potcar_object, reference_values, request):
    """Test if the attributes are correctly setup"""
    potcar_object = request.getfixturevalue(potcar_object)
    reference_values = request.getfixturevalue(reference_values)

    extra_values = {'symbol': 'In_d', 'element': 'In', 'functional': 'Perdew-Zunger81', 'functional_class': 'LDA'}

    reference_values.update(extra_values)

    for key, value in reference_values.items():
        assert hasattr(potcar_object, key.lower()), f'attribute {key} not found in potcar'
        assert getattr(
            potcar_object, key.lower()
        ) == value, f'value of attribute {key} {getattr(potcar_object, key.lower())} does not match reference {value}'
