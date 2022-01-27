"""Test eigenval."""
import os

import numpy as np
import pytest

from parsevasp.eigenval import Eigenval

compare_eigenvalues = np.array([[[
    -1.439825, 2.964373, 2.964373, 2.964373, 7.254542, 7.254542, 7.254542, 11.451811, 11.670398, 11.670398
]]])
compare_kpoints = np.array([[0.25, 0.25, 0.25, 1.0]])
compare_metadata = {
    0: [4, 4, 1, 1],
    1: [16.48482, 4.04e-10, 4.04e-10, 4.04e-10, 1e-16],
    2: 0.0001,
    'n_ions': 4,
    'n_atoms': 4,
    'p00': 1,
    'nspin': 1,
    'cartesian': True,
    'name': 'unknown system',
    'some_num': 12,
    'n_bands': 10,
    'n_kp': 1
}


@pytest.fixture
def eigenval_parser(request):
    """Load EIGENVAL file."""
    try:
        name = request.param
    except AttributeError:
        # Test not parametrized
        name = 'EIGENVAL'
    testdir = os.path.dirname(__file__)
    eigenvalfile = testdir + '/' + name
    eigenval = Eigenval(file_path=eigenvalfile)

    return eigenval


@pytest.fixture
def eigenval_parser_file_object(request):
    """Load EIGENVAL file from a file object."""
    try:
        name = request.param
    except AttributeError:
        # Test not parametrized
        name = 'EIGENVAL'
    testdir = os.path.dirname(__file__)
    eigenvalfile = testdir + '/' + name
    eigenval = None
    with open(eigenvalfile) as file_handler:
        eigenval = Eigenval(file_handler=file_handler)

    return eigenval


def test_eigenval(eigenval_parser):
    """Test that the content returned by the EIGENVAL parser returns correct eigenvalues, kpoints and metadata."""
    eigenvalues = eigenval_parser.get_eigenvalues()
    kpoints = eigenval_parser.get_kpoints()
    metadata = eigenval_parser.get_metadata()
    assert np.allclose(eigenvalues, compare_eigenvalues)
    assert np.allclose(kpoints, compare_kpoints)
    assert metadata == compare_metadata
