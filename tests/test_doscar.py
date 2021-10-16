import os
import pytest
import numpy as np
from parsevasp.doscar import Doscar

compare_dos = np.array([(-3.44, -1.10400000e-43, -2.09900000e-43),
                        (-1.539, 1.40000000e-01, 2.66100000e-01),
                        (0.362, -3.62400000e-73, 2.00000000e+00),
                        (2.264, -1.33800000e-05, 2.00000000e+00),
                        (4.165, 3.15600000e+00, 8.00000000e+00),
                        (6.066, -2.41200000e-15, 8.00000000e+00),
                        (7.967, 3.15600000e+00, 1.40000000e+01),
                        (9.868, -1.38100000e-27, 1.40000000e+01),
                        (11.769, 2.90100000e+00, 1.95200000e+01),
                        (13.67, 0.00000000e+00, 2.00000000e+01)],
                       dtype=[('energy', '<f8'), ('total', '<f8'), ('integrated', '<f8')])  # yapf: disable

@pytest.fixture
def doscar_parser(request):
    """Load DOSCAR file.

    """
    try: 
        name = request.param
    except AttributeError:
        # Test not parametrized
        name = 'DOSCAR'
    testdir = os.path.dirname(__file__)
    doscarfile = testdir + '/' + name
    doscar = Doscar(file_path=doscarfile)

    return doscar

@pytest.fixture
def doscar_parser_file_object(request):
    """Load DOSCAR file from a file object.

    """
    try: 
        name = request.param
    except AttributeError:
        # Test not parametrized
        name = 'DOSCAR'
    testdir = os.path.dirname(__file__)
    doscarfile = testdir + '/' + name
    doscar = None
    with open(doscarfile) as file_handler:
        doscar = Doscar(file_handler=file_handler)

    return doscar

def test_doscar(doscar_parser):
    """Test that the content returned by the DOSCAR parser returns total density of states."""
    dos = doscar_parser.get_dos()
    for i in range(dos.size):
        assert dos[i] == compare_dos[i]

def test_doscar_file_objcet(doscar_parser_file_object):
    """
    Test that the content returned by the DOSCAR parser returns total density of states
    using file handler.
    """
    dos = doscar_parser_file_object.get_dos()
    for i in range(dos.size):
        assert dos[i] == compare_dos[i]


@pytest.mark.parametrize('doscar_parser', ['DOSCAR.nopdos'], indirect=['doscar_parser'])
def test_doscar(doscar_parser):
    """Test that the content returned by the DOSCAR parser returns total density of states."""
    dos = doscar_parser.get_dos()
    for i in range(dos.size):
        assert dos[i] == compare_dos[i]
        
def test_doscar_partial(doscar_parser):
    """Test that the content returned by the DOSCAR parser returns the partial density of states."""
    pdos = doscar_parser.get_pdos()
    metadata = doscar_parser.get_metadata()
    n_atoms = metadata['n_atoms']
    compare_pdos = np.array([[(-3.44 , -6.543e-45, -3.040e-46, -3.040e-46, -3.040e-46, 0., 0., 0., 0., 0.),
                              (-1.539,  8.295e-03,  3.853e-04,  3.853e-04,  3.853e-04, 0., 0., 0., 0., 0.),
                              ( 0.362, -1.336e-74, -3.400e-75, -3.400e-75, -3.400e-75, 0., 0., 0., 0., 0.),
                              ( 2.264, -4.932e-07, -1.255e-07, -1.255e-07, -1.255e-07, 0., 0., 0., 0., 0.),
                              ( 4.165,  1.163e-01,  2.961e-02,  2.961e-02,  2.961e-02, 0., 0., 0., 0., 0.),
                              ( 6.066, -5.870e-17, -2.973e-17, -2.973e-17, -2.973e-17, 0., 0., 0., 0., 0.),
                              ( 7.967,  7.681e-02,  3.890e-02,  3.890e-02,  3.890e-02, 0., 0., 0., 0., 0.),
                              ( 9.868, -1.988e-29, -1.936e-29, -1.923e-29, -1.944e-29, 0., 0., 0., 0., 0.),
                              (11.769,  5.486e-02,  3.117e-02,  3.938e-02,  3.827e-02, 0., 0., 0., 0., 0.),
                              (13.67 ,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 0., 0., 0., 0., 0.)],
                             [(-3.44 , -6.543e-45, -3.040e-46, -3.040e-46, -3.040e-46, 0., 0., 0., 0., 0.),
                              (-1.539,  8.295e-03,  3.853e-04,  3.853e-04,  3.853e-04, 0., 0., 0., 0., 0.),
                              ( 0.362, -1.336e-74, -3.400e-75, -3.400e-75, -3.400e-75, 0., 0., 0., 0., 0.),
                              ( 2.264, -4.932e-07, -1.255e-07, -1.255e-07, -1.255e-07, 0., 0., 0., 0., 0.),
                              ( 4.165,  1.163e-01,  2.961e-02,  2.961e-02,  2.961e-02, 0., 0., 0., 0., 0.),
                              ( 6.066, -5.870e-17, -2.973e-17, -2.973e-17, -2.973e-17, 0., 0., 0., 0., 0.),
                              ( 7.967,  7.681e-02,  3.890e-02,  3.890e-02,  3.890e-02, 0., 0., 0., 0., 0.),
                              ( 9.868, -1.988e-29, -1.936e-29, -1.923e-29, -1.944e-29, 0., 0., 0., 0., 0.),
                              (11.769,  5.486e-02,  3.117e-02,  3.938e-02,  3.827e-02, 0., 0., 0., 0., 0.),
                              (13.67 ,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 0., 0., 0., 0., 0.)],
                             [(-3.44 , -6.543e-45, -3.040e-46, -3.040e-46, -3.040e-46, 0., 0., 0., 0., 0.),
                              (-1.539,  8.295e-03,  3.853e-04,  3.853e-04,  3.853e-04, 0., 0., 0., 0., 0.),
                              ( 0.362, -1.336e-74, -3.400e-75, -3.400e-75, -3.400e-75, 0., 0., 0., 0., 0.),
                              ( 2.264, -4.932e-07, -1.255e-07, -1.255e-07, -1.255e-07, 0., 0., 0., 0., 0.),
                              ( 4.165,  1.163e-01,  2.961e-02,  2.961e-02,  2.961e-02, 0., 0., 0., 0., 0.),
                              ( 6.066, -5.870e-17, -2.973e-17, -2.973e-17, -2.973e-17, 0., 0., 0., 0., 0.),
                              ( 7.967,  7.681e-02,  3.890e-02,  3.890e-02,  3.890e-02, 0., 0., 0., 0., 0.),
                              ( 9.868, -1.988e-29, -1.936e-29, -1.923e-29, -1.944e-29, 0., 0., 0., 0., 0.),
                              (11.769,  5.486e-02,  3.117e-02,  3.938e-02,  3.827e-02, 0., 0., 0., 0., 0.),
                              (13.67 ,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 0., 0., 0., 0., 0.)],
                             [(-3.44 , -6.543e-45, -3.040e-46, -3.040e-46, -3.040e-46, 0., 0., 0., 0., 0.),
                              (-1.539,  8.295e-03,  3.853e-04,  3.853e-04,  3.853e-04, 0., 0., 0., 0., 0.),
                              ( 0.362, -1.336e-74, -3.400e-75, -3.400e-75, -3.400e-75, 0., 0., 0., 0., 0.),
                              ( 2.264, -4.932e-07, -1.255e-07, -1.255e-07, -1.255e-07, 0., 0., 0., 0., 0.),
                              ( 4.165,  1.163e-01,  2.961e-02,  2.961e-02,  2.961e-02, 0., 0., 0., 0., 0.),
                              ( 6.066, -5.870e-17, -2.973e-17, -2.973e-17, -2.973e-17, 0., 0., 0., 0., 0.),
                              ( 7.967,  7.681e-02,  3.890e-02,  3.890e-02,  3.890e-02, 0., 0., 0., 0., 0.),
                              ( 9.868, -1.988e-29, -1.936e-29, -1.923e-29, -1.944e-29, 0., 0., 0., 0., 0.),
                              (11.769,  5.486e-02,  3.117e-02,  3.938e-02,  3.827e-02, 0., 0., 0., 0., 0.),
                              (13.67 ,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00, 0., 0., 0., 0., 0.)]],
                            dtype=[('energy', float), ('s', float), ('py', float), ('px', float), ('pz', float),
                                   ('dxy', float), ('dyz', float), ('dz2', float), ('dxz', float),
                                   ('dx2-y2', float)])  # yapf: disable
    for i in range(n_atoms):
        for j in range(pdos[i].size):
            assert pdos[i,j] == compare_pdos[i,j]

@pytest.mark.parametrize('doscar_parser', ['DOSCAR.spin'], indirect=['doscar_parser'])
def test_doscar_spin(doscar_parser):
    """
    Test that the content returned by the DOSCAR parser returns 
    the correct dimensions for spin density of states
    """
    dos = doscar_parser.get_dos()
    assert len(dos.dtype) == 3
    assert dos['total'].shape == (301, 2)

    pdos = doscar_parser.get_pdos()
    assert len(pdos.dtype) == 10
    assert pdos[0]['px'].shape == (301, 2)

@pytest.mark.parametrize('doscar_parser', ['DOSCAR.ncl'], indirect=['doscar_parser'])
def test_doscar_spin(doscar_parser):
    """
    Test that the content returned by the DOSCAR parser returns 
    the correct dimensions for density of states for colinear calculations.
    """    
    dos = doscar_parser.get_dos()
    assert len(dos.dtype) == 3
    assert dos['total'].shape == (301,)

    pdos = doscar_parser.get_pdos()
    assert len(pdos.dtype) == 10
    assert pdos[0]['px'].shape == (301, 4)
