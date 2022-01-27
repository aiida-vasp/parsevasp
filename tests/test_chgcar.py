"""Test chgcar."""
import os

import numpy as np
import pytest

from parsevasp.chgcar import Chgcar

compare_charge_density = np.array([[[0.09329446, 0.18658892, 0.27988338], [0.37317784, 0.4664723, 0.55976676],
                                    [0.65306122, 0.74635569, 0.83965015], [0.93294461, 1.02623907, 1.11953353]],
                                   [[1.21282799, 1.30612245, 1.39941691], [1.49271137, 1.58600583, 1.67930029],
                                    [1.77259475, 1.86588921, 1.95918367], [2.05247813, 2.14577259, 2.23906706]],
                                   [[2.33236152, 2.42565598, 2.51895044], [2.6122449, 2.70553936, 2.79883382],
                                    [2.89212828, 2.98542274, 3.0787172], [3.17201166, 3.26530612, 3.35860058]],
                                   [[3.45189504, 3.5451895, 3.63848397], [3.73177843, 3.82507289, 3.91836735],
                                    [4.01166181, 4.10495627, 4.19825073], [4.29154519, 4.38483965, 4.47813411]],
                                   [[4.57142857, 4.66472303, 4.75801749], [4.85131195, 4.94460641, 5.03790087],
                                    [5.13119534, 5.2244898, 5.31778426], [5.41107872, 5.50437318, 5.59766764]]])


@pytest.fixture
def chgcar_parser(request):
    """Load CHGCAR file."""
    try:
        name = request.param[0]
        handler = request.param[1]
    except AttributeError:
        # Test not parametrized
        name = 'CHGCAR'
        handler = True
    testdir = os.path.dirname(__file__)
    chgcarfile = testdir + '/' + name
    chgcar = None
    if handler:
        with open(chgcarfile) as file_handler:
            chgcar = Chgcar(file_handler=file_handler)
    else:
        chgcar = Chgcar(file_path=chgcarfile)

    return chgcar


@pytest.mark.parametrize('chgcar_parser', [('CHGCAR', True), ('CHGCAR', False)], indirect=True)
def test_charge_density(chgcar_parser):
    """Test that the content returned by the CHGCAR parser returns the correct charge density."""
    charge_density = chgcar_parser.charge_density
    assert np.allclose(charge_density, compare_charge_density)


@pytest.mark.parametrize('chgcar_parser', [('CHGCAR.spin', True)], indirect=True)
def test_magnetization_density(chgcar_parser):
    """Test that the content returned by the CHGCAR parser
    returns the correct charge and magnetization density.

    """
    charge_density = chgcar_parser.charge_density
    magnetization_density = chgcar_parser.magnetization_density
    assert np.allclose(charge_density, compare_charge_density)
    assert np.allclose(magnetization_density, compare_charge_density)


@pytest.mark.parametrize('chgcar_parser', [('CHGCAR.ncl', True)], indirect=True)
def test_magnetization_density_ncl(chgcar_parser):
    """Test that the content returned by the CHGCAR parser
    returns the correct charge and magnetization density for non-collinear calculations.

    """
    charge_density = chgcar_parser.charge_density
    magnetization_density = chgcar_parser.magnetization_density
    assert np.allclose(charge_density, compare_charge_density)
    assert set(['x', 'y', 'z']) == set(magnetization_density.keys())
    for key, item in magnetization_density.items():
        assert np.allclose(item, compare_charge_density)
