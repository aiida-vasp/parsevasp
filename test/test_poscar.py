import os
import pytest
import numpy as np
from parsevasp.poscar import Poscar, Site

@pytest.fixture(scope = 'module')
def poscar_parser():
    """Load POSCAR file.

    """
    
    testdir = os.path.dirname(__file__)
    poscarfile = testdir + "/POSCAR"
    poscar = Poscar(file_path = poscarfile)
    
    return poscar

@pytest.fixture(scope = 'module')
def poscar_parser_vel():
    """Load POSCAR file.

    """
    
    testdir = os.path.dirname(__file__)
    poscarfile = testdir + "/POSCARVEL"
    poscar = Poscar(file_path = poscarfile)
    
    return poscar

def test_poscar_exist(poscar_parser):
    """Check if poscar_parser exists.

    """

    assert poscar_parser.get_dict()

def test_poscar_entries(poscar_parser):
    """Check POSCAR entries.

    """
    
    poscar = poscar_parser.get_dict()
    assert poscar['comment'] == 'Compound: Co7Sb24.'
    unitcell = poscar['unitcell']
    test = np.array([[ 9.0164589999999993,  0.      ,  0.      ],
                     [ 0.      ,  9.0164589999999993,  0.      ],
                     [ 0.      ,  0.      ,  9.0164589999999993]])
    np.testing.assert_allclose(unitcell, test)
    sites = poscar['sites']
    assert len(sites) == 32
    test = ['Co', np.array([ 0.24999947,  0.24999947,  0.24999947]),
            [True, True, True], None, None, True]
    np.testing.assert_allclose(sites[0][1], test[1])
    assert sites[0][0] == test[0]
    assert sites[0][2] == [True, True, True]
    assert sites[0][3] == None
    np.testing.assert_allclose(sites[0][4], np.array([0.0, 0.0, 0.0]))
    assert sites[0][5]
    test = ['Co', np.array([ 0.74999953,  0.24999947,  0.24999947]),
            [True, True, True], None, None, True]
    np.testing.assert_allclose(sites[7][1], test[1])
    assert sites[7][0] == test[0]
    assert sites[7][2] == [True, True, True]
    assert sites[7][3] == None
    np.testing.assert_allclose(sites[7][4], np.array([8.0, 0.0, 0.0]))
    assert sites[7][5]
    test = ['Sb', np.array([ 0.        ,  0.33510051,  0.15804985]),
            [True, True, True], None, None, True]
    np.testing.assert_allclose(sites[8][1], test[1])
    assert sites[8][0] == test[0]
    assert sites[8][2] == [True, True, True]
    assert sites[8][3] == None
    np.testing.assert_allclose(sites[8][4], np.array([0.0, 9.0, 0.0]))
    assert sites[8][5]

def test_poscar_entries(poscar_parser_vel):
    """Check POSCAR entries including velocities.

    """
    
    poscar = poscar_parser_vel.get_dict()
    #print(poscar)
    #sys.exit(1)
    assert poscar['comment'] == 'Compound: Co7Sb24.'
    unitcell = poscar['unitcell']
    test = np.array([[ 9.0164589999999993,  0.      ,  0.      ],
                     [ 0.      ,  9.0164589999999993,  0.      ],
                     [ 0.      ,  0.      ,  9.0164589999999993]])
    np.testing.assert_allclose(unitcell, test)
    sites = poscar['sites']
    assert len(sites) == 32
    test = ['Co', np.array([ 0.24999947,  0.24999947,  0.24999947]),
            [True, True, True], None, None, True]
    np.testing.assert_allclose(sites[0][1], test[1])
    assert sites[0][0] == test[0]
    assert sites[0][2] == [True, True, True]
    np.testing.assert_allclose(sites[0][3], np.array([1.0, 0.0, 0.0]))
    np.testing.assert_allclose(sites[0][4], np.array([0.0, 0.0, 0.0]))
    assert sites[0][5]
    test = ['Co', np.array([ 0.74999953,  0.24999947,  0.24999947]),
            [True, True, True], None, None, True]
    np.testing.assert_allclose(sites[7][1], test[1])
    assert sites[7][0] == test[0]
    assert sites[7][2] == [True, True, True]
    np.testing.assert_allclose(sites[7][3], np.array([10.0, 1.0, 0.0]))
    np.testing.assert_allclose(sites[7][4], np.array([8.0, 0.0, 0.0]))
    assert sites[7][5]
    test = ['Sb', np.array([ 0.        ,  0.33510051,  0.15804985]),
            [True, True, True], None, None, True]
    np.testing.assert_allclose(sites[8][1], test[1])
    assert sites[8][0] == test[0]
    assert sites[8][2] == [True, True, True]
    np.testing.assert_allclose(sites[8][3], np.array([0.0, 0.0, 14.0]))
    np.testing.assert_allclose(sites[8][4], np.array([0.0, 9.0, 0.0]))
    assert sites[8][5]
    

def test_poscar_entries_string():
    """Test to check inititialization using string.

    """

    unitcell = np.array([[ 9.0164589999999993,  0.      ,  0.      ],
                         [ 0.      ,  9.0164589999999993,  0.      ],
                         [ 0.      ,  0.      ,  9.0164589999999993]])
    sites = []
    sites.append(Site('Co', np.array([0.0, 0.0, 0.0]),
                        [True, True, True], None, None, True))
    poscar_dict = {'comment': 'Example file',
                   'unitcell': unitcell,
                   'sites': sites}
    poscar_parser = Poscar(poscar_dict = poscar_dict)
    poscar = poscar_parser.get_dict()
    assert poscar['comment'] == 'Example file'
    unitcell = poscar['unitcell']
    test = np.array([[ 9.0164589999999993,  0.      ,  0.      ],
                     [ 0.      ,  9.0164589999999993,  0.      ],
                     [ 0.      ,  0.      ,  9.0164589999999993]])
    np.testing.assert_allclose(unitcell, test)
    sites = poscar['sites']
    assert len(sites) == 1
    test = ['Co', np.array([ 0.0,  0.0,  0.0]),
            [True, True, True], None, None, True]
    np.testing.assert_allclose(sites[0][1], test[1])
    assert sites[0][0] == test[0]
    assert sites[0][2] == [True, True, True]
    assert sites[0][3] == None
    assert sites[0][4] == None
    assert sites[0][5]
