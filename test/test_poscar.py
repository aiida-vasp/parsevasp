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
    print sites
    assert len(sites) == 32
    test = {'specie': 'Co', 'position': np.array([ 0.24999947,  0.24999947,  0.24999947]),
            'selective': [True, True, True], 'velocities': None,
            'predictors': None, 'direct': True}
    np.testing.assert_allclose(sites[0]['position'], test['position'])
    assert sites[0]['specie'] == test['specie']
    assert sites[0]['selective'] == [True, True, True]
    assert sites[0]['velocities'] == None
    np.testing.assert_allclose(sites[0]['predictors'], np.array([0.0, 0.0, 0.0]))
    assert sites[0]['direct']
    test = {'specie': 'Co', 'position': np.array([ 0.74999953,  0.24999947,  0.24999947]),
            'selective': [True, True, True], 'velocties': None,
            'predictors': None, 'direct': True}
    np.testing.assert_allclose(sites[7]['position'], test['position'])
    assert sites[7]['specie'] == test['specie']
    assert sites[7]['selective'] == [True, True, True]
    assert sites[7]['velocities'] == None
    np.testing.assert_allclose(sites[7]['predictors'], np.array([8.0, 0.0, 0.0]))
    assert sites[7]['direct']
    test = {'specie': 'Sb', 'position': np.array([ 0.        ,  0.33510051,  0.15804985]),
            'selective': [True, True, True], 'velocities': None,
            'predictors': None, 'direct': True}
    np.testing.assert_allclose(sites[8]['position'], test['position'])
    assert sites[8]['specie'] == test['specie']
    assert sites[8]['selective'] == [True, True, True]
    assert sites[8]['velocities'] == None
    np.testing.assert_allclose(sites[8]['predictors'], np.array([0.0, 9.0, 0.0]))
    assert sites[8]['direct']


def test_poscar_entries_vel(poscar_parser_vel):
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
    test = {'specie': 'Co', 'position': np.array([ 0.24999947,  0.24999947,  0.24999947]),
            'selective': [True, True, True], 'velocities': None,
            'predictors': None, 'direct': True}
    np.testing.assert_allclose(sites[0]['position'], test['position'])
    assert sites[0]['specie'] == test['specie']
    assert sites[0]['selective'] == [True, True, True]
    np.testing.assert_allclose(sites[0]['velocities'], np.array([1.0, 0.0, 0.0]))
    np.testing.assert_allclose(sites[0]['predictors'], np.array([0.0, 0.0, 0.0]))
    assert sites[0]['direct']
    test = {'specie': 'Co', 'position': np.array([ 0.74999953,  0.24999947,  0.24999947]),
            'selective': [True, True, True], 'velocities': None,
            'predictors': None, 'direct': True}
    np.testing.assert_allclose(sites[7]['position'], test['position'])
    assert sites[7]['specie'] == test['specie']
    assert sites[7]['selective'] == [True, True, True]
    np.testing.assert_allclose(sites[7]['velocities'], np.array([10.0, 1.0, 0.0]))
    np.testing.assert_allclose(sites[7]['predictors'], np.array([8.0, 0.0, 0.0]))
    assert sites[7]['direct']
    test = {'specie': 'Sb', 'position': np.array([ 0.        ,  0.33510051,  0.15804985]),
            'selective': [True, True, True], 'velocties': None,
            'predictors': None, 'direct': True}
    np.testing.assert_allclose(sites[8]['position'], test['position'])
    assert sites[8]['specie'] == test['specie']
    assert sites[8]['selective'] == [True, True, True]
    np.testing.assert_allclose(sites[8]['velocities'], np.array([0.0, 0.0, 14.0]))
    np.testing.assert_allclose(sites[8]['predictors'], np.array([0.0, 9.0, 0.0]))
    assert sites[8]['direct']
    

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
    test = {'specie': 'Co', 'position': np.array([ 0.0,  0.0,  0.0]),
            'selective': [True, True, True], 'velocities': None,
            'predictors': None, 'direct': True}
    np.testing.assert_allclose(sites[0]['position'], test['position'])
    assert sites[0]['specie'] == test['specie']
    assert sites[0]['selective'] == [True, True, True]
    assert sites[0]['velocities'] == None
    assert sites[0]['predictors'] == None
    assert sites[0]['direct']
