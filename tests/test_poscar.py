import os
import pytest
import numpy as np
from parsevasp.poscar import Poscar, Site


@pytest.fixture(scope='module')
def poscar_parser(request):
    """Load POSCAR file.

    """

    testdir = os.path.dirname(__file__)
    poscarfile = testdir + '/POSCAR'
    poscar = Poscar(file_path=poscarfile, write_direct=request.param)

    return poscar


@pytest.fixture(scope='module')
def poscar_parser_file_object():
    """Load POSCAR file using a file object.

    """

    testdir = os.path.dirname(__file__)
    poscarfile = testdir + '/POSCAR'
    poscar = None
    with open(poscarfile) as file_handler:
        poscar = Poscar(file_handler=file_handler)

    return poscar


@pytest.fixture(scope='module')
def poscar_parser_names():
    """Load POSCAR file.

    """

    testdir = os.path.dirname(__file__)
    poscarfile = testdir + '/POSCARNAMES'
    poscar = Poscar(file_path=poscarfile)

    return poscar


@pytest.fixture(scope='module')
def poscar_parser_vel():
    """Load POSCAR file.

    """

    testdir = os.path.dirname(__file__)
    poscarfile = testdir + '/POSCARVEL'
    poscar = Poscar(file_path=poscarfile)

    return poscar

@pytest.mark.parametrize('poscar_parser', [[True]], indirect=True)
def test_poscar_exist(poscar_parser):
    """Check if poscar_parser exists.

    """

    assert poscar_parser.get_dict()

@pytest.mark.parametrize('poscar_parser', [[True]], indirect=True)
def test_poscar_write(tmp_path, poscar_parser):
    """Test that the POSCAR writes correctly."""
    poscar = poscar_parser.get_dict()
    # Write it
    poscar_write_path = tmp_path / "POSCAR"
    poscar_parser.write(poscar_write_path)
    # Then load same content again, but from the newly written file
    poscar_reloaded = Poscar(file_path=poscar_write_path).get_dict()
    # Comment entry is modified so need to remove that before comparing
    del poscar['comment']
    del poscar_reloaded['comment']
    # Compare
    assert np.allclose(poscar['unitcell'], poscar_reloaded['unitcell'])
    for index, site in enumerate(poscar['sites']):
        reloaded = poscar_reloaded['sites'][index]
        assert site['specie'] == reloaded['specie']
        assert np.allclose(site['position'], reloaded['position'])
        assert site['selective'] == reloaded['selective']
        assert site['velocities'] == reloaded['velocities']
        assert np.allclose(site['predictors'], reloaded['predictors'])
        assert site['direct'] == reloaded['direct']

@pytest.mark.parametrize('poscar_parser', [False], indirect=True)
def test_poscar_cartesian(tmp_path, poscar_parser):
    """Test that the POSCAR writes and reads positional cartesian coordinates correctly."""
    poscar = poscar_parser.get_dict()
    # Write it
    poscar_write_path = tmp_path / "POSCAR"
    poscar_parser.write(poscar_write_path)
    poscar_cartesian = None
    with open(poscar_write_path) as file_object:
        poscar_cartesian = file_object.readlines()
    assert poscar_cartesian[7] == 'Cartesian\n'
    assert poscar_cartesian[8] == '  2.254110000000   2.254110000000   2.254110000000\n'
    assert poscar_cartesian[9] == '  6.762340000000   6.762340000000   6.762340000000\n'
    assert poscar_cartesian[39] == '  7.529650000000   3.083180000000   4.508230000000\n'
    # Then load same content again, but from the newly written file, which should now be on cartesian form
    poscar_reloaded = Poscar(file_path=poscar_write_path).get_dict()
    # Comment entry is modified so need to remove that before comparing
    del poscar['comment']
    del poscar_reloaded['comment']
    # Compare
    assert np.allclose(poscar['unitcell'], poscar_reloaded['unitcell'])
    for index, site in enumerate(poscar['sites']):
        reloaded = poscar_reloaded['sites'][index]
        assert site['specie'] == reloaded['specie']
        assert np.allclose(site['position'], reloaded['position'])
        assert site['selective'] == reloaded['selective']
        assert site['velocities'] == reloaded['velocities']
        assert np.allclose(site['predictors'], reloaded['predictors'])
        assert site['direct'] == reloaded['direct']    


def test_poscar_entries(poscar_parser_names):
    """Check that the POSCAR do not have to have
    elemental names.

    """

    poscar = poscar_parser_names.get_dict()
    sites = poscar['sites']
    test = {
        'specie': 'Hamburger',
        'position': np.array([0.24999947, 0.24999947, 0.24999947]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    assert sites[0]['specie'] == test['specie']
    test = {
        'specie': 'Pizza',
        'position': np.array([0., 0.33510051, 0.15804985]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    assert sites[8]['specie'] == test['specie']

@pytest.mark.parametrize('poscar_parser', [[True]], indirect=True)
def test_poscar_entries(poscar_parser):
    """Check POSCAR entries.

    """

    poscar = poscar_parser.get_dict()
    assert poscar['comment'] == 'Compound: Co7Sb24.'
    unitcell = poscar['unitcell']
    test = np.array([[9.0164589999999993, 0.,
                      0.], [0., 9.0164589999999993, 0.],
                     [0., 0., 9.0164589999999993]])
    np.testing.assert_allclose(unitcell, test)
    sites = poscar['sites']
    assert len(sites) == 32
    test = {
        'specie': 'Co',
        'position': np.array([0.24999947, 0.24999947, 0.24999947]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[0]['position'], test['position'])
    assert sites[0]['specie'] == test['specie']
    assert sites[0]['selective'] == [True, True, True]
    assert sites[0]['velocities'] == None
    np.testing.assert_allclose(sites[0]['predictors'],
                               np.array([0.0, 0.0, 0.0]))
    assert sites[0]['direct']
    test = {
        'specie': 'Co',
        'position': np.array([0.74999953, 0.24999947, 0.24999947]),
        'selective': [True, True, True],
        'velocties': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[7]['position'], test['position'])
    assert sites[7]['specie'] == test['specie']
    assert sites[7]['selective'] == [True, True, True]
    assert sites[7]['velocities'] == None
    np.testing.assert_allclose(sites[7]['predictors'],
                               np.array([8.0, 0.0, 0.0]))
    assert sites[7]['direct']
    test = {
        'specie': 'Sb',
        'position': np.array([0., 0.33510051, 0.15804985]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[8]['position'], test['position'])
    assert sites[8]['specie'] == test['specie']
    assert sites[8]['selective'] == [True, True, True]
    assert sites[8]['velocities'] == None
    np.testing.assert_allclose(sites[8]['predictors'],
                               np.array([0.0, 9.0, 0.0]))
    assert sites[8]['direct']


def test_poscar_entries_file_object(poscar_parser_file_object):
    """Check POSCAR entries.

    """

    poscar = poscar_parser_file_object.get_dict()
    assert poscar['comment'] == 'Compound: Co7Sb24.'
    unitcell = poscar['unitcell']
    test = np.array([[9.0164589999999993, 0.,
                      0.], [0., 9.0164589999999993, 0.],
                     [0., 0., 9.0164589999999993]])
    np.testing.assert_allclose(unitcell, test)
    sites = poscar['sites']
    assert len(sites) == 32
    test = {
        'specie': 'Co',
        'position': np.array([0.24999947, 0.24999947, 0.24999947]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[0]['position'], test['position'])
    assert sites[0]['specie'] == test['specie']
    assert sites[0]['selective'] == [True, True, True]
    assert sites[0]['velocities'] == None
    np.testing.assert_allclose(sites[0]['predictors'],
                               np.array([0.0, 0.0, 0.0]))
    assert sites[0]['direct']
    test = {
        'specie': 'Co',
        'position': np.array([0.74999953, 0.24999947, 0.24999947]),
        'selective': [True, True, True],
        'velocties': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[7]['position'], test['position'])
    assert sites[7]['specie'] == test['specie']
    assert sites[7]['selective'] == [True, True, True]
    assert sites[7]['velocities'] == None
    np.testing.assert_allclose(sites[7]['predictors'],
                               np.array([8.0, 0.0, 0.0]))
    assert sites[7]['direct']
    test = {
        'specie': 'Sb',
        'position': np.array([0., 0.33510051, 0.15804985]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[8]['position'], test['position'])
    assert sites[8]['specie'] == test['specie']
    assert sites[8]['selective'] == [True, True, True]
    assert sites[8]['velocities'] == None
    np.testing.assert_allclose(sites[8]['predictors'],
                               np.array([0.0, 9.0, 0.0]))
    assert sites[8]['direct']


def test_poscar_entries_vel(poscar_parser_vel):
    """Check POSCAR entries including velocities.

    """

    poscar = poscar_parser_vel.get_dict()
    assert poscar['comment'] == 'Compound: Co7Sb24.'
    unitcell = poscar['unitcell']
    test = np.array([[9.0164589999999993, 0.,
                      0.], [0., 9.0164589999999993, 0.],
                     [0., 0., 9.0164589999999993]])
    np.testing.assert_allclose(unitcell, test)
    sites = poscar['sites']
    assert len(sites) == 32
    test = {
        'specie': 'Co',
        'position': np.array([0.24999947, 0.24999947, 0.24999947]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[0]['position'], test['position'])
    assert sites[0]['specie'] == test['specie']
    assert sites[0]['selective'] == [True, True, True]
    np.testing.assert_allclose(sites[0]['velocities'],
                               np.array([1.0, 0.0, 0.0]))
    np.testing.assert_allclose(sites[0]['predictors'],
                               np.array([0.0, 0.0, 0.0]))
    assert sites[0]['direct']
    test = {
        'specie': 'Co',
        'position': np.array([0.74999953, 0.24999947, 0.24999947]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[7]['position'], test['position'])
    assert sites[7]['specie'] == test['specie']
    assert sites[7]['selective'] == [True, True, True]
    np.testing.assert_allclose(sites[7]['velocities'],
                               np.array([10.0, 1.0, 0.0]))
    np.testing.assert_allclose(sites[7]['predictors'],
                               np.array([8.0, 0.0, 0.0]))
    assert sites[7]['direct']
    test = {
        'specie': 'Sb',
        'position': np.array([0., 0.33510051, 0.15804985]),
        'selective': [True, True, True],
        'velocties': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[8]['position'], test['position'])
    assert sites[8]['specie'] == test['specie']
    assert sites[8]['selective'] == [True, True, True]
    np.testing.assert_allclose(sites[8]['velocities'],
                               np.array([0.0, 0.0, 14.0]))
    np.testing.assert_allclose(sites[8]['predictors'],
                               np.array([0.0, 9.0, 0.0]))
    assert sites[8]['direct']


def test_poscar_entries_dict():
    """Test to check inititialization using dict.

    """

    unitcell = np.array([[9.0164589999999993, 0., 0.],
                         [0., 9.0164589999999993, 0.],
                         [0., 0., 9.0164589999999993]])
    sites = []
    sites.append(
        Site('Co', np.array([0.0, 0.0, 0.0]), [True, True, True], None, None,
             True))
    poscar_dict = {
        'comment': 'Example file',
        'unitcell': unitcell,
        'sites': sites
    }
    poscar_parser = Poscar(poscar_dict=poscar_dict)
    poscar = poscar_parser.get_dict()
    assert poscar['comment'] == 'Example file'
    unitcell = poscar['unitcell']
    test = np.array([[9.0164589999999993, 0.,
                      0.], [0., 9.0164589999999993, 0.],
                     [0., 0., 9.0164589999999993]])
    np.testing.assert_allclose(unitcell, test)
    sites = poscar['sites']
    assert len(sites) == 1
    test = {
        'specie': 'Co',
        'position': np.array([0.0, 0.0, 0.0]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[0]['position'], test['position'])
    assert sites[0]['specie'] == test['specie']
    assert sites[0]['selective'] == [True, True, True]
    assert sites[0]['velocities'] == None
    assert sites[0]['predictors'] == None
    assert sites[0]['direct']


@pytest.mark.parametrize('poscar_parser', [[True]], indirect=True)
def test_poscar_cartesian(poscar_parser):
    """Check that get_dict can return cartesian coordinates.

    """

    poscar = poscar_parser.get_dict(direct=False)
    sites = poscar['sites']
    test = {
        'specie': 'Co',
        'position': np.array([2.25411, 2.25411, 2.25411]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[0]['position'], test['position'])
    test = {
        'specie': 'Co',
        'position': np.array([6.76234, 2.25411, 2.25411]),
        'selective': [True, True, True],
        'velocties': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[7]['position'], test['position'])
    test = {
        'specie': 'Sb',
        'position': np.array([0., 3.02142, 1.42505]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': None,
        'direct': True
    }
    np.testing.assert_allclose(sites[8]['position'], test['position'])
    assert not sites[8]['direct']
