import os
import pytest
import numpy as np
from parsevasp.kpoints import Kpoints, Kpoint
import utils

@pytest.fixture(scope = 'module')
def kpoints_parser_auto():
    """Load KPOINTS file.

    """

    testdir = os.path.dirname(__file__)
    kpointsfile = testdir + "/KPOINTS"
    kpoints = Kpoints(file_path = kpointsfile)
    
    return kpoints

@pytest.fixture(scope = 'module')
def kpoints_parser_auto_file_object():
    """Load KPOINTS file.

    """

    testdir = os.path.dirname(__file__)
    kpointsfile = testdir + "/KPOINTS"
    kpoints = None
    with open(kpointsfile) as file_handler:
        kpoints = Kpoints(file_handler=file_handler)
    
    return kpoints


@pytest.fixture(scope = 'module')
def kpoints_parser_explicit():
    """Load KPOINTS file.

    """

    testdir = os.path.dirname(__file__)
    kpointsfile = testdir + "/KPOINTSEXP"
    kpoints = Kpoints(file_path = kpointsfile)
    
    return kpoints

@pytest.fixture(scope = 'module')
def kpoints_parser_line():
    """Load KPOINTS file.

    """

    testdir = os.path.dirname(__file__)
    kpointsfile = testdir + "/KPOINTSLINE"
    kpoints = Kpoints(file_path = kpointsfile)
    
    return kpoints

def test_kpoints_exist(kpoints_parser_auto):
    """Check if kpoints_parser exists.

    """

    assert kpoints_parser_auto.get_dict()

def test_kpoints_params_auto(kpoints_parser_auto):
    """Check parameters in KPOINTS for automatic generation.

    """
    
    kpoints = kpoints_parser_auto.get_dict()
    assert kpoints['mode'] == 'automatic'
    assert kpoints['comment'] == 'Example file'
    assert kpoints['divisions'] == [4, 4, 4]
    assert kpoints['shifts'] == None
    assert kpoints['points'] == None
    assert kpoints['centering'] == 'Gamma'
    assert kpoints['tetra'] == None
    assert kpoints['tetra_volume'] == None
    assert kpoints['num_kpoints'] == 0    

def test_kpoints_params_auto_file_object(kpoints_parser_auto_file_object):
    """Check parameters in KPOINTS for automatic generation using a file object.

    """
    
    kpoints = kpoints_parser_auto_file_object.get_dict()
    assert kpoints['mode'] == 'automatic'
    assert kpoints['comment'] == 'Example file'
    assert kpoints['divisions'] == [4, 4, 4]
    assert kpoints['shifts'] == None
    assert kpoints['points'] == None
    assert kpoints['centering'] == 'Gamma'
    assert kpoints['tetra'] == None
    assert kpoints['tetra_volume'] == None
    assert kpoints['num_kpoints'] == 0
    
def test_kpoints_params_explicit(kpoints_parser_explicit):
    """Check parameters in KPOINTS for explicit generation.

    """
    
    kpoints = kpoints_parser_explicit.get_dict()
    assert kpoints['mode'] == 'explicit'
    assert kpoints['comment'] == 'Example file'
    assert kpoints['divisions'] == None
    assert kpoints['shifts'] == None
    assert kpoints['centering'] == None
    points = kpoints['points']
    assert len(points) == 4
    np.testing.assert_allclose(points[0][0], np.array([0.0, 0.0, 0.0]))
    np.testing.assert_allclose(points[1][0], np.array([0.0, 0.0, 0.5]))
    np.testing.assert_allclose(points[2][0], np.array([0.0, 0.5, 0.5]))
    np.testing.assert_allclose(points[3][0], np.array([0.5, 0.5, 0.5]))
    assert utils.isclose(points[0][1], 1.0)
    assert utils.isclose(points[1][1], 1.0)
    assert utils.isclose(points[2][1], 2.0)
    assert utils.isclose(points[3][1], 4.0)
    assert points[0][2]
    assert points[1][2]
    assert points[2][2]
    assert points[3][2]
    assert kpoints['tetra'] == [[6, 1, 2, 3, 4]]
    assert utils.isclose(kpoints['tetra_volume'], 0.183333333333333)

def test_kpoints_params_line(kpoints_parser_line):
    """Check parameters in KPOINTS for line generation.

    """
    
    kpoints = kpoints_parser_line.get_dict()
    assert kpoints['mode'] == 'line'
    assert kpoints['comment'] == 'k-points along high symmetry lines'
    assert kpoints['divisions'] == None
    assert kpoints['shifts'] == None
    assert kpoints['centering'] == None
    assert kpoints['num_kpoints'] == 40
    points = kpoints['points']
    np.testing.assert_allclose(points[0][0], np.array([0.0, 0.0, 0.0]))
    np.testing.assert_allclose(points[1][0], np.array([0.5, 0.5, 0.0]))
    np.testing.assert_allclose(points[2][0], np.array([0.5, 0.5, 0.0]))
    np.testing.assert_allclose(points[3][0], np.array([0.5, 0.75, 0.25]))
    np.testing.assert_allclose(points[4][0], np.array([0.5, 0.75, 0.25]))
    np.testing.assert_allclose(points[5][0], np.array([0.0, 0.0, 0.0]))
    assert utils.isclose(points[0][1], 1.0)
    assert utils.isclose(points[1][1], 1.0)
    assert utils.isclose(points[2][1], 1.0)
    assert utils.isclose(points[3][1], 1.0)
    assert utils.isclose(points[4][1], 1.0)
    assert utils.isclose(points[5][1], 1.0)
    assert points[0][2]
    assert points[1][2]
    assert points[2][2]
    assert points[3][2]
    assert points[4][2]
    assert points[5][2]
    
def test_kpoints_write_auto(kpoints_parser_auto, tmpdir):
    """Test read, write and read KPOINTS in auto mode.

    """
    
    kpoints = kpoints_parser_auto.get_dict()
    temp_file = str(tmpdir.join('KPOINTS'))
    kpoints_parser_auto.write(file_path = temp_file)
    kpoints_parser_auto_temp = Kpoints(file_path = temp_file)
    kpoints_temp = kpoints_parser_auto_temp.get_dict()
    assert kpoints_temp['mode'] == 'automatic'
    assert kpoints_temp['comment'] == 'Example file'
    assert kpoints_temp['divisions'] == [4, 4, 4]
    assert kpoints_temp['shifts'] == [0.0, 0.0, 0.0]
    assert kpoints_temp['points'] == None
    assert kpoints_temp['centering'] == 'Gamma'
    assert kpoints_temp['tetra'] == None
    assert kpoints_temp['tetra_volume'] == None
    assert kpoints_temp['num_kpoints'] == 0

def test_kpoints_write_explicit(kpoints_parser_explicit, tmpdir):
    """Test read, write and read KPOINTS in explicit mode.

    """
    
    kpoints = kpoints_parser_explicit.get_dict()
    temp_file = str(tmpdir.join('KPOINTSEXP'))
    kpoints_parser_explicit.write(file_path = temp_file)
    kpoints_parser_explicit_temp = Kpoints(file_path = temp_file)
    kpoints_temp = kpoints_parser_explicit_temp.get_dict()
    assert kpoints_temp['mode'] == 'explicit'
    assert kpoints_temp['comment'] == 'Example file'
    assert kpoints_temp['divisions'] == None
    assert kpoints_temp['shifts'] == None
    assert kpoints_temp['centering'] == None
    points = kpoints_temp['points']
    assert len(points) == 4
    np.testing.assert_allclose(points[0][0], np.array([0.0, 0.0, 0.0]))
    np.testing.assert_allclose(points[1][0], np.array([0.0, 0.0, 0.5]))
    np.testing.assert_allclose(points[2][0], np.array([0.0, 0.5, 0.5]))
    np.testing.assert_allclose(points[3][0], np.array([0.5, 0.5, 0.5]))
    assert utils.isclose(points[0][1], 1.0)
    assert utils.isclose(points[1][1], 1.0)
    assert utils.isclose(points[2][1], 2.0)
    assert utils.isclose(points[3][1], 4.0)
    assert points[0][2]
    assert points[1][2]
    assert points[2][2]
    assert points[3][2]
    assert kpoints_temp['tetra'] == [[6, 1, 2, 3, 4]]
    assert utils.isclose(kpoints_temp['tetra_volume'], 0.183333333333333)

def test_kpoints_write_line(kpoints_parser_line, tmpdir):
    """Test read, write and read KPOINTS in line mode.

    """
    
    kpoints = kpoints_parser_line.get_dict()
    temp_file = str(tmpdir.join('KPOINTSLINE'))
    kpoints_parser_line.write(file_path = temp_file)
    kpoints_parser_line_temp = Kpoints(file_path = temp_file)
    kpoints_temp = kpoints_parser_line_temp.get_dict()
    assert kpoints_temp['mode'] == 'line'
    assert kpoints_temp['comment'] == 'k-points along high symmetry lines'
    assert kpoints_temp['divisions'] == None
    assert kpoints_temp['shifts'] == None
    assert kpoints_temp['centering'] == None
    assert kpoints_temp['num_kpoints'] == 40
    points = kpoints_temp['points']
    np.testing.assert_allclose(points[0][0], np.array([0.0, 0.0, 0.0]))
    np.testing.assert_allclose(points[1][0], np.array([0.5, 0.5, 0.0]))
    np.testing.assert_allclose(points[2][0], np.array([0.5, 0.5, 0.0]))
    np.testing.assert_allclose(points[3][0], np.array([0.5, 0.75, 0.25]))
    np.testing.assert_allclose(points[4][0], np.array([0.5, 0.75, 0.25]))
    np.testing.assert_allclose(points[5][0], np.array([0.0, 0.0, 0.0]))
    assert utils.isclose(points[0][1], 1.0)
    assert utils.isclose(points[1][1], 1.0)
    assert utils.isclose(points[2][1], 1.0)
    assert utils.isclose(points[3][1], 1.0)
    assert utils.isclose(points[4][1], 1.0)
    assert utils.isclose(points[5][1], 1.0)
    assert points[0][2]
    assert points[1][2]
    assert points[2][2]
    assert points[3][2]
    assert points[4][2]
    assert points[5][2]

def test_kpoints_modify_auto(kpoints_parser_auto, tmpdir):
    """Test read, modify, write and read KPOINTS in auto mode.

    """
    
    kpoints = kpoints_parser_auto.get_dict()
    assert kpoints['comment'] == 'Example file'
    assert kpoints['divisions'] == [4, 4, 4]
    kpoints_parser_auto.modify('comment', 'No comment')
    kpoints_parser_auto.modify('divisions', [5, 5, 5])
    temp_file = str(tmpdir.join('KPOINTS'))
    kpoints_parser_auto.write(file_path = temp_file)
    kpoints_parser_auto_temp = Kpoints(file_path = temp_file)
    kpoints_temp = kpoints_parser_auto_temp.get_dict()
    assert kpoints_temp['comment'] == 'No comment'
    assert kpoints_temp['divisions'] == [5, 5, 5]
    
def test_kpoints_modify_explicit(kpoints_parser_explicit, tmpdir):
    """Test read, modify, write and read KPOINTS in explicit mode.

    """
    
    kpoints = kpoints_parser_explicit.get_dict()
    assert kpoints['comment'] == 'Example file'
    points = kpoints['points']
    assert len(points) == 4
    np.testing.assert_allclose(points[0][0], np.array([0.0, 0.0, 0.0]))
    np.testing.assert_allclose(points[1][0], np.array([0.0, 0.0, 0.5]))
    np.testing.assert_allclose(points[2][0], np.array([0.0, 0.5, 0.5]))
    np.testing.assert_allclose(points[3][0], np.array([0.5, 0.5, 0.5]))
    assert utils.isclose(points[0][1], 1.0)
    assert utils.isclose(points[1][1], 1.0)
    assert utils.isclose(points[2][1], 2.0)
    assert utils.isclose(points[3][1], 4.0)
    assert points[0][2]
    assert points[1][2]
    assert points[2][2]
    assert points[3][2]
    kpoints_parser_explicit.modify('comment', 'Nada comment')
    point = Kpoint(np.array([0.0, 0.0, 0.0]), 1.0)
    kpoints_parser_explicit.modify('points', point, point_number = 3)
    temp_file = str(tmpdir.join('KPOINTSEXP'))
    kpoints_parser_explicit.write(file_path = temp_file)
    kpoints_parser_explicit_temp = Kpoints(file_path = temp_file)
    kpoints_temp = kpoints_parser_explicit_temp.get_dict()
    assert kpoints_temp['comment'] == 'Nada comment'
    points = kpoints_temp['points']
    assert len(points) == 4
    np.testing.assert_allclose(points[0][0], np.array([0.0, 0.0, 0.0]))
    np.testing.assert_allclose(points[1][0], np.array([0.0, 0.0, 0.5]))
    np.testing.assert_allclose(points[2][0], np.array([0.0, 0.5, 0.5]))
    np.testing.assert_allclose(points[3][0], np.array([0.0, 0.0, 0.0]))
    assert utils.isclose(points[0][1], 1.0)
    assert utils.isclose(points[1][1], 1.0)
    assert utils.isclose(points[2][1], 2.0)
    assert utils.isclose(points[3][1], 1.0)
    assert points[0][2]
    assert points[1][2]
    assert points[2][2]
    assert points[3][2]
    
def test_kpoints_modify_line(kpoints_parser_line, tmpdir):
    """Test read, modify, write and read KPOINTS in line mode.

    """
    
    kpoints = kpoints_parser_line.get_dict()
    assert kpoints['comment'] == 'k-points along high symmetry lines'
    kpoints_parser_line.modify('comment', 'No comment')
    point = Kpoint(np.array([0.5, 0.5, 0.25]), 1.0)
    kpoints_parser_line.modify('points', point, point_number = 3)
    kpoints_parser_line.modify('points', point, point_number = 4)    
    temp_file = str(tmpdir.join('KPOINTSLINE'))
    kpoints_parser_line.write(file_path = temp_file)
    kpoints_parser_line_temp = Kpoints(file_path = temp_file)
    kpoints_temp = kpoints_parser_line_temp.get_dict()
    assert kpoints_temp['comment'] == 'No comment'
    points = kpoints_temp['points']
    np.testing.assert_allclose(points[0][0], np.array([0.0, 0.0, 0.0]))
    np.testing.assert_allclose(points[1][0], np.array([0.5, 0.5, 0.0]))
    np.testing.assert_allclose(points[2][0], np.array([0.5, 0.5, 0.0]))
    np.testing.assert_allclose(points[3][0], np.array([0.5, 0.5, 0.25]))
    np.testing.assert_allclose(points[4][0], np.array([0.5, 0.5, 0.25]))
    np.testing.assert_allclose(points[5][0], np.array([0.0, 0.0, 0.0]))
    assert utils.isclose(points[0][1], 1.0)
    assert utils.isclose(points[1][1], 1.0)
    assert utils.isclose(points[2][1], 1.0)
    assert utils.isclose(points[3][1], 1.0)
    assert utils.isclose(points[4][1], 1.0)
    assert utils.isclose(points[5][1], 1.0)
    assert points[0][2]
    assert points[1][2]
    assert points[2][2]
    assert points[3][2]
    assert points[4][2]
    assert points[5][2]

def test_kpoints_string(tmpdir):
    """Test to initialize KPOINTS in auto mode using string.

    """

    kpoints_str = "# Example file\n0\nG\n4 4 4\n"
    temp_file = str(tmpdir.join('KPOINTS'))
    kpoints_parser_auto_temp = Kpoints(kpoints_string = kpoints_str)
    kpoints_temp = kpoints_parser_auto_temp.get_dict()
    assert kpoints_temp['mode'] == 'automatic'
    assert kpoints_temp['comment'] == 'Example file'
    assert kpoints_temp['divisions'] == [4, 4, 4]
    assert kpoints_temp['shifts'] == None
    assert kpoints_temp['points'] == None
    assert kpoints_temp['centering'] == 'Gamma'
    assert kpoints_temp['tetra'] == None
    assert kpoints_temp['tetra_volume'] == None
    assert kpoints_temp['num_kpoints'] == 0


def test_kpoints_dict(tmpdir):
    """Test to initialize KPOINTS in auto mode using dictionary.

    """
    
    kpoints_dict = {'comment': 'Example file', 'divisions': [5, 5, 5],
                    'mode': 'automatic', 'shifts': None, 'points': None,
                    'centering': 'Gamma', 'tetra': None,
                    'tetra_volume': None, 'num_kpoints': 0}
    temp_file = str(tmpdir.join('KPOINTS'))
    kpoints_parser_auto_temp = Kpoints(kpoints_dict = kpoints_dict)
    kpoints_temp = kpoints_parser_auto_temp.get_dict()
    assert kpoints_temp['mode'] == 'automatic'
    assert kpoints_temp['comment'] == 'Example file'
    assert kpoints_temp['divisions'] == [5, 5, 5]
    assert kpoints_temp['shifts'] == None
    assert kpoints_temp['points'] == None
    assert kpoints_temp['centering'] == 'Gamma'
    assert kpoints_temp['tetra'] == None
    assert kpoints_temp['tetra_volume'] == None
    assert kpoints_temp['num_kpoints'] == 0
