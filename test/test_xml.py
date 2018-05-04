import os
import pytest
import numpy as np
from parsevasp import vasprun

@pytest.fixture(scope = 'module', params=[0, 1])
def xml_parser(request, tmpdir_factory):
    """Load XML file.

    """
    testdir = os.path.dirname(__file__)
    xmlfile = testdir + "/basic.xml"
    tmpfile = str(tmpdir_factory.mktemp('data').join('basic_trunc.xml'))
    xml_truncate(request.param, xmlfile, tmpfile)
    xml  = vasprun.Xml(tmpfile)
    
    return xml

def xml_truncate(index, original, tmp):
    """Truncate the XML file.

    """
    
    with open(original, 'r') as xmlfile:
        content = xmlfile.read().splitlines()
    truncated_content = '\n'.join(content[:-index or None])
    with open(tmp, 'w') as xmlfile:
        xmlfile.write(str(truncated_content))

    return

def test_xml_exist(xml_parser):
    """Check if xml_parser exists.

    """
    assert xml_parser.get_dict()
    
def test_xml_energies(xml_parser):
    """Check energies.

    """
    
    assert xml_parser.get_energies('initial') == [-43.312106219999997]
    assert xml_parser.get_energies('final') == [-43.312106219999997]
    assert xml_parser.get_energies('all') == [-43.312106219999997, -43.312106219999997]

def test_xml_forces(xml_parser):
    """Check forces.

    """

    forces = xml_parser.get_forces('initial')
    assert np.all(forces == np.array([[ 0., -0., 0.],
                                      [ 0.,  0., -0.],
                                      [-0.,  0., -0.],
                                      [-0.,  0.,  0.],
                                      [-0., -0.,  0.],
                                      [-0., -0.,  0.],
                                      [ 0., -0., -0.],
                                      [-0.,  0., -0.]]))
    forces = xml_parser.get_forces('initial')
    assert np.all(forces == np.array([[ 0., -0., 0.],
                                      [ 0.,  0., -0.],
                                      [-0.,  0., -0.],
                                      [-0.,  0.,  0.],
                                      [-0., -0.,  0.],
                                      [-0., -0.,  0.],
                                      [ 0., -0., -0.],
                                      [-0.,  0., -0.]]))
    forces = xml_parser.get_forces('all')
    assert np.all(forces[1] == np.array([[ 0., -0., 0.],
                                       [ 0.,  0., -0.],
                                       [-0.,  0., -0.],
                                       [-0.,  0.,  0.],
                                       [-0., -0.,  0.],
                                       [-0., -0.,  0.],
                                       [ 0., -0., -0.],

                                         [-0.,  0., -0.]]))
    assert np.all(forces[2] == np.array([[ 0., -0., 0.],
                                       [ 0.,  0., -0.],
                                       [-0.,  0., -0.],
                                       [-0.,  0.,  0.],
                                       [-0., -0.,  0.],
                                       [-0., -0.,  0.],
                                       [ 0., -0., -0.],
                                       [-0.,  0., -0.]]))           

def test_xml_stress(xml_parser):
    """Check stress.

    """

    assert xml_parser.get_stress('initial') == None
    assert xml_parser.get_stress('final') == None
    assert xml_parser.get_stress('all') == {1: None, 2: None}

def test_xml_hessian(xml_parser):
    """Check hessian matrix.

    """

    assert xml_parser.get_hessian() == None

def test_xml_dynmat(xml_parser):
    """Check the dynamical metrix.

    """
    
    assert xml_parser.get_dynmat() == None

def test_xml_dielectrics(xml_parser):
    """Check the dielectric functions.

    """

    assert xml_parser.get_dielectrics() == None

def test_xml_born(xml_parser):
    """Check the born effective masses.

    """
    
    assert xml_parser.get_born() == None

def test_xml_fermi_level(xml_parser):
    """Check the Fermi level.

    """

    assert xml_parser.get_fermi_level() == 5.92134456

def test_xml_unitcell(xml_parser):
    """Check the unitcell.

    """
    
    unitcell = xml_parser.get_unitcell('initial')
    assert np.all(unitcell == np.array([[ 5.46900498, 0.        , 0.        ],
                                        [ 0.        , 5.46900498, 0.        ],
                                        [ 0.        , 0.        , 5.46900498]]))
    unitcell = xml_parser.get_unitcell('final')
    assert np.all(unitcell == np.array([[ 5.46900498, 0.        , 0.        ],
                                        [ 0.        , 5.46900498, 0.        ],
                                        [ 0.        , 0.        , 5.46900498]]))
    unitcell = xml_parser.get_unitcell('all')
    assert np.all(unitcell[1] == np.array([[ 5.46900498, 0.        , 0.        ],
                                           [ 0.        , 5.46900498, 0.        ],
                                           [ 0.        , 0.        , 5.46900498]]))
    assert np.all(unitcell[2] == np.array([[ 5.46900498, 0.        , 0.        ],
                                           [ 0.        , 5.46900498, 0.        ],
                                           [ 0.        , 0.        , 5.46900498]]))

def test_xml_positions(xml_parser):
    """Check the positions.

    """
    
    positions = xml_parser.get_positions('initial')
    assert np.all(positions == np.array([[ 0.        ,  0.        ,  0.        ],
                                         [ 0.        ,  0.50000092,  0.50000092],
                                         [ 0.50000092,  0.50000092,  0.        ],
                                         [ 0.50000092,  0.        ,  0.50000092],
                                         [ 0.75000046,  0.24999954,  0.75000046],
                                         [ 0.24999954,  0.24999954,  0.24999954],
                                         [ 0.24999954,  0.75000046,  0.75000046],
                                         [ 0.75000046,  0.75000046,  0.24999954]]))
    positions = xml_parser.get_positions('final')
    assert np.all(positions == np.array([[ 0.        ,  0.        ,  0.        ],
                                         [ 0.        ,  0.50000092,  0.50000092],
                                         [ 0.50000092,  0.50000092,  0.        ],
                                         [ 0.50000092,  0.        ,  0.50000092],
                                         [ 0.75000046,  0.24999954,  0.75000046],
                                         [ 0.24999954,  0.24999954,  0.24999954],
                                         [ 0.24999954,  0.75000046,  0.75000046],
                                         [ 0.75000046,  0.75000046,  0.24999954]]))
    positions = xml_parser.get_positions('all')
    assert np.all(positions[1] == np.array([[ 0.        ,  0.        ,  0.        ],
                                            [ 0.        ,  0.50000092,  0.50000092],
                                            [ 0.50000092,  0.50000092,  0.        ],
                                            [ 0.50000092,  0.        ,  0.50000092],
                                            [ 0.75000046,  0.24999954,  0.75000046],
                                            [ 0.24999954,  0.24999954,  0.24999954],
                                            [ 0.24999954,  0.75000046,  0.75000046],
                                            [ 0.75000046,  0.75000046,  0.24999954]]))
    assert np.all(positions[2] == np.array([[ 0.        ,  0.        ,  0.        ],
                                            [ 0.        ,  0.50000092,  0.50000092],
                                            [ 0.50000092,  0.50000092,  0.        ],
                                            [ 0.50000092,  0.        ,  0.50000092],
                                            [ 0.75000046,  0.24999954,  0.75000046],
                                            [ 0.24999954,  0.24999954,  0.24999954],
                                            [ 0.24999954,  0.75000046,  0.75000046],
                                            [ 0.75000046,  0.75000046,  0.24999954]]))

def test_xml_species(xml_parser):
    """Check the species.

    """

    species = xml_parser.get_species()
    assert np.all(species == [14, 14, 14, 14, 14, 14, 14, 14])

def test_xml_kpoints(xml_parser):
    """Check the kpoints.

    """
    
    kpoints = xml_parser.get_kpoints()
    assert np.all(kpoints == np.array([[ 0.        ,  0.        ,  0.        ],
                                       [ 0.16666667,  0.        ,  0.        ],
                                       [ 0.33333333,  0.        ,  0.        ],
                                       [ 0.5       ,  0.        ,  0.        ],
                                       [ 0.16666667,  0.16666667,  0.        ],
                                       [ 0.33333333,  0.16666667,  0.        ],
                                       [ 0.5       ,  0.16666667,  0.        ],
                                       [ 0.33333333,  0.33333333,  0.        ],
                                       [ 0.5       ,  0.33333333,  0.        ],
                                       [ 0.5       ,  0.5       ,  0.        ],
                                       [ 0.16666667,  0.16666667,  0.16666667],
                                       [ 0.33333333,  0.16666667,  0.16666667],
                                       [ 0.5       ,  0.16666667,  0.16666667],
                                       [ 0.33333333,  0.33333333,  0.16666667],
                                       [ 0.5       ,  0.33333333,  0.16666667],
                                       [ 0.5       ,  0.5       ,  0.16666667],
                                       [ 0.33333333,  0.33333333,  0.33333333],
                                       [ 0.5       ,  0.33333333,  0.33333333],
                                       [ 0.5       ,  0.5       ,  0.33333333],
                                       [ 0.5       ,  0.5       ,  0.5       ]]))
    
def test_xml_kpointsw(xml_parser):
    """Check the kpoint weights.

    """
    
    kpointsw = xml_parser.get_kpointsw()
    assert np.all(kpointsw == np.array([0.00462963,  0.02777778,  0.02777778,
                                        0.01388889,  0.05555556,  0.11111111,
                                        0.05555556,  0.05555556,  0.05555556,
                                        0.01388889,  0.03703704,  0.11111111,
                                        0.05555556,  0.11111111,  0.11111111,
                                        0.02777778,  0.03703704,  0.05555556,
                                        0.02777778,  0.00462963]))
