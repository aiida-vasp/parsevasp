import os
import pytest
import numpy as np
from parsevasp.outcar import Outcar


@pytest.fixture(scope='module', params=[0])
def outcar_parser(request, tmpdir_factory):
    """Load OUTCAR file.

    """
    testdir = os.path.dirname(__file__)
    outcarfile = testdir + '/OUTCAR'
    tmpfile = str(tmpdir_factory.mktemp('data').join('OUTCAR'))
    outcar_truncate(request.param, outcarfile, tmpfile)
    outcar = Outcar(file_path=tmpfile)

    return outcar


@pytest.fixture(scope='module', params=[0])
def outcar_parser_file_objects(request, tmpdir_factory):
    """Load OUTCAR file from a file object.

    """
    testdir = os.path.dirname(__file__)
    outcarfile = testdir + '/OUTCAR'
    tmpfile = str(tmpdir_factory.mktemp('data').join('OUTCAR'))
    outcar_truncate(request.param, outcarfile, tmpfile)
    outcar = None
    with open(tmpfile) as file_handler:
        outcar = Outcar(file_handler=file_handler)

    return outcar


@pytest.fixture(scope='module', params=[0])
def outcar_parser_magnetization(request, tmpdir_factory):
    """Load OUTCAR file.

    """
    testdir = os.path.dirname(__file__)
    outcarfile = testdir + '/OUTCAR_MAG'
    tmpfile = str(tmpdir_factory.mktemp('data').join('OUTCAR_MAG'))
    outcar_truncate(request.param, outcarfile, tmpfile)
    outcar = Outcar(file_path=tmpfile)

    return outcar


@pytest.fixture(scope='module', params=[0])
def outcar_parser_magnetization_single(request, tmpdir_factory):
    """Load OUTCAR file.

    """
    testdir = os.path.dirname(__file__)
    outcarfile = testdir + '/OUTCAR_MAG_SINGLE'
    tmpfile = str(tmpdir_factory.mktemp('data').join('OUTCAR_MAG_SINGLE'))
    outcar_truncate(request.param, outcarfile, tmpfile)
    outcar = Outcar(file_path=tmpfile)

    return outcar


def outcar_truncate(index, original, tmp):
    """Truncate the OUTCAR file.

    """

    with open(original, 'r') as outcarfile:
        content = outcarfile.read().splitlines()
    truncated_content = '\n'.join(content[:-index or None])
    with open(tmp, 'w') as outcarfile:
        outcarfile.write(str(truncated_content))

    return


def test_outcar_symmetry(outcar_parser):
    """Check if outcar_parser returns correct symmetry entries.

    """

    symmetry = outcar_parser.get_symmetry()
    test = [
        'T_d', 'D_2d.', 'D_2d.', 'D_2d.', 'D_2d.', 'D_2d.', 'D_2d.', 'C_2',
        'C_2', 'C_2', 'C_2', 'C_2', 'C_2', 'C_2v.', 'C_2v.', 'T_d'
    ]
    assert symmetry['site_symmetry_at_origin']['static'] == test
    assert symmetry['site_symmetry_at_origin']['dynamic'] == test
    test = [
        'O_h', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'C_2h.',
        'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'D_2h.', 'D_2h.', 'O_h'
    ]
    assert symmetry['point_group']['static'] == test
    assert symmetry['point_group']['dynamic'] == test
    test = [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48]
    assert symmetry['num_space_group_operations']['static'] == test
    assert symmetry['num_space_group_operations']['dynamic'] == test
    test = [
        'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
        'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
        'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
        'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell'
    ]
    assert symmetry['original_cell_type']['static'] == test
    assert symmetry['original_cell_type']['dynamic'] == test
    test = [
        'face centered cubic supercell.',
        'body centered tetragonal supercell.',
        'body centered tetragonal supercell.',
        'body centered tetragonal supercell.',
        'body centered tetragonal supercell.',
        'body centered tetragonal supercell.',
        'body centered tetragonal supercell.',
        'base centered monoclinic supercell.',
        'base centered monoclinic supercell.',
        'base centered monoclinic supercell.',
        'base centered monoclinic supercell.',
        'base centered monoclinic supercell.',
        'base centered monoclinic supercell.',
        'face centered cubic supercell.', 'face centered cubic supercell.',
        'face centered cubic supercell.'
    ]
    assert symmetry['symmetrized_cell_type']['static'] == test
    assert symmetry['symmetrized_cell_type']['dynamic'] == test
    test = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    assert symmetry['primitive_translations'] == test


def test_outcar_elastic(outcar_parser):
    """Check if outcar_parser returns correct elastic moduli entries.

    """

    elastic = outcar_parser.get_elastic_moduli()
    test = np.array([[
        1.6740702e+03, 7.0419980e+02, 7.0419980e+02, -0.0000000e+00,
        0.0000000e+00, 0.0000000e+00
    ],
                     [
                         7.0502380e+02, 1.6748491e+03, 7.0502380e+02,
                         -0.0000000e+00, -0.0000000e+00, 0.0000000e+00
                     ],
                     [
                         7.0499350e+02, 7.0499350e+02, 1.6748165e+03,
                         0.0000000e+00, -0.0000000e+00, 0.0000000e+00
                     ],
                     [
                         8.2260000e-01, 8.7980000e-01, 1.2896000e+00,
                         1.1225901e+03, -0.0000000e+00, 0.0000000e+00
                     ],
                     [
                         -7.8000000e-03, -4.9500000e-02, 1.4700000e-02,
                         0.0000000e+00, 1.1230829e+03, -0.0000000e+00
                     ],
                     [
                         -2.9200000e-02, -5.3200000e-02, -2.1970000e-01,
                         -0.0000000e+00, 0.0000000e+00, 1.1223147e+03
                     ]])
    np.testing.assert_allclose(elastic['non_symmetrized'], test)
    test = np.array([[1674.5786, 704.739, 704.739, -0., 0., 0.],
                     [704.739, 1674.5786, 704.739, -0., 0., 0.],
                     [704.739, 704.739, 1674.5786, -0., -0., 0.],
                     [-0., -0., -0., 1122.6622, 0., -0.],
                     [0., 0., -0., 0., 1122.6622, -0.],
                     [0., 0., 0., -0., -0., 1122.6622]])
    np.testing.assert_allclose(elastic['symmetrized'], test)
    test = np.array([[1674.5786, 704.739, 704.739, -0., 0., 0.],
                     [704.739, 1674.5786, 704.739, -0., 0., 0.],
                     [704.739, 704.739, 1674.5786, -0., -0., 0.],
                     [-0., -0., -0., 775.8054, 0., -0.],
                     [0., 0., -0., 0., 775.8054, -0.],
                     [0., 0., 0., -0., -0., 775.8054]])
    np.testing.assert_allclose(elastic['total'], test)


def test_outcar_elastic_file_object(outcar_parser_file_objects):
    """Check if outcar_parser returns correct elastic moduli entries.

    """

    elastic = outcar_parser_file_objects.get_elastic_moduli()
    test = np.array([[
        1.6740702e+03, 7.0419980e+02, 7.0419980e+02, -0.0000000e+00,
        0.0000000e+00, 0.0000000e+00
    ],
                     [
                         7.0502380e+02, 1.6748491e+03, 7.0502380e+02,
                         -0.0000000e+00, -0.0000000e+00, 0.0000000e+00
                     ],
                     [
                         7.0499350e+02, 7.0499350e+02, 1.6748165e+03,
                         0.0000000e+00, -0.0000000e+00, 0.0000000e+00
                     ],
                     [
                         8.2260000e-01, 8.7980000e-01, 1.2896000e+00,
                         1.1225901e+03, -0.0000000e+00, 0.0000000e+00
                     ],
                     [
                         -7.8000000e-03, -4.9500000e-02, 1.4700000e-02,
                         0.0000000e+00, 1.1230829e+03, -0.0000000e+00
                     ],
                     [
                         -2.9200000e-02, -5.3200000e-02, -2.1970000e-01,
                         -0.0000000e+00, 0.0000000e+00, 1.1223147e+03
                     ]])
    np.testing.assert_allclose(elastic['non_symmetrized'], test)
    test = np.array([[1674.5786, 704.739, 704.739, -0., 0., 0.],
                     [704.739, 1674.5786, 704.739, -0., 0., 0.],
                     [704.739, 704.739, 1674.5786, -0., -0., 0.],
                     [-0., -0., -0., 1122.6622, 0., -0.],
                     [0., 0., -0., 0., 1122.6622, -0.],
                     [0., 0., 0., -0., -0., 1122.6622]])
    np.testing.assert_allclose(elastic['symmetrized'], test)
    test = np.array([[1674.5786, 704.739, 704.739, -0., 0., 0.],
                     [704.739, 1674.5786, 704.739, -0., 0., 0.],
                     [704.739, 704.739, 1674.5786, -0., -0., 0.],
                     [-0., -0., -0., 775.8054, 0., -0.],
                     [0., 0., -0., 0., 775.8054, -0.],
                     [0., 0., 0., -0., -0., 775.8054]])
    np.testing.assert_allclose(elastic['total'], test)


def test_outcar_magnetization(outcar_parser_magnetization):
    """Check if outcar_magnetization_parser returns the correct magnetization
    """

    magnetization = outcar_parser_magnetization.get_magnetization()
    test = {
        'sphere': {
            'x': {
                'site_moment': {
                    1: {
                        's': -0.014,
                        'p': -0.051,
                        'd': 1.687,
                        'tot': 1.621
                    },
                    2: {
                        's': -0.015,
                        'p': -0.052,
                        'd': 1.686,
                        'tot': 1.619
                    },
                    3: {
                        's': -0.014,
                        'p': -0.053,
                        'd': 1.708,
                        'tot': 1.64
                    },
                    4: {
                        's': -0.014,
                        'p': -0.053,
                        'd': 1.708,
                        'tot': 1.64
                    }
                },
                'total_magnetization': {
                    's': -0.057,
                    'p': -0.21,
                    'd': 6.788,
                    'tot': 6.521
                }
            },
            'y': {
                'site_moment': {},
                'total_magnetization': {}
            },
            'z': {
                'site_moment': {},
                'total_magnetization': {}
            }
        },
        'full_cell': np.asarray([6.4424922])
    }

    for _proj in ['x', 'y', 'z']:

        for _key, _val in test['sphere'][_proj]['site_moment'].items():
            _test = np.asarray(list(_val.values()))
            _mag = np.asarray(
                list(magnetization['sphere'][_proj]['site_moment']
                     [_key].values()))
            np.testing.assert_allclose(_mag, _test)

        _test = np.asarray(
            list(test['sphere'][_proj]['total_magnetization'].values()))
        _mag = np.asarray(
            list(magnetization['sphere'][_proj]
                 ['total_magnetization'].values()))
        np.testing.assert_allclose(_mag, _test)
    _mag = np.asarray(list(magnetization['full_cell']))
    _test = np.asarray(list(test['full_cell']))
    np.testing.assert_allclose(_mag, _test)


def test_outcar_magnetization_single(outcar_parser_magnetization_single):
    """Check if outcar_magnetization_parser returns the correct magnetization
    for a single atom in the unit cell
    """

    magnetization = outcar_parser_magnetization_single.get_magnetization()

    test = {
        'sphere': {
            'x': {
                'site_moment': {
                    1: {
                        's': -0.012,
                        'p': -0.043,
                        'd': 2.49,
                        'tot': 2.434
                    },
                },
                'total_magnetization': {
                    's': -0.012,
                    'p': -0.043,
                    'd': 2.49,
                    'tot': 2.434
                }
            },
            'y': {
                'site_moment': {},
                'total_magnetization': {}
            },
            'z': {
                'site_moment': {},
                'total_magnetization': {}
            }
        },
        'full_cell': np.asarray([2.4077611])
    }

    for _proj in ['x', 'y', 'z']:

        for _key, _val in test['sphere'][_proj]['site_moment'].items():
            _test = np.asarray(list(_val.values()))
            _mag = np.asarray(
                list(magnetization['sphere'][_proj]['site_moment']
                     [_key].values()))
            np.testing.assert_allclose(_mag, _test)

        _test = np.asarray(
            list(test['sphere'][_proj]['total_magnetization'].values()))
        _mag = np.asarray(
            list(magnetization['sphere'][_proj]
                 ['total_magnetization'].values()))
        np.testing.assert_allclose(_mag, _test)
    _mag = np.asarray(list(magnetization['full_cell']))
    _test = np.asarray(list(test['full_cell']))
    np.testing.assert_allclose(_mag, _test)
