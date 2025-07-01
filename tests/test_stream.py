"""Test stream."""

import os
import pathlib
import re

import pytest

from parsevasp.stream import Stream

cwd = pathlib.Path(__file__).parent


@pytest.fixture()
def stream_parser(request, tmpdir_factory):
    """Load a stream."""
    try:
        name = request.param
    except AttributeError:
        # Test not parametrized
        name = 'stdout'
    testdir = os.path.dirname(__file__)
    stream_file = testdir + '/' + name
    stream = Stream(file_path=stream_file)

    return stream


@pytest.fixture()
def stream_parser_file_objects(request, tmpdir_factory):
    """Load stream file from a file object."""
    try:
        name = request.param
    except AttributeError:
        # Test not parametrized
        name = 'stdout'
    testdir = os.path.dirname(__file__)
    stream_file = testdir + '/' + name
    stream = None
    with open(stream_file) as file_handler:
        stream = Stream(file_handler=file_handler)

    return stream


def test_stream(stream_parser):
    """Check if stream_parser returns expected results."""

    entries = stream_parser.entries
    assert stream_parser.configured_streams
    print(entries)
    assert stream_parser.number_of_entries == 1
    assert stream_parser.has_entries
    assert str(entries[0]) == '(ERROR) ibzkpt: Error with the k-points.'


def test_stream_objects(stream_parser_file_objects):
    """Check if stream_parser_file_objects returns expected results passing an object"""

    entries = stream_parser_file_objects.entries
    assert stream_parser_file_objects.configured_streams
    assert stream_parser_file_objects.number_of_entries == 1
    assert stream_parser_file_objects.has_entries
    assert str(entries[0]) == '(ERROR) ibzkpt: Error with the k-points.'


@pytest.mark.parametrize('stream_parser', (['stdout_nostart']), indirect=['stream_parser'])
def test_executed(stream_parser):
    """Check if stream_parser returns expected results for execution checks."""

    entries = stream_parser.entries
    assert entries[0].shortname == 'nostart'


def test_stream_override(stream_parser):
    """Check that the stream override works."""
    testdir = os.path.dirname(__file__)
    stream_file = testdir + '/stdout'
    stream = Stream(
        file_path=stream_file,
        config={
            'ibzkpt': {
                'kind': 'WARNING',
                'regex': 'internal error',
                'message': 'some error',
                'suggestion': 'none',
                'location': 'STDOUT',
                'recover': True,
            }
        },
    )
    assert len(stream.entries) == 1
    assert stream.entries[0].kind == 'WARNING'
    assert stream.entries[0].regex == re.compile('internal error')
    assert stream.entries[0].message == 'some error'
    assert stream.entries[0].suggestion == 'none'
    assert stream.entries[0].location == stream_parser.entries[0].location
    assert stream.entries[0].recover == stream_parser.entries[0].recover


def test_stream_zbrent():
    """Check parsing ZBRENT error."""
    stream_file = cwd / 'stdout_ZBRENT'
    stream = Stream(file_path=stream_file)
    assert stream.entries[0].kind == 'ERROR'
    assert stream.entries[0].message == 'Error in ZBRENT'


def test_stream_inconsistent_lattice():
    """Check parsing inconsistent lattice type error in VAPS6 style"""
    stream_file = cwd / 'stdout_inconsistent_lattice'
    stream = Stream(file_path=stream_file)
    assert stream.entries[0].kind == 'ADVICE'
    assert stream.entries[1].kind == 'ERROR'
    assert 'Inconsistent Bravais lattice types found for crystalline and' in stream.entries[1].message
    assert len(stream.entries) == 2


def test_generic_box_error():
    """Check parsing error box in VAPS6 style"""
    stream_file = cwd / 'stdout_generic_box_error'
    stream = Stream(file_path=stream_file)
    assert stream.entries[1].kind == 'ERROR'
    assert stream.entries[1].regex.pattern == 'Bla bla bla'


def test_generic_warning_error():
    """Check parsing warning box in VAPS6 style"""
    stream_file = cwd / 'stdout_warning'
    stream = Stream(file_path=stream_file)
    assert stream.entries[0].kind == 'ADVICE'
    assert stream.entries[0].regex.pattern == 'You enforced a specific xc type in the INCAR file but a different'
    assert stream.entries[0].shortname == 'xc_enforced'
    assert stream.entries[1].kind == 'WARNING'
    assert stream.entries[1].regex.pattern == 'You use a magnetic or noncollinear calculation, but did not specify'


def test_generic_bug_error():
    """Check parsing bug box in VAPS6 style"""
    stream_file = cwd / 'stdout_bug'
    stream = Stream(file_path=stream_file)
    assert any(entry.kind == 'BUG' for entry in stream.entries)
    bug_entry = next(entry for entry in stream.entries if entry.kind == 'BUG')
    assert bug_entry.regex.pattern == 'internal error in: radial.F  at line: 844'
