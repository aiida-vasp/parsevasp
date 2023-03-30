"""Test stream."""
import os

import numpy as np
import pytest

from parsevasp.stream import Stream


@pytest.fixture(scope='module', params=[0])
def stream_parser(request, tmpdir_factory):
    """Load a stream.

    """
    testdir = os.path.dirname(__file__)
    stream_file = testdir + '/stdout'
    stream = Stream(file_path=stream_file)

    return stream


@pytest.fixture(scope='module', params=[0])
def stream_parser_file_objects(request, tmpdir_factory):
    """Load stream file from a file object.

    """
    testdir = os.path.dirname(__file__)
    stream_file = testdir + '/stdout'
    stream = None
    with open(stream_file) as file_handler:
        stream = Stream(file_handler=file_handler)

    return stream


def test_stream(stream_parser):
    """Check if stream_parser returns expected results.

    """

    entries = stream_parser.entries
    assert stream_parser.configured_streams
    assert stream_parser.number_of_entries == 2
    assert stream_parser.has_entries
    assert str(entries[1]) == '(ERROR) ibzkpt: Error with the k-points.'


def test_executed(stream_parser):
    """Check if stream_parser returns expected results for execution checks.

    """

    entries = stream_parser.entries
    assert entries[0].shortname == 'started'


def test_stream_override(stream_parser):
    """Check that the stream override works."""
    import re
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
                'recover': True
            }
        }
    )
    assert len(stream.entries) == 2
    assert stream.entries[1].kind == 'WARNING'
    assert stream.entries[1].regex == re.compile('internal error')
    assert stream.entries[1].message == 'some error'
    assert stream.entries[1].suggestion == 'none'
    assert stream.entries[1].location == stream_parser.entries[1].location
    assert stream.entries[1].recover == stream_parser.entries[1].recover
