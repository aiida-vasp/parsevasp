import os
import pytest
import numpy as np
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
    assert stream_parser.number_of_entries == 1
    assert stream_parser.has_entries
    print(entries[0])
    assert str(entries[0]) == '(ERROR) ibzkpt: Error with the k-points'


def test_stream_override(stream_parser):
    """Check that the stream override works."""
    import re
    testdir = os.path.dirname(__file__)
    stream_file = testdir + '/stdout'
    stream = Stream(file_path=stream_file, config={'ibzkpt': {'kind': 'WARNING',
                                                              'regex': 'internal error',
                                                              'message': 'some error',
                                                              'suggestion': 'none',
                                                              'location': 'STDOUT',
                                                              'recover': True}})
    print(stream.entries)
    assert len(stream.entries) == 1
    assert stream.entries[0].kind == 'WARNING'
    assert stream.entries[0].regex == re.compile('internal error')
    assert stream.entries[0].message == 'some error'
    assert stream.entries[0].suggestion == 'none' 
    assert stream.entries[0].location == stream_parser.entries[0].location
    assert stream.entries[0].recover == stream_parser.entries[0].recover
