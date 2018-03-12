#!/usr/bin/env python
# python specifics
import sys
import logging
import os

def readlines_from_file(filename, contains=None):
    """ Read a file and return the whole file or specific lines

    Parameters
    ----------
    filename : str
        The location and filename to be read.
    contains : list of str
        A list of string of identifiers for the lines that is to be
        returned. If None, the whole file is returned.

    Returns
    -------
    lines : list of str
        The list of strings containing the whole or specific
        lines from a file.

    """

    inputfile = file_handler(filename, status='r')
    file_data = inputfile.readlines()
    file_handler(file_handler=inputfile)
    lines = []

    # first check if contains is a list
    is_list = is_sequence(contains)

    if contains is not None:
        # this can be a bit faster (comprehension), but do not care for this
        # now
        for line_index, line in enumerate(file_data):
            if is_list:
                for element in contains:
                    if element in line:
                        lines.append(line)
            else:
                if contains in line:
                    lines = line
    else:
        lines = file_data

    return lines


def file_handler(filename="", file_handler=None, status=None):
    """ Open and close files

    Parameters
    ----------
    filename : str, optional
        The name of the file to be handled (defaults to '').
    file_handler : object, optional
        An existing `file` object. If not supplied a file is
        created. Needed for file close, otherwise not.
    status : str, optional
        The string containing the status to write, read, append etc.
        If not supplied, assume file close and `file_handler` need
        to be supplied.

    Returns
    -------
    file_handler : object
        If `status` is supplied
        A `file` object

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)

    if status is None:
        if file_handler is None:
            logger.error("Could not close an empty file handler. Exiting.")
            sys.exit(1)
        file_handler.close()
    else:
        try:
            file_handler = open(filename, status)
            return file_handler
        except:
            logger.error("Could not open " + filename + ". Exiting.")
            sys.exit(1)


def is_sequence(arg):
    """ Checks to see if something is a sequence (list)

    Parameters
    ----------
    arg : str
        The string to be examined.

    Returns
    -------
    sequence : bool
        Is True if `arg` is a list.

    """

    sequence = (not hasattr(arg, "strip") and
                hasattr(arg, "__getitem__") or
                hasattr(arg, "__iter__"))

    return sequence


def is_number(s):
    """ Check if a string is a number

    Parameters
    ----------
    s: str
        The input string

    Returns
    -------
    is_number: bool
        Is True if the input string is a number, otherwise False

    """

    try:
        float(s)
        is_number = True
    except ValueError:
        is_number = False

    return is_number

def remove_newline(fobj):
    """Removes the newline at the end of a file. Usefull
    to run after a for loop that writes a newline character
    at each step. Other solutions cannot handle very large files.

    Parameters
    ----------
    fobj : object
        A file object.

    """

    # remove last newline, check number of chars, different
    # for each OS
    remove_chars = len(os.linesep)
    fobj.truncate(fobj.tell() - remove_chars)

def get_gcd(lst):
    """Get the greater common divider for a list of numbers.
    Float numbers are ignored, but divided by the located gcd
    from the integers in the end

    """

    return
