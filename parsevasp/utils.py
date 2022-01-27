"""Utiles."""
import logging
import math
import os
import re
import sys

import numpy as np

from parsevasp.base import open_close_file_handler


def read_from_file(file_name, input_file_handler, contains=None, lines=True, logger=None):
    """
    Read a file and return the whole file or specific lines.

    Parameters
    ----------
    file_name : str
        The location and file name to be read.
    file_handler : object
        A valid file handler. If both file name and file_handler is set,
        the file handler takes presence.
    contains : list of str
        A list of string of identifiers for the lines that is to be
        returned. If None, the whole file is returned.
    lines : bool, optional
        If set to False, this method will just return the read() from supplied path or handler.
        Defaults to True.
    logger : object, optional
        A logger object to use.

    Returns
    -------
    parsed : list of str or a str
        If `lines` is True, the list of strings containing the whole or specific
        lines from a file is returned. If `lines` is False, a string of all the content
        is returned.

    """

    if logger is None:
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=W0212

    if input_file_handler is not None:
        inputfile = input_file_handler
        if not lines:
            # Only want a string of file content, return.
            return inputfile.read()
        file_data = inputfile.readlines()
    else:
        inputfile = open_close_file_handler(file_name=file_name, status='r', logger=logger)
        if not lines:
            # Only want a string of file content, return.
            file_data = inputfile.read()
        else:
            file_data = inputfile.readlines()
        open_close_file_handler(logger, file_handler=inputfile, logger=logger)
        if not lines:
            return file_data

    parsed = []

    # first check if contains is a list
    is_list = is_sequence(contains)

    if contains is not None:
        # this can be a bit faster (comprehension), but do not care for this
        # now
        for _, line in enumerate(file_data):
            if is_list:
                for element in contains:
                    if element in line:
                        parsed.append(line)
            else:
                if contains in line:
                    parsed = line
    else:
        parsed = file_data

    return parsed


def file_exists(file_path, logger=None):
    """
    Check if the file exists.

    Parameters
    ----------
    file_path : string
        The file path to be checked.
    logger : object, optional
        A logger object to use.

    Returns
    -------
    status : bool
        If file does not exists or `file_path` empty, else False.
    """
    from parsevasp.base import BaseParser

    if logger is None:
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=W0212

    if not file_path:
        logger.error(BaseParser.ERROR_MESSAGES[BaseParser.ERROR_EMPTY_FILE_PATH])
        sys.exit(BaseParser.ERROR_EMPTY_FILE_PATH)

    status = True
    try:
        os.stat(file_path)
    except OSError:
        logger.error(
            f'{BaseParser.ERROR_MESSAGES[BaseParser.ERROR_FILE_NOT_FOUND]} The file in question is: {file_path}'
        )
        status = False

    return status


def is_sequence(arg):
    """
    Checks to see if something is a sequence (list).

    Parameters
    ----------
    arg : str
        The string to be examined.

    Returns
    -------
    sequence : bool
        Is True if `arg` is a list.

    """

    if not hasattr(arg, 'strip') and hasattr(arg, '__getitem__'):
        return True
    if hasattr(arg, '__iter__'):
        return True
    return False


def test_string_content(string):
    """
    Detects if string is integer, float or string.

    Parameters
    ----------
    string : string
        An input string to be tested.

    Returns
    -------
    string
        A string with value 'int' if input is an integer,
        'float' if the input is a float and 'string' if it
        is just a regular string.

    """
    try:
        float(string)
        return 'int' if ((string.count('.') == 0) and \
                         ('e' not in string) and \
                         ('E' not in string)) else 'float'
    except ValueError:
        return 'string'


def is_numbers(string, splitter=' '):
    """
    Check if a string only contains numbers.

    Parameters
    ----------
    s: str
        The input string
    splitter : string, optional
        The splitting character to be used, defaults to blank spaces.

    Returns
    -------
    is_nums: bool
        Is True if all entries in the input string is a numbers,
        otherwise False.

    """

    entries = string.split(splitter)
    is_nums = True
    for entry in entries:
        if not is_number(entry):
            is_nums = False
            return is_nums

    return is_nums


def is_number(string):
    """
    Check if a string is a number.

    Parameters
    ----------
    s: str
        The input string

    Returns
    -------
    is_num: bool
        Is True if the input string is a number, otherwise False

    """

    try:
        float(string)
        is_num = True
    except ValueError:
        is_num = False

    return is_num


def remove_newline(fobj, num_newlines=1):
    """
    Removes the newline at the end of a file.

    Usefull to run after a for loop that writes a newline character
    at each step. Other solutions cannot handle very large files.

    Parameters
    ----------
    fobj : object
        A file object.
    num_newlines : int, optional
        The number of newlines to remove. Defaults to 1.

    """

    # remove last newline, check number of chars, different
    # for each OS
    remove_chars = len(os.linesep) + num_newlines - 1
    fobj.truncate(fobj.tell() - remove_chars)


def dir_to_cart(vector, lattice):
    """
    Convert direct coordinates to cartesian.

    Parameters
    ----------
    vector : ndarray
        | Dimension: (3)

        The direct vector to be converted.
    lattice : ndarray
        | Dimension: (3,3)

        The crystal lattice, where the first lattice vector is
        [0,:], the second, [1,:] etc.

    Returns
    -------
    cart : ndarray
        | Dimension: (3)

        The cartesian vector.

    """

    cart = np.dot(vector, lattice)

    return cart


def cart_to_dir(vector, lattice):
    """
    Convert cartesian coordinates to direct.

    Parameters
    ----------
    vectir : ndarray
        | Dimension: (3)

        The cartersian vector.
    lattice : ndarray
        | Dimension: (3,3)
        The crystal lattice, where the first lattice vector is
        (0,:), the second, (1,:) etc.

    Returns
    -------
    direct : ndarray
        | Dimension: (3)
        The direct vector.

    """

    direct = np.dot(vector, np.linalg.inv(lattice))

    return direct


def lat_to_reclat(lattice):
    r"""
    Convert the lattice to the reciprocal lattice.

    Parameters
    ----------
    lattice : ndarray
        | Dimension: (3,3)
        The crystal lattice, where the first lattice vector is
        (0,:), the second, (1,:) etc.

    Returns
    -------
    lattice_rec : ndarray
        | Dimension: (3,3)
        Reciprocal lattice including the 2:math:`\pi` factor,
        see `lattice` for layout.

    Notes
    -----
    In general, `lattice_rec`=2pi*(lattice.T)^-1

    """

    lattice_trans = np.transpose(lattice)
    lattice_rec = 2 * math.pi * np.linalg.inv(lattice_trans)

    return lattice_rec


def match_integer_param(inputs, key, string):
    """
    Search a string for a given parameter and set its values, assuming integer.

    Parameters
    ----------
    inputs : dict
      A dictionary containing the parameters we want to set in lowercase and default values.
    key : string
      A string containing the matching parameter we are looking for.
    string: string
      A string to be searched for key.

    Returns
    -------
    value : integer
      The located value.

    """

    match = re.match(r'^ +' + key + r' *= *([-0-9]+)', string)
    if match:
        inputs[key.lower()] = int(match.group(1))


def line_to_type(fobject_or_string, d_type=str, no_split=False):
    """
    Grab a line from a file like object or string and convert it to d_type (default: str).

    Parameters
    ----------
    fobject_or_string : object
        A file like object or a string containing something that is to be converted to a specified type
    dtype : object
        The dtype one want to convert to. The standard Python dtypes are supported.
    no_splot : bool
        If True do not split a string. Useful for comments etc. We still strip.

    """
    if isinstance(fobject_or_string, str):
        line = fobject_or_string
    else:
        line = fobject_or_string.readline()
    # Previously this was map instead of list comprehension
    if not no_split:
        result = [d_type(item) for item in line.split()]
    else:
        result = line.strip()
    if len(result) == 1:
        return result[0]
    return result


empty_line = re.compile(r'[\r\n]\s*[\r\n]')
