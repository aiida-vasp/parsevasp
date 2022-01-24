#!/usr/bin/python
import sys
import logging
import os

from parsevasp import utils


class BaseParser:

    ERROR_USE_ONE_ARGUMENT = 10
    ERROR_NO_ENTRIES = 11
    ERROR_NO_KEY = 12
    ERROR_KEY_INVALID_TYPE = 13
    ERROR_FILE_NOT_FOUND = 14
    ERROR_EMPTY_HANDLER = 15
    ERROR_EMPTY_FILE_PATH = 16
    ERROR_MESSAGES = {
        ERROR_USE_ONE_ARGUMENT:
        'Supply only one argument when initializing the parser class.',
        ERROR_NO_ENTRIES: "There is no 'entries' class attribute.",
        ERROR_NO_KEY: "The correct key in 'entries' is missing.",
        ERROR_KEY_INVALID_TYPE: 'The key has a wrong type.',
        ERROR_FILE_NOT_FOUND: 'The path did not contain a file.',
        ERROR_EMPTY_HANDLER: 'The supplied file handler is empty.',
        ERROR_EMPTY_FILE_PATH: 'The supplied file path is empty.'
    }

    def __init__(self, file_path=None, file_handler=None, logger=None):
        """Initialize a general parser object. Used as a base class for the specific parser classes."

        Parameters
        ----------
        file_path : string, optional
            The file path in which the INCAR is read.
        file_hander: object
            A valid file handler object.
        logger : object, optional
            A standard Python logger object.

        """

        self._file_path = file_path
        self._file_handler = file_handler

        # set logger
        if logger is not None:
            self._logger = logger
        else:
            logging.basicConfig(level=logging.DEBUG)
            self._logger = logging.getLogger('ParsevaspParser')


    def write(self, **kwargs):
        """Write respective content as files using a path or handler.

        Parameters
        ----------
        file_path : str, optional
            A string containing the file path to the file that is going to be parsed.
        file_handler : object, optional
            A file like object that acts as a handler for the content to be parsed.
        other : optional
            Any other argument than file path or handler is passed to the specific
            `_write` function.

        One has to provide either a file path or a file handler.

        """

        # Check that we only supply either or of path and handler.
        if ('file_path' in kwargs and 'file_handler' in kwargs) or \
           ('file_path' not in kwargs and 'file_handler' not in kwargs):
            self._logger.error(
                self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        file_path = kwargs.pop('file_path', '')
        file_handler = kwargs.pop('file_handler', '')
        if file_path:
            # Open file
            file_handler = utils.file_handler(file_path, status='w', logger=self._logger)

        # Do the write for each specific content parser _write function using handler, also
        # bring any extra arguments.
        self._write(file_handler, **kwargs)

        if file_path:
            # Close file
            utils.file_handler(file_handler=file_handler, logger=self._logger)


    def _check_file(self, file_path=None):
        """
        Check if a file exists

        Parameters
        ----------
        file_path : string, optional
            The path of the file to be checked. If not supplied, the file path set will be used.

        Returns
        -------
        None

        """
        if file_path is None:
            file_path = self._file_path

        if not os.path.isfile(file_path):
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_FILE_NOT_FOUND] +
                               ' The file requested from '
                               'path ' + file_path + ' was not found.')
            sys.exit(self.ERROR_FILE_NOT_FOUND)
