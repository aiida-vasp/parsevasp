#!/usr/bin/python
import sys
import logging
import os

class BaseParser(object):

    ERROR_USE_ONE_ARGUMENT = 10
    ERROR_NO_ENTRIES = 11
    ERROR_NO_KEY = 12
    ERROR_KEY_INVALID_TYPE = 13
    ERROR_FILE_NOT_FOUND = 14
    ERROR_EMPTY_HANDLER = 15
    ERROR_EMPTY_FILE_PATH = 16
    ERROR_MESSAGES = {ERROR_USE_ONE_ARGUMENT: "Supply only one argument when initializing the parser class.",
                      ERROR_NO_ENTRIES: "There is no 'entries' class attribute.",
                      ERROR_NO_KEY: "The correct key in 'entries' is missing.",
                      ERROR_KEY_INVALID_TYPE: "The key has a wrong type.",
                      ERROR_FILE_NOT_FOUND: "The path did not contain an file.",
                      ERROR_EMPTY_HANDLER: "The supplied file handler is empty.",
                      ERROR_EMPTY_FILE_PATH: "The supplied file path is empty."
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
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_FILE_NOT_FOUND] + " The file requested from "
                               "path " + file_path + " was not found.")
            sys.exit(self.ERROR_FILE_NOT_FOUND)
