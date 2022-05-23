"""Base class to handle VASP files."""
# pylint: disable=consider-using-with
import logging
import os
import sys
from abc import ABC, abstractmethod


class BaseParser(ABC):  # pylint: disable=R0903
    """Base class to handle VASP files."""

    ERROR_USE_ONE_ARGUMENT = 10
    ERROR_NO_ENTRIES = 11
    ERROR_NO_KEY = 12
    ERROR_KEY_INVALID_TYPE = 13
    ERROR_FILE_NOT_FOUND = 14
    ERROR_EMPTY_HANDLER = 15
    ERROR_EMPTY_FILE_PATH = 16
    ERROR_MESSAGES = {
        ERROR_USE_ONE_ARGUMENT: 'Supply only one argument when initializing the parser class.',
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
            self._logger = self._setup_logger(logging.DEBUG)

    def _setup_logger(self, level):
        """Setup a logger for this class"""
        logger = logging.getLogger(self.__module__ + '.' + self.__class__.__name__)
        logger.setLevel(level)
        if not logger.handlers:
            handler = logging.StreamHandler()
            handler.setLevel(level)
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        return logger

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
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        file_path = kwargs.pop('file_path', '')
        file_handler = kwargs.pop('file_handler', '')
        if file_path:
            # Open file
            file_handler = open_close_file_handler(file_path, status='w', logger=self._logger)

        # Do the write for each specific content parser _write function using handler, also
        # bring any extra arguments.
        self._write(file_handler, **kwargs)

        if file_path:
            # Close file
            open_close_file_handler(file_handler=file_handler, logger=self._logger)

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
            self._logger.error(
                f'{self.ERROR_MESSAGES[self.ERROR_FILE_NOT_FOUND]} The file requested from '
                f'path {file_path} was not found.'
            )
            sys.exit(self.ERROR_FILE_NOT_FOUND)

    @abstractmethod
    def _write(self, file_handler, **kwargs):
        pass


def open_close_file_handler(file_name='', file_handler=None, status=None, encoding='utf8', logger=None):
    """
    Open and close files.

    Parameters
    ----------
    file_name : str, optional
        The name of the file to be handled (defaults to '').
    file_handler : object, optional
        An existing `file` object. If not supplied a file is
        created. Needed for file close, otherwise not.
    status : str, optional
        The string containing the status to write, read, append etc.
        If not supplied, assume file close and `file_handler` need
        to be supplied.
    encoding : str, optional
        Specify the encoding. Defaults to utf8.
    logger : object, optional
        A logger object to use.

    Returns
    -------
    file_handler : object
        If `status` is supplied
        A `file` object

    """

    if logger is None:
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=W0212

    if status is None:
        if file_handler is None:
            logger.error(BaseParser.ERROR_MESSAGES[BaseParser.ERROR_EMPTY_HANDLER])
            sys.exit(BaseParser.ERROR_EMPTY_HANDLER)
        file_handler.close()
    else:
        try:
            file_handler = open(file_name, status, encoding=encoding)
            return file_handler
        except IOError:
            logger.error(
                f'{BaseParser.ERROR_MESSAGES[BaseParser.ERROR_FILE_NOT_FOUND]} The file in question is: {file_name}'
            )
            sys.exit(BaseParser.ERROR_FILE_NOT_FOUND)
