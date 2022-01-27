"""Handle CHGCAR."""
import sys

import numpy as np

from parsevasp import utils
from parsevasp.base import BaseParser


class Chgcar(BaseParser):
    """Class to handle CHGCAR."""

    def __init__(self, file_path=None, file_handler=None, logger=None):
        """
        Initialize an CHGCAR object and set content as a dictionary.

        file_path : string, optional
            A string containing the file path to the file that is going to be parsed.
        file_handler : object, optional
            A file like object that acts as a handler for the content to be parsed.
        density : bool, optional
            If True, we divide the read CHGCAR data by the unit cell volume. (not implemented.)
        logger : object, optional
            A logger object if you would like to use an external logger for messages
            ejected inside this parser.

        """

        super(Chgcar, self).__init__(file_path=file_path, file_handler=file_handler, logger=logger)

        # check that at least one is supplied
        if self._file_path is None and self._file_handler is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        if self._file_path is None and self._file_handler is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        self._data = {
            'total': None,
            'magnetization': None,
        }

        # parse parse parse
        self._parse()

    def _write(self, *args, **kwargs):
        """Write not supported for CHGCAR."""
        raise NotImplementedError('Writing CHGCAR files is not supported.')

    def _parse(self):
        """Perform the actual parsing by parsing from file like content."""
        self._from_file()

    def _from_file(self):
        """
        Load CHGCAR into NumPy arrays.

        This method is presently not optimized to use as little memory as possible.

        """
        content = utils.read_from_file(self._file_path, self._file_handler, lines=False)
        # Extract header
        temp = content.split('\n\n', 1)
        header = temp[0]
        content = temp[1]
        header = header.split('\n')
        # comment = header[0]
        scaling = float(header[1])
        lattice_vectors = np.zeros((3, 3))
        for i in range(3):
            # Read and scale lattice vectors
            lattice_vectors[i] = scaling * np.array([float(item) for item in header[i + 2].split()])
            # Calculate volume for later scaling
        volume = 1.0
        if volume:
            volume = np.dot(lattice_vectors[0], np.cross(lattice_vectors[1], lattice_vectors[2]))
        # First line of content should now be NGXF, NGYF, NGZF
        temp = content.split('\n', 1)
        ngf_string = temp[0]
        content = temp[1]
        ngf = [int(item) for item in ngf_string.split()]
        # Need to reverse as CHGCAR is x fastest, while we want
        # to comply with z fastest (C order).
        ngf.reverse()
        # Check how many datasets we have
        content = content.split(ngf_string)
        num_datasets = len(content)
        # First dataset is always there
        self._data['total'] = np.fromstring(content[0].split('augmentation occupancies')[0], dtype=float,
                                            sep=' ').reshape(ngf) / volume
        if num_datasets == 2:
            # Collinear spin
            self._data['magnetization'] = np.fromstring(
                content[1].split('augmentation occupancies')[0], dtype=float, sep=' '
            ).reshape(ngf) / volume
        elif num_datasets == 4:
            # Non-collinear spin
            self._data['magnetization'] = {}
            self._data['magnetization']['x'] = np.fromstring(
                content[1].split('augmentation occupancies')[0], dtype=float, sep=' '
            ).reshape(ngf) / volume
            self._data['magnetization']['y'] = np.fromstring(
                content[2].split('augmentation occupancies')[0], dtype=float, sep=' '
            ).reshape(ngf) / volume
            self._data['magnetization']['z'] = np.fromstring(
                content[3].split('augmentation occupancies')[0], dtype=float, sep=' '
            ).reshape(ngf) / volume

    @property
    def charge_density(self):
        """Return the charge density."""
        return self._data['total']

    @property
    def magnetization_density(self):
        """Return the magnetization density."""
        return self._data['magnetization']
