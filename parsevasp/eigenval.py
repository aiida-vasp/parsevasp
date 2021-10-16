#!/usr/bin/python
import sys
import logging
import re
import numpy as np

from parsevasp import utils
from parsevasp.base import BaseParser


class Eigenval(BaseParser):
    def __init__(self,
                 file_path=None,
                 file_handler=None,
                 logger=None,
                 prec=None,
                 conserve_order=False):
        """Initialize an EIGENVAL object and set content as a dictionary."""

        super(Eigenval, self).__init__(file_path=file_path,
                                     file_handler=file_handler,
                                     logger=logger)

        self._conserve_order = conserve_order

        # check that at least one is supplied
        if self._file_path is None and self._file_handler is None:
            self._logger.error(
                self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        if self._file_path is None and self._file_handler is None:
            self._logger.error(
                self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        self._data = {
            'bands': None,
            'kpoints': None,
            'metadata': None
        }

        # parse parse parse
        self._parse()

    def _parse(self):
        """Perform the actual parsing

        Parameters
        ----------
        None

        Returns
        -------
        None

        """

        if self._file_path is None and self._file_handler is None:
            return

        # create dictionary from a file
        self._from_file()

    def _from_file(self):
        """Create a dictionary of entries from a
        file and store them in the this instance's data dictionary.

        """

        eigenval = utils.readlines_from_file(self._file_path, self._file_handler)
        self._from_list(eigenval)

    def _from_list(self, eigenval):
        """Go through the list and extract bands, kpoints and metadata

        Parameters
        ----------
        eigenval : list
            A list of strings containing each line in the EIGENVAL file.

        Returns
        -------
        eigenval_dict : dictionary
            A dictionary containing each EIGENVAL tag as a key with the
            associated element.

        Notes
        -----
        No checking for consistency is done here. We do this at a later step
        in order to be able to keep the input methods as clean as possible.

        """

        # Read metadata
        line_0 = utils.line_to_type(eigenval[0], int)
        line_1 = utils.line_to_type(eigenval[1], float)
        line_2 = utils.line_to_type(eigenval[2], float)
        coord_type = utils.line_to_type(eigenval[3])

        # Read name
        name = utils.line_to_type(eigenval[4])

        # Read number of kpoints and bands
        param_0, num_kp, num_bands = utils.line_to_type(eigenval[5], int)

        # Read the rest of the data
        # Here we convert back to a string since the remainder of this parser
        # segment was taken from an older parser that worked on the string.
        data = ''.join(eigenval[6:])

        # Set rest of metadata
        num_ions, num_atoms, p00, num_spins = line_0

        # Datablocks
        data = re.split(utils.empty_line, data)
        data = [[line.split() for line in block.splitlines()] for block in data]
        kpoints = np.zeros((num_kp, 4))
        bands = np.zeros((num_spins, num_kp, num_bands))
        # Iterate over blocks, pr. k-point
        for k, field in enumerate(data):
            # Remove empty lines
            kpbs = [x for x in field if x]
            # First line in the data block is the kpoint coordinates and weight
            kpi = [float(x) for x in kpbs.pop(0)]
            kpoints[k] = kpi
            # The rest is the band energies
            for point in kpbs:
                bands[:, k, int(point[0]) - 1] = point[1:num_spins + 1]

        # Create the metadata dict
        metadata = {}
        metadata[0] = line_0
        metadata[1] = line_1
        metadata[2] = line_2
        metadata['n_ions'] = num_ions
        metadata['n_atoms'] = num_atoms
        metadata['p00'] = p00
        metadata['nspin'] = num_spins
        metadata['cartesian'] = coord_type.startswith(('c', 'C'))
        metadata['name'] = name
        metadata['some_num'] = param_0
        metadata['n_bands'] = num_bands
        metadata['n_kp'] = num_kp

        # Store
        self._data['metadata'] = metadata
        self._data['bands'] = bands
        self._data['kpoints'] = kpoints

    def get_metadata(self):
        """Return the metadata.

        Parameters
        ----------
        None

        Returns
        -------
        metadata : dict
            A dictionary containing the number of number of atoms, ions, spin flag,
            coordinates etc.  

        """

        metadata = self._data['metadata']
        return metadata

    def get_bands(self):
        """Return the bands.

        Parameters
        ----------
        None

        Returns
        -------
        elastic : nparray
            A numpy array containing the bands. First index is spin, second k-points and the last,
            the band index.

        """

        bands = self._data['bands']
        return bands

    def get_kpoints(self):
        """Return the kpoints.

        Parameters
        ----------
        None

        Returns
        -------
        kpoints : nparray
            A NumPy array containing the kpoints. First index is the k-point number, last the direction plus the
            weight (four in total, last is the weight).

        """

        kpoints = self._data['kpoints']
        return kpoints
