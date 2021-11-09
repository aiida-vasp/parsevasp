#!/usr/bin/python
import sys
import logging
import re
import numpy as np

from parsevasp import utils
from parsevasp.base import BaseParser

# Map from number of columns in DOSCAR to dtype for the total density of states.
DTYPES_DOS = {
    3: np.dtype([('energy', float), ('total', float), ('integrated', float)]),
    5: np.dtype([('energy', float), ('total', float, (2,)), ('integrated', float, (2,))]),
}

# Map from the number of columns in DDSCAR to dtype for the partial density of states.
DTYPES_PDOS = {
    # l-decomposed
    4:
        np.dtype([('energy', float), ('s', float), ('p', float), ('d', float)]),
    7:
        np.dtype([('energy', float), ('s', float, (2,)), ('p', float, (2,)), ('d', float, (2,))]),
    13:
        np.dtype([('energy', float), ('s', float, (4,)), ('p', float, (4,)), ('d', float, (4,))]),
    5:
        np.dtype([('energy', float), ('s', float), ('p', float), ('d', float), ('f', float)]),
    9:
        np.dtype([('energy', float), ('s', float, (2,)), ('p', float, (2,)), ('d', float, (2,)), ('f', float, (2,))]),
    17:
        np.dtype([('energy', float), ('s', float, (4,)), ('p', float, (4,)), ('d', float, (4,)), ('f', float, (4,))]),
    # lm-decomposed
    10:
        np.dtype([('energy', float), ('s', float), ('py', float), ('px', float), ('pz', float), ('dxy', float), ('dyz', float),
                  ('dz2', float), ('dxz', float), ('dx2-y2', float)]),
    18:
        np.dtype([('energy', float), ('s', float), ('py', float), ('px', float), ('pz', float), ('dxy', float), ('dyz', float),
                  ('dz2', float), ('dxz', float), ('dx2-y2', float), ('fy(3x2-y2)', float), ('fxyz', float), ('fyz2', float),
                  ('fz3', float), ('fxz2', float), ('fz(x2-y2)', float), ('fx(x2-3y2)', float)]),
    19:
        np.dtype([('energy', float), ('s', float, (2,)), ('py', float, (2,)), ('px', float, (2,)), ('pz', float, (2,)),
                  ('dxy', float, (2,)), ('dyz', float, (2,)), ('dz2', float, (2,)), ('dxz', float, (2,)), ('dx2-y2', float, (2,))]),
    35:
        np.dtype([('energy', float), ('s', float, (2,)), ('py', float, (2,)), ('px', float, (2,)), ('pz', float, (2,)),
                  ('dxy', float, (2,)), ('dyz', float, (2,)), ('dz2', float, (2,)), ('dxz', float, (2,)), ('dx2-y2', float, (2,)),
                  ('fy(3x2-y2)', float, (2,)), ('fxyz', float, (2,)), ('fyz2', float, (2,)), ('fz3', float, (2,)), ('fxz2', float, (2,)),
                  ('fz(x2-y2)', float, (2,)), ('fx(x2-3y2)', float, (2,))]),
    37:
        np.dtype([('energy', float), ('s', float, (4,)), ('py', float, (4,)), ('px', float, (4,)), ('pz', float, (4,)),
                  ('dxy', float, (4,)), ('dyz', float, (4,)), ('dz2', float, (4,)), ('dxz', float, (4,)), ('x2-y2', float, (4,))]),
    69:
        np.dtype([('energy', float), ('s', float, (4,)), ('py', float, (4,)), ('px', float, (4,)), ('pz', float, (4,)),
                  ('dxy', float, (4,)), ('dyz', float, (4,)), ('dz2', float, (4,)), ('dxz', float, (4,)), ('dx2-y2', float, (4,)),
                  ('fy(3x2-y2)', float, (4,)), ('fxyz', float, (4,)), ('fyz2', float, (4,)), ('fz3', float, (4,)), ('fxz2', float, (4,)),
                  ('fz(x2-y2)', float, (4,)), ('fx(x2-3y2)', float, (4,))]),
}

# Mapping between the number of columns to the number of spins.
COLSPIN_MAP = {7: 2, 9: 2, 19: 2, 35: 2, 13: 4, 17: 4, 37: 4, 69: 4, 4: 1, 5: 1, 10: 1, 18: 1}


class Doscar(BaseParser):
    
    def __init__(self,
                 file_path=None,
                 file_handler=None,
                 logger=None):
        """
        Initialize an DOSCAR object and set content as a dictionary.

        file_path : string
            A string containing the file path to the file that is going to be parsed.
        file_handler : object
            A file like object that acts as a handler for the content to be parsed.
        logger : object
            A logger object if you would like to use an external logger for messages
            ejected inside this parser.

        """

        super(Doscar, self).__init__(file_path=file_path,
                                     file_handler=file_handler,
                                     logger=logger)

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
            'dos': None,
            'pdos': None,
            'metadata': None
        }

        # parse parse parse
        self._parse()

    def _write(self, *args, **kwargs):
        """Write not supported for DOSCAR."""
        raise NotImplementedError('Writing DOSCAR files is not supported.')

    def _parse(self):
        """Perform the actual parsing by parsing from file like content."""
        self._from_file()

    def _from_file(self):
        """
        Create a dictionary of entries from a
        file and store them in the this instance's data dictionary.

        """

        doscar = utils.read_from_file(self._file_path, self._file_handler)
        self._from_list(doscar)

    def _from_list(self, doscar):
        """
        Go through the list and extract total and partial density of states
        and some metadata.

        Parameters
        ----------
        doscar : list
            A list of strings containing each line in the DOSCAR file.

        """

        # Set some metadata
        num_ions, num_atoms, part, ncdij = utils.line_to_type(doscar[0], int)

        # Figure out if we have a partial density of states
        partial = bool(int(part))
        
        # Volume of cell (AA^3), length of basis vectors (meters) and POTIMS
        line_0 = utils.line_to_type(doscar[1], float)

        # The initial temperature
        line_1 = utils.line_to_type(doscar[2], float)

        # Fetch coordinates used
        coord_type = utils.line_to_type(doscar[3])

        # Name of system
        system = utils.line_to_type(doscar[4], no_split=True)

        # Energy min, energy max, number of points between, fermi level, weight
        line_2 = utils.line_to_type(doscar[5], float)
        emax, emin, ndos, efermi, weight = line_2
        ndos = int(ndos)

        # The rest of the file is density of states data, convert to float
        data = [utils.line_to_type(line, d_type=float) for line in doscar[6:]]
        
        # Get the number of columns for the total dos section to figure out
        # if data is spin decomposed
        count = len(data[ndos - 1])

        num_spin = 1
        if count == 5:
            num_spin = 2

        # Total density of states
        dos_data = np.array(data[:ndos])
        dos = np.zeros((dos_data.shape[0]), DTYPES_DOS[count])
        dos['energy'] = dos_data[:, 0]
        for i, name in enumerate(DTYPES_DOS[count].names[1:]):
            if num_spin == 1:
                dos[name] = np.squeeze(dos_data[:, i + 1:i + 1 + num_spin], axis=1)
            else:
                dos[name] = dos_data[:, i + 1:i + 1 + num_spin]

        # Partial density of states
        pdos_items = []
        pdos = np.array([])  # Partial dos is empty by default
        if line_2 in data:
            for _ in range(num_ions):
                start = data.index(line_2) + 1
                pdos_items += [data[start:start + ndos]]

            # Get the number of columns for the pdos section.
            count = len(pdos_items[-1][-1])
            pdos_data = np.array(pdos_items)

            # Adjust the spin according to the column definitions
            num_spin = COLSPIN_MAP.get(count)
            if num_spin is None:
                raise ValueError(f'Unkown column count: {count} in DOSCAR')

            pdos = np.zeros((pdos_data.shape[0], pdos_data.shape[1]), DTYPES_PDOS[count])
            pdos['energy'] = pdos_data[:, :, 0]
            for i, name in enumerate(DTYPES_PDOS[count].names[1:]):
                if num_spin == 1:  # Only squeeze if there is only one spin component
                    pdos[name] = np.squeeze(pdos_data[:, :, i + 1:i + 1 + num_spin], axis=2)
                else:
                    pdos[name] = pdos_data[:, :, i + 1:i + 1 + num_spin]

        metadata = {}
        metadata['n_ions'] = num_ions
        metadata['n_atoms'] = num_atoms
        metadata['cartesian'] = coord_type.startswith(('c', 'C'))
        metadata['name'] = system
        metadata['emax'] = emax
        metadata['emin'] = emin
        metadata['n_dos'] = ndos
        metadata['efermi'] = efermi
        metadata['weight'] = weight

        # Store
        self._data['metadata'] = metadata
        self._data['pdos'] = pdos
        self._data['dos'] = dos

    def get_metadata(self):
        """
        Return the metadata.

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

    def get_dos(self):
        """
        Return the total density of states.

        Parameters
        ----------
        None

        Returns
        -------
        dos : nparray
            A numpy array containing the total density of states. First index is the
            energy samples, while the last index if composed of the energy sample, total 
            density of states and integrated density of states at that energy sample, respectively.

        """

        dos = self._data['dos']
        return dos

    def get_pdos(self):
        """
        Return the partial density of states.

        Parameters
        ----------
        None

        Returns
        -------
        pdos : nparray

        """

        pdos = self._data['pdos']
        return pdos
