#!/usr/bin/python
import sys
import logging
import numpy as np
from collections import Counter
import StringIO
from itertools import groupby

import utils


class Outcar(object):

    def __init__(self, file_path=None, logger=None, prec=None, conserve_order=False):
        """Initialize an OUTCAR object and set content as a dictionary.

        Parameters
        ----------
        file_path : string
            The file path in which the OUTCAR is read.
        logger : object, optional
            A standard Python logger object.
        prec : int, optional
            An integer describing how many decimals the users wants
            when printing files.

        """

        self._file_path = file_path
        self._conserve_order = conserve_order

        # check that at least one is suplpied
        if self._file_path is None:
            self._logger.error("Please supply file_path when "
                               "initializing Outcar. Exiting.")
            sys.exit(1)

        # set logger
        if logger is not None:
            self._logger = logger
        else:
            logging.basicConfig(level=logging.DEBUG)
            self._logger = logging.getLogger('OutcarParser')

        # set precision
        if prec is None:
            self._prec = 12
        else:
            self._prec = prec
        self._width = self._prec + 4

        self._data = {'elastic_moduli': {'non-symmetrized': None,
                                         'symmetrized': None,
                                         'total': None},
                      'symmetry': {'original_cell_type': {'static': [], 'dynamic': []},
                                   'symmetrized_cell_type': {'static': [], 'dynamic': []},
                                   'num_space_group_operations': {'static': [], 'dynamic': []},
                                   'num_point_group_operations': {'static': [], 'dynamic': []},
                                   'point_group': {'static': [], 'dynamic': []},
                                   'space_group': {'static': [], 'dynamic': []}}
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

        if self._file_path is None:
            return

        # create dictionary from a file
        self._from_file()

    def _from_file(self):
        """Create a dictionary of entries from a
        file and store them in the this instance's data dictionary.

        """

        outcar = utils.readlines_from_file(self._file_path)
        self._from_list(outcar)

    def _from_list(self, outcar):
        """Go through the list and extract what is not present in the
        XML file.

        Parameters
        ----------
        outcar : list
            A list of strings containing each line in the OUTCAR file.

        Returns
        -------
        outcar_dict : dictionary
            A dictionary containing each OUTCAR tag as a key with the
            associated element.

        Notes
        -----
        No checking for consistency is done here. We do this at a later step
        in order to be able to keep the input methods as clean as posible.

        """


        config = ''
        for index, line in enumerate(outcar):
            
            # first, fetch the symmetry
            if line.strip().startswith('Analysis of symmetry for initial positions (statically)'):
                config = 'static'
            if line.strip().startswith('Analysis of symmetry for dynamics'):
                config = 'dynamic'
            if config:
                if line.strip().startswith('Subroutine PRICEL returns'):
                    text = outcar[index+1].strip().lower()
                    if text:
                        self._data['symmetry']['original_cell_type'][config].append('primitive cell')
                if ('primitive cells build up your supercell') in line :
                    text = '{} primitive cells'.format(line.strip().split())
                    self._data['symmetry']['original_cell_type'][config].append(text)
                if line.strip().startswith('Routine SETGRP: Setting up the symmetry group for a'):
                    self._data['symmetry']['symmetrized_cell_type'][config].append(outcar[index+1].strip().lower())
                if line.strip().startswith('Subroutine GETGRP returns'):
                    self._data['symmetry']['num_space_group_operations'][config].append(int(line.strip().split()[4]))
                    self._data['symmetry']['num_point_group_operations'][config].append(int(outcar[index+1].strip().split()[1]))
                if ('configuration has the point symmetry') in line :
                    self._data['symmetry']['point_group'][config].append(line.strip().split()[7])
                    self._data['symmetry']['space_group'][config].append(outcar[index+1].strip().split()[10])
                    config=''

            # then the elastic tensors etc. in kBar
            if line.strip().startswith('ELASTIC MODULI  (kBar)'):
                tensor = []
                for idx in range(3, 9):
                    tensor.append([float(item) for item in outcar[index+idx].strip().split()[1:]])
                self._data['elastic_moduli']['non_symmetrized'] = np.asarray(tensor)
            if line.strip().startswith('SYMMETRIZED ELASTIC MODULI'):
                tensor = []
                for idx in range(3, 9):
                    tensor.append([float(item) for item in outcar[index+idx].strip().split()[1:]])
                self._data['elastic_moduli']['symmetrized'] = np.asarray(tensor)
            if line.strip().startswith('TOTAL ELASTIC MODULI'):
                tensor = []
                for idx in range(3, 9):
                    tensor.append([float(item) for item in outcar[index+idx].strip().split()[1:]])
                self._data['elastic_moduli']['total'] = np.asarray(tensor)

                
    def get_symmetry(self):
        """Return the symmetry.

        Parameters
        ----------
        None

        Returns
        -------
        symmetry : dict
            A dictionary containing the symmetry information.

        """

        symmetry = self._data['symmetry']
        return symmetry

    def get_elastic_moduli(self):
        """Return the elastic moduli in kBar.

        Parameters
        ----------
        None

        Returns
        -------
        elastic : dict
            A dictionary containing the elastic moduli.

        """

        elastic = self._data['elastic_moduli']
        return elastic
