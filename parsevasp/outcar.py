#!/usr/bin/python
import sys
import logging
import numpy as np

from parsevasp import utils
from parsevasp.base import BaseParser


class Outcar(BaseParser):
    def __init__(self,
                 file_path=None,
                 file_handler=None,
                 logger=None,
                 prec=None,
                 conserve_order=False):
        """Initialize an OUTCAR object and set content as a dictionary.

        Parameters
        ----------
        prec : int, optional
            An integer describing how many decimals the users wants
            when printing files.

        """

        super(Outcar, self).__init__(file_path=file_path,
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

        # set precision
        if prec is None:
            self._prec = 12
        else:
            self._prec = prec
        self._width = self._prec + 4

        self._data = {
            'elastic_moduli': {
                'non-symmetrized': None,
                'symmetrized': None,
                'total': None
            },
            'symmetry': {
                'original_cell_type': {
                    'static': [],
                    'dynamic': []
                },
                'symmetrized_cell_type': {
                    'static': [],
                    'dynamic': []
                },
                'num_space_group_operations': {
                    'static': [],
                    'dynamic': []
                },
                'site_symmetry_at_origin': {
                    'static': [],
                    'dynamic': []
                },
                'primitive_translations': [],
                'point_group': {
                    'static': [],
                    'dynamic': []
                }
            },
            'magnetization': {
                'sphere': {
                    'x': {
                        'site_moment': {},
                        'total_magnetization': {}
                    },
                    'y': {
                        'site_moment': {},
                        'total_magnetization': {}
                    },
                    'z': {
                        'site_moment': {},
                        'total_magnetization': {}
                    },
                },
                'full_cell': {},
            }
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

        outcar = utils.readlines_from_file(self._file_path, self._file_handler)
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
        in order to be able to keep the input methods as clean as possible.

        """

        config = ''
        s_orb = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
        for index, line in enumerate(outcar):

            # first, fetch the symmetry
            if line.strip().startswith(
                    'Analysis of symmetry for initial positions (statically)'):
                config = 'static'
            if line.strip().startswith('Analysis of symmetry for dynamics'):
                config = 'dynamic'
            if config:
                if line.strip().startswith('Subroutine PRICEL returns'):
                    text = outcar[index + 1].strip().lower()
                    if text:
                        self._data['symmetry']['original_cell_type'][
                            config].append('primitive cell')
                if ('primitive cells build up your supercell') in line:
                    text = '{} primitive cells'.format(line.strip().split())
                    self._data['symmetry']['original_cell_type'][
                        config].append(text)
                if line.strip().startswith(
                        'Routine SETGRP: Setting up the symmetry group for a'):
                    self._data['symmetry']['symmetrized_cell_type'][
                        config].append(outcar[index + 1].strip().lower())
                if line.strip().startswith('Subroutine GETGRP returns'):
                    self._data['symmetry']['num_space_group_operations'][
                        config].append(int(line.strip().split()[4]))
                if ('configuration has the point symmetry') in line:
                    self._data['symmetry']['site_symmetry_at_origin'][
                        config].append(line.strip().split()[7])
                    next_line = outcar[index + 1].strip().split()
                    if next_line:
                        # Point group only available in recent versions
                        self._data['symmetry']['point_group'][config].append(
                            next_line[10])
                    else:
                        self._data['symmetry']['point_group'][config].append(
                            None)
                    config = ''

            if line.strip().startswith('Subroutine INISYM returns'):
                prim_line = outcar[index + 2].strip().split()
                self._data['symmetry']['primitive_translations'].append(
                    int(prim_line[2]))

            # then the elastic tensors etc. in kBar
            if line.strip().startswith('ELASTIC MODULI  (kBar)'):
                tensor = []
                for idx in range(3, 9):
                    tensor.append([
                        float(item)
                        for item in outcar[index + idx].strip().split()[1:]
                    ])
                self._data['elastic_moduli']['non_symmetrized'] = np.asarray(
                    tensor)
            if line.strip().startswith('SYMMETRIZED ELASTIC MODULI'):
                tensor = []
                for idx in range(3, 9):
                    tensor.append([
                        float(item)
                        for item in outcar[index + idx].strip().split()[1:]
                    ])
                self._data['elastic_moduli']['symmetrized'] = np.asarray(
                    tensor)
            if line.strip().startswith('TOTAL ELASTIC MODULI'):
                tensor = []
                for idx in range(3, 9):
                    tensor.append([
                        float(item)
                        for item in outcar[index + idx].strip().split()[1:]
                    ])
                self._data['elastic_moduli']['total'] = np.asarray(tensor)
            for _proj in ['x', 'y', 'z']:
                if line.strip().startswith('magnetization ({})'.format(_proj)):
                    _counter = 0
                    mag_found = False
                    while not mag_found:
                        if outcar[index + 4 + _counter].strip().split():
                            if not outcar[index + 4 + _counter].strip(
                            ).startswith('-') and not outcar[
                                    index + 4 +
                                    _counter].strip().startswith('tot'):
                                mag_line = outcar[index + 4 + _counter].split()
                                self._data['magnetization']['sphere'][
                                    '{}'.format(_proj)]['site_moment'][int(
                                        mag_line[0])] = dict()
                                for _count, orb in enumerate(mag_line[1:-1]):
                                    self._data['magnetization']['sphere'][
                                        '{}'.format(_proj)]['site_moment'][int(
                                            mag_line[0])][
                                                s_orb[_count]] = float(orb)
                                self._data['magnetization']['sphere'][
                                    '{}'.format(_proj)]['site_moment'][int(
                                        mag_line[0])]['tot'] = float(
                                            mag_line[-1])
                            if outcar[index + 4 +
                                      _counter].strip().startswith('tot'):
                                mag_line = outcar[index + 4 + _counter].split()
                                self._data['magnetization']['sphere'][
                                    '{}'.format(
                                        _proj)]['total_magnetization'] = dict(
                                        )
                                for _count, orb in enumerate(mag_line[1:-1]):
                                    self._data['magnetization']['sphere'][
                                        '{}'.format(
                                            _proj)]['total_magnetization'][
                                                s_orb[_count]] = float(orb)
                                self._data['magnetization']['sphere'][
                                    '{}'.format(_proj)]['total_magnetization'][
                                        'tot'] = float(mag_line[-1])
                                mag_found = True
                        else:
                            self._data['magnetization']['sphere']['{}'.format(
                                _proj)]['total_magnetization'] = dict()
                            self._data['magnetization']['sphere']['{}'.format(_proj)]['total_magnetization'] =\
                                self._data['magnetization']['sphere']['{}'.format(_proj)]['site_moment'][list(self._data['magnetization']['sphere']['{}'.format(_proj)]['site_moment'].keys())[0]]
                            mag_found = True
                        _counter = _counter + 1
            if line.strip().startswith('number of electron'):
                self._data['magnetization']['full_cell'] = [
                    float(_val) for _val in outcar[index].strip().split()[5:]
                ]

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

    def get_magnetization(self):
        """Return the magnetization of the cell.

        Parameters
        ----------
        None

        Returns
        -------
        magnetic : dict
            A dictionary containing the magnetization of the cell.

        """

        magnetic = self._data['magnetization']
        return magnetic
