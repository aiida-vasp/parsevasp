#!/usr/bin/python
import sys
import os
import numpy as np
import logging
import mmap
import copy

from parsevasp import constants
from parsevasp import utils
from parsevasp.base import BaseParser

from lxml import etree

# Try to import lxml, if not present fall back to
# intrinsic ElementTree
lxml = False
try:
    from lxml import etree
    lxml = True
except ImportError:
    try:
        # Python 2.5
        import xml.etree.cElementTree as etree
    except ImportError:
        try:
            # Python 2.5
            import xml.etree.ElementTree as etree
        except ImportError:
            try:
                # normal cElementTree
                import cElementTree as etree
            except ImportError:
                try:
                    # normal ElementTree
                    import elementtree.ElementTree as etree
                except ImportError:
                    logging.error('Failed to import ElementTree.')
                    sys.exit('Failed to import ElementTree.')


_SUPPORTED_TOTAL_ENERGIES = {'energy_extrapolated': 'e_0_energy',
                             'energy_free': 'e_fr_energy',
                             'energy_no_entropy': 'e_wo_entrp'}


class Xml(BaseParser):

    ERROR_MULTIPLE_ENTRIES = 500
    ERROR_NO_SPECIES = 501
    ERROR_NO_ISPIN = 502
    ERROR_NO_NBANDS = 503
    ERROR_NO_KPOINTS = 504
    ERROR_MISMATCH_KPOINTS_NBANDS = 505
    ERROR_UNKNOWN_ELEMENT = 506
    ERROR_UNSUPPORTED_STATUS = 507
    ERROR_NO_SIZE = 508
    ERROR_OVERFLOW = 509
    BaseParser.ERROR_MESSAGES.update({
        ERROR_NO_SPECIES:
        'Please extract the species first.',
        ERROR_MULTIPLE_ENTRIES:
        'Multiple entries of were located.',
        ERROR_NO_ISPIN:
        'Please extract ISPIN first.',
        ERROR_NO_NBANDS:
        'Please extract NBANDS first.',
        ERROR_NO_KPOINTS:
        'Please extract the kpoints first.',
        ERROR_MISMATCH_KPOINTS_NBANDS:
        'The number of kpoints and bands for the entries does not match the number of '
        'located kpoints and NBANDS.',
        ERROR_UNKNOWN_ELEMENT:
        'There is an atomic element present in the XML file that is unknown.',
        ERROR_UNSUPPORTED_STATUS:
        'The supplied status is not supported.',
        ERROR_NO_SIZE:
        'Can not calculate size.',
        ERROR_OVERFLOW:
        'Overflow detected in the XML file.'
    })
    ERROR_MESSAGES = BaseParser.ERROR_MESSAGES

    def __init__(self,
                 file_path=None,
                 file_handler=None,
                 k_before_band=False,
                 extract_all=True,
                 logger=None,
                 event=False):
        """Initialize the XmlParser by first trying the lxml and
        fall back to the standard ElementTree if that is not present.

        Parameters
        ----------
        k_before_band : bool
            If True the kpoint index runs before the bands
            index.
        extract_all : bool
            Extract data from all calculation (i.e. ionic steps)
        event : bool
            If True, force event based method.

        Notes
        -----
        lxml should be used and is required for large files
        """

        super(Xml, self).__init__(file_path=file_path,
                                  file_handler=file_handler,
                                  logger=logger)

        self._sizecutoff = 500
        self._event = event

        if self._file_path is None and self._file_handler is None:
            self._logger.error(
                self.ERROR_MESSAGES[self.ERROR_ONLY_ONE_ARGUMENT])
            sys.exit(self.ERROR_ONLY_ONE_ARGUMENT)

        # extract data from all calculations (e.g. ionic steps)
        self._extract_all = extract_all

        # kpoint index before band index (for instance for the ordering
        # of the eigenvalue data etc.)?
        self._k_before_band = k_before_band

        # version
        self._version = None
        
        # dictionaries that contain the output of the parsing
        self._parameters = {
            'symprec': None,
            'ismear': None,
            'sigma': None,
            'ispin': None,
            'nbands': None,
            'nelect': None,
            'system': None,
            'nelm': None,
            'nsw': None,
        }
        self._lattice = {
            'unitcell': None,
            'species': None,
            'positions': None,
            'kpoints': None,
            'kpointsw': None,
            'kpointdiv': None
        }
        self._data = {
            'eigenvalues': None,
            'eigenvalues_specific': None,
            'eigenvelocities': None,
            'kpoints': None,
            'kpointsw': None,
            'occupancies': None,
            'dos': None,
            'dos_specific': None,
            'totens': None,
            'forces': None,
            'stress': None,
            'dielectrics': None,
            'projectors': None,
            'hessian': None,
            'dynmat': None,
            'born': None
        }

        if lxml:
            self._logger.info('We are utilizing lxml!')
        else:
            self._logger.info('We are not uitilizing lxml!')

        # parse parse parse
        self._parse()

    @property
    def truncated(self):
        """Return True of the xml parsed is truncated."""
        return self._xml_recover
        
    def _parse(self):
        """Perform the actual parsing

        Parameters
        ----------
        None

        Returns
        -------
        None

        """

        # Check size of the XML file. For large files we need to
        # perform event driven parsing. For smaller files this is
        # not necessary and is too slow.
        file_size = self._file_size()
        if file_size is None:
            return None

        # Do a quick check to see if the XML file is not truncated
        self._xml_recover = self._check_xml()
            
        if ((file_size < self._sizecutoff) or self._xml_recover) and \
           not self._event:
            # run regular method (loads file into memory) and
            # enable recovery mode if necessary
            self._parsew(self._xml_recover)
        else:
            # event based, saves a bit of memory
            self._parsee()

    def _parsew(self, xml_recover):
        """Performs parsing on the whole XML files. For smaller files

        """

        self._logger.debug('Running parsew.')

        # now open the complete file, we have already checked its presence when we checked the recover.
        if self._file_handler is None:
            filer = self._file_path
        else:
            filer = self._file_handler

        # make sure we enable the recovery mode
        # pretty sure there is a performance bottleneck running this
        # enabled at all times, so consider to add check for
        # truncated XML files and then enable

        if lxml and xml_recover:
            if xml_recover:
                self._logger.debug('Running LXML in recovery mode.')
            parser = etree.XMLParser(recover=True)
            vaspxml = etree.parse(filer, parser=parser)
        else:
            vaspxml = etree.parse(filer)

        # do we want to extract data from all calculations (e.g. ionic steps)
        extract_all = self._extract_all

        # let us start to parse the content
        self._version = self._fetch_versionw(vaspxml)
        self._parameters['symprec'] = self._fetch_symprecw(vaspxml)
        self._parameters['sigma'] = self._fetch_sigmaw(vaspxml)
        self._parameters['ismear'] = self._fetch_ismearw(vaspxml)
        self._parameters['ispin'] = self._fetch_ispinw(vaspxml)
        self._parameters['nbands'] = self._fetch_nbandsw(vaspxml)
        self._parameters['nelect'] = self._fetch_nelectw(vaspxml)
        self._parameters['system'] = self._fetch_systemw(vaspxml)
        self._parameters['nelm'] = self._fetch_nelmw(vaspxml)
        self._parameters['nsw'] = self._fetch_nsww(vaspxml)
        self._lattice['species'] = self._fetch_speciesw(vaspxml)
        self._lattice['unitcell'], self._lattice['positions'], \
            self._data['forces'], self._data['stress'] = \
            self._fetch_upfsw(vaspxml, extract_all=extract_all)
        self._lattice['kpoints'] = self._fetch_kpointsw(vaspxml)
        self._lattice['kpointsw'] = self._fetch_kpointsww(vaspxml)
        self._lattice['kpointdiv'] = self._fetch_kpointdivw(vaspxml)
        self._data['eigenvalues'], self._data[
            'occupancies'] = self._fetch_eigenvaluesw(vaspxml)
        self._data['eigenvalues_specific'] = self._fetch_eigenvalues_specificw(
            vaspxml)
        self._data['eigenvelocities'] = self._fetch_eigenvelocitiesw(vaspxml)
        self._data['dos'], self._data['dos_specific'] = self._fetch_dosw(
            vaspxml)
        self._data['totens'] = self._fetch_totensw(vaspxml)
        self._data['dielectrics'] = self._fetch_dielectricsw(vaspxml)
        self._data['projectors'] = self._fetch_projectorsw(vaspxml)
        self._data['hessian'] = self._fetch_hessian(vaspxml)
        self._data['dynmat'] = self._fetch_dynmatw(vaspxml)
        self._data['born'] = self._fetch_bornw(vaspxml)

    def _parsee(self):
        """Performs parsing in an event driven fashion on the XML file.
        Slower, but suitable for bigger files.

        """

        # set logger
        self._logger.debug('Running parsee.')

        # helper lists
        data = []
        data2 = []
        data3 = []
        data4 = []
        data5 = []
        data6 = []

        # dicts
        cell = {}
        pos = {}
        force = {}
        stress = {}
        dos = {}
        totens = {}
        dynmat = {}
        _dos = {}
        _dos2 = {}

        # bool to control extraction of content
        extract_generator = False
        extract_parameters = False
        extract_calculation = False
        extract_latticedata = False
        extract_unitcell = False
        extract_positions = False
        extract_species = False
        extract_kpointdata = False
        extract_kpoints = False
        extract_kpointsw = False
        extract_kpointdiv = False
        extract_kpoints_specific = False
        extract_kpointsw_specific = False
        extract_eigenvalues = False
        extract_eigenvalues_specific = False
        extract_eigenvalues_spin1 = False
        extract_eigenvalues_spin2 = False
        extract_eigenvalues_specific_spin1 = False
        extract_eigenvalues_specific_spin2 = False
        extract_eigenvelocities = False
        extract_eigenvelocities_spin1 = False
        extract_eigenvelocities_spin2 = False
        extract_dos = False
        extract_dos_specific = False
        extract_total_dos = False
        extract_partial_dos = False
        extract_dos_ispin1 = False
        extract_dos_ispin2 = False
        extract_dos_specific_ispin1 = False
        extract_dos_specific_ispin2 = False
        extract_projected = False
        extract_forces = False
        extract_stress = False
        extract_energies = False
        extract_e_0_energy = False
        extract_e_fr_energy = False
        extract_e_wo_entrp = False
        extract_scstep = False
        extract_dielectrics = False
        extract_eig_proj = False
        extract_eig_proj_ispin1 = False
        extract_eig_proj_ispin2 = False
        extract_dynmat = False
        extract_dynmat_eigen = False
        extract_hessian = False
        extract_born = False

        # do we want to extract data from all calculations (e.g. ionic steps)
        extract_all = self._extract_all

        if self._file_handler is None:
            filer = self._file_path
        else:
            filer = self._file_handler

        # index that control the calculation step (e.g. ionic step)
        calc = 1
        for event, element in etree.iterparse(filer, events=('start', 'end')):
            # set extraction points (what to read and when to read it)
            # here we also set the relevant data elements when the tags
            # close when they contain more than one element
            if event == 'start' and element.tag == 'generator':
                extract_generator = True
            if event == 'end' and element.tag == 'generator':
                extract_generator = False
            if event == 'start' and element.tag == 'parameters':
                extract_parameters = True
            if event == 'end' and element.tag == 'parameters':
                extract_parameters = False
            if event == 'start' and element.tag == 'calculation':
                # Instead of needing to check the nested dicts we initialize them here
                # so that we can use update for each calculation.
                totens[calc] = {}
                extract_calculation = True
            if event == 'end' and element.tag == 'calculation':
                data3 = self._convert_array1D_f(data3)
                data4 = self._convert_array1D_f(data4)
                data5 = self._convert_array1D_f(data5)
                totens[calc].update({'energy_extrapolated': data3,
                                     'energy_free': data4,
                                     'energy_no_entropy': data5
                })
                data3 = []
                data4 = []
                data5 = []
                # update index for the calculation
                calc = calc + 1
                extract_calculation = False
            if event == 'start' and element.tag == 'array' \
               and element.attrib.get('name') == 'atoms':
                extract_species = True
            if event == 'end' and element.tag == 'array' \
               and element.attrib.get('name') == 'atoms':
                # only need every other element (element, not atomtype)
                self._lattice['species'] = self._convert_species(data[::2])
                data = []
                extract_species = False
            if event == 'start' and element.tag == 'kpoints' and not extract_calculation:
                extract_kpointdata = True
            if event == 'end' and element.tag == 'kpoints' and not extract_calculation:
                extract_kpointdata = False
            if event == 'start' and element.tag == 'projected':
                extract_projected = True
            if event == 'end' and element.tag == 'projected':
                extract_projected = False

            # now fetch the data
            if extract_generator:
                try:
                    if event == 'start' and element.attrib['name'] == 'version':
                        self._version = element.text
                except KeyError:
                    pass
            if extract_parameters:
                try:
                    if event == 'start' and element.attrib['name'] == 'SYMPREC':
                        self._parameters['symprec'] = self._convert_f(element)
                except KeyError:
                    pass
                try:
                    if event == 'start' and element.attrib['name'] == 'ISPIN':
                        self._parameters['ispin'] = self._convert_i(element)
                except KeyError:
                    pass
                try:
                    if event == 'start' and element.attrib['name'] == 'ISMEAR':
                        self._parameters['ismear'] = self._convert_i(element)
                except KeyError:
                    pass
                try:
                    if event == 'start' and element.attrib['name'] == 'SIGMA':
                        self._parameters['sigma'] = self._convert_f(element)
                except KeyError:
                    pass
                try:
                    if event == 'start' and element.attrib['name'] == 'NBANDS':
                        self._parameters['nbands'] = self._convert_i(element)
                except KeyError:
                    pass
                try:
                    if event == 'start' and element.attrib['name'] == 'NELECT':
                        self._parameters['nelect'] = self._convert_f(element)
                except KeyError:
                    pass
                try:
                    if event == 'start' and element.attrib['name'] == 'SYSTEM':
                        self._parameters['system'] = element.text
                except KeyError:
                    pass
                try:
                    if event == 'start' and element.attrib['name'] == 'NELM' \
                        and element.getparent().attrib['name'] == 'electronic convergence':
                        self._parameters['nelm'] = self._convert_i(element)
                except KeyError:
                    pass
                try:
                    if event == 'start' and element.attrib['name'] == 'NSW':
                        self._parameters['nsw'] = self._convert_i(element)
                except KeyError:
                    pass

            if extract_calculation:
                # it would be very tempting just to fill the data and disect
                # it later, would be faster, but it is not so easy since
                # we do not know how many calculations have been performed
                # or how many scteps there are per calculation
                if event == 'start' and element.tag == 'dos' and element.attrib.get(
                        'comment') is None:
                    extract_dos = True
                if event == 'end' and element.tag == 'total' and extract_dos:
                    if data2:
                        # only store energy for one part as
                        # this is the same for both
                        dos_ispin = self._convert_array2D_f(data, 3)
                        _dos['energy'] = dos_ispin[:, 0]
                        _dos['total'] = dos_ispin[:, 1]
                        _dos['integrated'] = dos_ispin[:, 2]
                        _dos['partial'] = None
                        dos_ispin = self._convert_array2D_f(data2, 3)
                        _dos2['total'] = dos_ispin[:, 1]
                        _dos2['integrated'] = dos_ispin[:, 2]
                        _dos2['partial'] = None
                    else:
                        dos_ispin = self._convert_array2D_f(data, 3)
                        _dos['energy'] = dos_ispin[:, 0]
                        _dos['total'] = dos_ispin[:, 1]
                        _dos['integrated'] = dos_ispin[:, 2]
                        _dos['partial'] = None
                    data = []
                    data2 = []
                if event == 'end' and element.tag == 'partial' and extract_dos:
                    num_atoms = 0
                    if self._lattice['species'] is not None:
                        num_atoms = self._lattice['species'].shape[0]
                    else:
                        self._logger.error(
                            self.ERROR_MESSAGES[self.ERROR_NO_SPECIES])
                        sys.exit(self.ERROR_NO_SPECIES)
                    if data2:
                        dos_ispin = self._convert_array2D_f(data, 10)
                        # do not need the energy term (similar to total)
                        _dos['partial'] = np.asarray(
                            np.split(dos_ispin[:, 1:10], num_atoms))
                        dos_ispin = self._convert_array2D_f(data2, 10)
                        # do not need the energy term (similar to total)
                        _dos2['partial'] = np.asarray(
                            np.split(dos_ispin[:, 1:10], num_atoms))
                    else:
                        dos_ispin = self._convert_array2D_f(data, 10)
                        # do not need the energy term (similar to total)
                        _dos['partial'] = np.asarray(
                            np.split(dos_ispin[:, 1:10], num_atoms))
                    data = []
                    data2 = []
                if event == 'end' and element.tag == 'dos' and extract_dos:
                    # check the Fermi level
                    if len(data6) == 1:
                        fermi_level = self._convert_f(data6[0])
                    elif len(data6) > 1:
                        self._logger.error(
                            self.ERROR_MESSAGES[self.ERROR_MULTIPLE_ENTRIES] +
                            " The tag in question is 'efermi'.")
                        sys.exit(self.ERROR_MULTIPLE_ENTRIES)
                    else:
                        fermi_level = None

                    if _dos2:
                        dos['up'] = _dos
                        dos['down'] = _dos2
                        dos['total'] = {
                            'fermi_level': fermi_level,
                            'energy': _dos['energy']
                        }
                        del dos['up']['energy']
                    else:
                        _dos['fermi_level'] = fermi_level
                        dos['total'] = _dos
                    self._data['dos'] = copy.deepcopy(dos)
                    data = []
                    data2 = []
                    data6 = []
                    _dos = {}
                    _dos2 = {}
                    extract_dos = False
                if event == 'start' and element.tag == 'dos' and element.attrib.get(
                        'comment') == 'interpolated':
                    extract_dos_specific = True
                if event == 'end' and element.tag == 'total' and extract_dos_specific:
                    if data2:
                        # only store energy for one part as
                        # this is the same for both
                        dos_ispin = self._convert_array2D_f(data, 3)
                        _dos['energy'] = dos_ispin[:, 0]
                        _dos['total'] = dos_ispin[:, 1]
                        _dos['integrated'] = dos_ispin[:, 2]
                        _dos['partial'] = None
                        dos_ispin = self._convert_array2D_f(data2, 3)
                        _dos2['total'] = dos_ispin[:, 1]
                        _dos2['integrated'] = dos_ispin[:, 2]
                        _dos2['partial'] = None
                    else:
                        dos_ispin = self._convert_array2D_f(data, 3)
                        _dos['energy'] = dos_ispin[:, 0]
                        _dos['total'] = dos_ispin[:, 1]
                        _dos['integrated'] = dos_ispin[:, 2]
                        _dos['partial'] = None
                    data = []
                    data2 = []
                if event == 'end' and element.tag == 'partial' and extract_dos_specific:
                    num_atoms = 0
                    if self._lattice['species'] is not None:
                        num_atoms = self._lattice['species'].shape[0]
                    else:
                        self._logger.error(
                            self.ERROR_MESSAGES[self.ERROR_NO_SPECIES])
                        sys.exit(self.ERROR_NO_SPECIES)
                    if data2:
                        dos_ispin = self._convert_array2D_f(data, 10)
                        # do not need the energy term (similar to total)
                        _dos['partial'] = np.asarray(
                            np.split(dos_ispin[:, 1:10], num_atoms))
                        dos_ispin = self._convert_array2D_f(data2, 10)
                        # do not need the energy term (similar to total)
                        _dos2['partial'] = np.asarray(
                            np.split(dos_ispin[:, 1:10], num_atoms))
                    else:
                        dos_ispin = self._convert_array2D_f(data, 10)
                        # do not need the energy term (similar to total)
                        _dos['partial'] = np.asarray(
                            np.split(dos_ispin[:, 1:10], num_atoms))
                    data = []
                    data2 = []
                if event == 'end' and element.tag == 'dos' and extract_dos_specific:
                    # check the Fermi level
                    if len(data6) == 1:
                        fermi_level = self._convert_f(data6[0])
                    elif len(data6) > 1:
                        self._logger.error(
                            self.ERROR_MESSAGES[self.ERROR_MULTIPLE_ENTRIES] +
                            " The tag in question is 'efermi'.")
                        sys.exit(self.ERROR_MULTIPLE_ENTRIES)
                    else:
                        fermi_level = None

                    if _dos2:
                        dos['up'] = _dos
                        dos['down'] = _dos2
                        dos['total'] = {
                            'fermi_level': fermi_level,
                            'energy': _dos['energy']
                        }
                        del dos['up']['energy']
                    else:
                        _dos['fermi_level'] = fermi_level
                        dos['total'] = _dos
                    self._data['dos_specific'] = dos
                    data = []
                    data2 = []
                    data6 = []
                    _dos = {}
                    _dos2 = {}
                    extract_dos_specific = False

                if event == 'start' and element.tag == 'structure':
                    extract_latticedata = True
                if event == 'end' and element.tag == 'structure':
                    extract_latticedata = False
                if event == 'start' and element.tag == 'varray' and \
                   element.attrib['name'] == 'forces':
                    extract_forces = True
                if event == 'end' and element.tag == 'varray' and \
                   element.attrib['name'] == 'forces':
                    force[calc] = self._convert_array2D_f(data, 3)
                    data = []
                    extract_forces = False
                if event == 'start' and element.tag == 'varray' and \
                   element.attrib['name'] == 'stress':
                    extract_stress = True
                if event == 'end' and element.tag == 'varray' and \
                   element.attrib['name'] == 'stress':
                    stress[calc] = self._convert_array2D_f(data, 3)
                    data = []
                    extract_stress = False
                if event == 'start' and element.tag == 'energy' and not extract_scstep:
                    extract_energies = True
                if event == 'end' and element.tag == 'energy' and not extract_scstep:
                    extract_energies = False
                if event == 'start' and element.tag == 'scstep':
                    extract_scstep = True
                if event == 'end' and element.tag == 'scstep':
                    extract_scstep = False
                if event == 'start' and element.tag == 'eigenvalues' and not extract_eigenvelocities and element.attrib.get(
                        'comment'
                ) != 'interpolated' and not extract_eigenvalues_specific:
                    extract_eigenvalues = True
                if event == 'end' and element.tag == 'eigenvalues' and extract_eigenvalues:
                    num_kpoints = len(self._lattice['kpoints'])
                    if not data2:
                        eigenvalues, occupancies = self._extract_eigenvalues(
                            data, None, num_kpoints)
                    else:
                        eigenvalues, occupancies = self._extract_eigenvalues(
                            data, data2, num_kpoints)
                    self._data['eigenvalues'] = eigenvalues
                    self._data['occupancies'] = occupancies
                    data = []
                    data2 = []
                    extract_eigenvalues = False
                if event == 'start' and element.tag == 'eigenvalues' and element.attrib.get(
                        'comment') == 'interpolated':
                    extract_eigenvalues_specific = True

                if event == 'end' and element.tag == 'eigenvalues' and extract_eigenvalues_specific:
                    num_kpoints = len(self._data['kpoints'])
                    if not data2:
                        eigenvalues_specific, _ = self._extract_eigenvalues(
                            data, None, num_kpoints)
                    else:
                        eigenvalues_specific, _ = self._extract_eigenvalues(
                            data, data2, num_kpoints)
                    self._data['eigenvalues_specific'] = eigenvalues_specific
                    data = []
                    data2 = []
                    extract_eigenvalues_specific = False
                if event == 'start' and element.tag == 'eigenvelocities':
                    extract_eigenvelocities = True
                if event == 'end' and element.tag == 'eigenvelocities':
                    num_kpoints = len(self._data['kpoints'])
                    if not data2:
                        eigenvelocities = self._extract_eigenvelocities(
                            data, None, num_kpoints)
                    else:
                        eigenvelocities = self._extract_eigenvelocities(
                            data, data2, num_kpoints)
                    self._data['eigenvelocities'] = eigenvelocities
                    data = []
                    data2 = []
                    extract_eigenvelocities = False
                if event == 'start' and element.tag == 'dielectricfunction':
                    extract_dielectrics = True
                if event == 'end' and element.tag == 'dielectricfunction':
                    _diel = {}
                    diel = np.split(self._convert_array2D_f(data, 7), 2)
                    _diel['energy'] = diel[0][:, 0]
                    _diel['imag'] = diel[0][:, 1:7]
                    _diel['real'] = diel[1][:, 1:7]
                    self._data['dielectrics'] = _diel
                    data = []
                    extract_dielectrics = False
                if event == 'start' and element.tag == 'dynmat':
                    extract_dynmat = True
                if event == 'end' and element.tag == 'dynmat':
                    self._data['dynmat'] = dynmat
                    extract_dynmat = False
                if event == 'start' and element.tag == 'array':
                    # a bit of special threatment here as there is
                    # an array element without attributes, so we get
                    # KeyErrors
                    try:
                        if element.attrib['name'] == 'born_charges':
                            extract_born = True
                    except KeyError:
                        pass
                if event == 'end' and element.tag == 'array':
                    # again a bit special
                    try:
                        if element.attrib['name'] == 'born_charges':
                            num_atoms = 0
                            if self._lattice['species'] is not None:
                                num_atoms = self._lattice['species'].shape[0]
                            else:
                                self._logger.error(
                                    self.ERROR_MESSAGES[self.ERROR_NO_SPECIES])
                                sys.exit(self.ERROR_NO_SPECIES)
                            data = self._convert_array2D_f(data, 3)
                            data = np.split(data, num_atoms)
                            self._data['born'] = np.asarray(data)
                            data = []
                            extract_born = False
                    except KeyError:
                        pass

                # now extract data
                if extract_scstep:
                    # extrapolated energy
                    if event == 'start' and element.tag == 'i' and \
                       element.attrib['name'] == 'e_0_energy':
                        extract_e_0_energy = True
                    if event == 'end' and element.tag == 'i' and \
                       element.attrib['name'] == 'e_0_energy':
                        extract_e_0_energy = False
                    if extract_e_0_energy:
                        data3.append(element)
                    # free energy
                    if event == 'start' and element.tag == 'i' and \
                       element.attrib['name'] == 'e_fr_energy':
                        extract_e_fr_energy = True
                    if event == 'end' and element.tag == 'i' and \
                       element.attrib['name'] == 'e_fr_energy':
                        extract_e_fr_energy = False
                    if extract_e_fr_energy:
                        data4.append(element)
                    # energy without entropy
                    if event == 'start' and element.tag == 'i' and \
                       element.attrib['name'] == 'e_wo_entrp':
                        extract_e_wo_entrp = True
                    if event == 'end' and element.tag == 'i' and \
                       element.attrib['name'] == 'e_wo_entrp':
                        extract_e_wo_entrp = False
                    if extract_e_wo_entrp:
                        data5.append(element)
                if extract_latticedata:
                    if event == 'start' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'basis':
                        extract_unitcell = True
                    if event == 'end' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'basis':
                        cell[calc] = self._convert_array2D_f(data, 3)
                        data = []
                        extract_unitcell = False

                    if event == 'start' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'positions':
                        extract_positions = True
                    if event == 'end' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'positions':
                        pos[calc] = self._convert_array2D_f(data, 3)
                        data = []
                        extract_positions = False

                    if extract_unitcell:
                        if event == 'start' and element.tag == 'v':
                            data.append(element)
                    if extract_positions:
                        if event == 'start' and element.tag == 'v':
                            data.append(element)

                if extract_forces:
                    if event == 'start' and element.tag == 'v':
                        data.append(element)

                if extract_stress:
                    if event == 'start' and element.tag == 'v':
                        data.append(element)

                if extract_energies:
                    # extrapolated energy
                    if event == 'start' and element.tag == 'i' and \
                       element.attrib['name'] == 'e_0_energy':
                        totens[calc].update({'energy_extrapolated_final': float(element.text)})
                    # free energy
                    if event == 'start' and element.tag == 'i' and \
                       element.attrib['name'] == 'e_fr_energy':
                        totens[calc].update({'energy_free_final': float(element.text)})
                    # energy without entropy
                    if event == 'start' and element.tag == 'i' and \
                       element.attrib['name'] == 'e_wo_entrp':
                        totens[calc].update({'energy_no_entropy_final': float(element.text)})

                if extract_eigenvalues:
                    if event == 'start' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 1':
                        extract_eigenvalues_spin1 = True
                    if event == 'end' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 1':
                        extract_eigenvalues_spin1 = False
                    if event == 'start' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 2':
                        extract_eigenvalues_spin2 = True
                    if event == 'end' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 2':
                        extract_eigenvalues_spin2 = False
                    if extract_eigenvalues_spin1:
                        if event == 'start' and element.tag == 'r':
                            data.append(element)
                    if extract_eigenvalues_spin2:
                        if event == 'start' and element.tag == 'r':
                            data2.append(element)

                if extract_eigenvalues_specific:
                    if event == 'start' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'kpointlist':
                        extract_kpoints_specific = True
                    if event == 'end' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'kpointlist':
                        self._data['kpoints'] = self._convert_array2D_f(
                            data, 3)
                        data = []
                        extract_kpoints_specific = False
                    if event == 'start' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'weights':
                        extract_kpointsw_specific = True
                    if event == 'end' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'weights':
                        self._data['kpointsw'] = self._convert_array1D_f(data)
                        data = []
                        extract_kpointsw_specific = False
                    if event == 'start' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 1':
                        extract_eigenvalues_spin1 = True
                    if event == 'end' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 1':
                        extract_eigenvalues_spin1 = False
                    if event == 'start' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 2':
                        extract_eigenvalues_spin2 = True
                    if event == 'end' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 2':
                        extract_eigenvalues_spin2 = False
                    if extract_kpoints_specific:
                        if event == 'start' and element.tag == 'v':
                            data.append(element)
                    if extract_kpointsw_specific:
                        if event == 'start' and element.tag == 'v':
                            data.append(element)
                    if extract_eigenvalues_spin1:
                        if event == 'start' and element.tag == 'r':
                            data.append(element)
                    if extract_eigenvalues_spin2:
                        if event == 'start' and element.tag == 'r':
                            data2.append(element)

                if extract_eigenvelocities:
                    if event == 'start' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'kpointlist':
                        extract_kpoints_specific = True
                    if event == 'end' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'kpointlist':
                        self._data['kpoints'] = self._convert_array2D_f(
                            data, 3)
                        data = []
                        extract_kpoints_specific = False
                    if event == 'start' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'weights':
                        extract_kpointsw_specific = True
                    if event == 'end' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'weights':
                        self._data['kpointsw'] = self._convert_array1D_f(data)
                        data = []
                        extract_kpointsw_specific = False
                    if event == 'start' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 1':
                        extract_eigenvelocities_spin1 = True
                    if event == 'end' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 1':
                        extract_eigenvelocities_spin1 = False
                    if event == 'start' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 2':
                        extract_eigenvelocities_spin2 = True
                    if event == 'end' and element.tag == 'set' \
                       and element.attrib.get('comment') == 'spin 2':
                        extract_eigenvelocities_spin2 = False
                    if extract_kpoints_specific:
                        if event == 'start' and element.tag == 'v':
                            data.append(element)
                    if extract_kpointsw_specific:
                        if event == 'start' and element.tag == 'v':
                            data.append(element)
                    if extract_eigenvelocities_spin1:
                        if event == 'start' and element.tag == 'r':
                            data.append(element)
                    if extract_eigenvelocities_spin2:
                        if event == 'start' and element.tag == 'r':
                            data2.append(element)

                if extract_dielectrics:
                    if event == 'start' and element.tag == 'r':
                        data.append(element)

                if extract_projected:
                    # make sure we skip the first entry containing
                    # the eigenvalues (already stored at this point)
                    if event == 'end' and element.tag == 'eigenvalues':
                        extract_eig_proj = True
                    if event == 'end' and element.tag == 'array' and \
                       extract_eig_proj:
                        if not data2:
                            projectors = self._extract_projectors(data, None)
                        else:
                            projectors = self._extract_projectors(data, data2)
                        self._data['projectors'] = projectors
                        data = []
                        data2 = []
                        extract_eig_proj = False

                    if extract_eig_proj:
                        if event == 'start' and element.tag == 'set' \
                           and element.attrib.get('comment') == 'spin1':
                            extract_eig_proj_ispin1 = True
                        if event == 'end' and element.tag == 'set' \
                           and element.attrib.get('comment') == 'spin1':
                            extract_eig_proj_ispin1 = False
                        if event == 'start' and element.tag == 'set' \
                           and element.attrib.get('comment') == 'spin2':
                            extract_eig_proj_ispin2 = True
                        if event == 'end' and element.tag == 'set' \
                           and element.attrib.get('comment') == 'spin2':
                            extract_eig_proj_ispin2 = False
                        if extract_eig_proj_ispin1:
                            if event == 'start' and element.tag == 'r':
                                data.append(element)
                        if extract_eig_proj_ispin2:
                            if event == 'start' and element.tag == 'r':
                                data2.append(element)

                if extract_dynmat:
                    if event == 'start' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'hessian':
                        extract_hessian = True
                    if event == 'end' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'hessian':
                        num_atoms = 0
                        if self._lattice['species'] is not None:
                            num_atoms = self._lattice['species'].shape[0]
                        else:
                            self._logger.error(
                                self.ERROR_MESSAGES[self.ERROR_NO_SPECIES])
                            sys.exit(self.ERROR_NO_SPECIES)
                        hessian = self._convert_array2D_f(data, num_atoms * 3)
                        self._data['hessian'] = hessian
                        data = []
                        extract_hessian = False
                    if event == 'start' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'eigenvectors':
                        extract_dynmat_eigen = True
                    if event == 'end' and element.tag == 'varray' \
                       and element.attrib.get('name') == 'eigenvectors':
                        num_atoms = 0
                        if self._lattice['species'] is not None:
                            num_atoms = self._lattice['species'].shape[0]
                        else:
                            self._logger.error(
                                self.ERROR_MESSAGES[self.ERROR_NO_SPECIES])
                            sys.exit(self.ERROR_NO_SPECIES)
                        eigenvec = self._convert_array2D_f(data, num_atoms * 3)
                        dynmat['eigenvectors'] = eigenvec
                        data = []
                        extract_dynmat_eigen = False
                    if extract_hessian:
                        if event == 'start' and element.tag == 'v':
                            data.append(element)
                    if extract_dynmat_eigen:
                        if event == 'start' and element.tag == 'v':
                            data.append(element)
                    try:
                        if event == 'start' and \
                           element.attrib['name'] == 'eigenvalues':
                            dynmat['eigenvalues'] = self._convert_array_f(
                                element)
                    except KeyError:
                        pass

                if extract_born:
                    if event == 'start' and element.tag == 'v':
                        data.append(element)

            if extract_species:
                if event == 'start' and element.tag == 'c':
                    data.append(element)

            if extract_kpointdata:
                try:
                    if event == 'start' and element.tag == 'v' and \
                       element.attrib['name'] == 'divisions':
                        self._lattice['kpointdiv'] = self._convert_array_i(
                            element)
                except KeyError:
                    pass

                if event == 'start' and element.tag == 'varray' \
                   and element.attrib.get('name') == 'kpointlist':
                    extract_kpoints = True
                if event == 'end' and element.tag == 'varray' \
                   and element.attrib.get('name') == 'kpointlist':
                    self._lattice['kpoints'] = self._convert_array2D_f(data, 3)
                    data = []
                    extract_kpoints = False
                if event == 'start' and element.tag == 'varray' \
                   and element.attrib.get('name') == 'weights':
                    extract_kpointsw = True
                if event == 'end' and element.tag == 'varray' \
                   and element.attrib.get('name') == 'weights':
                    self._lattice['kpointsw'] = self._convert_array1D_f(data)
                    data = []
                    extract_kpointsw = False
                if extract_kpoints:
                    if event == 'start' and element.tag == 'v':
                        data.append(element)
                if extract_kpointsw:
                    if event == 'start' and element.tag == 'v':
                        data.append(element)

            if extract_dos:
                if event == 'start' and element.tag == 'i' and \
                   element.attrib.get('name') == 'efermi':
                    data6.append(element)
                if event == 'start' and element.tag == 'set' \
                   and element.attrib.get('comment') == 'spin 1':
                    extract_dos_ispin1 = True
                if event == 'end' and element.tag == 'set' \
                   and element.attrib.get('comment') == 'spin 1':
                    extract_dos_ispin1 = False
                if event == 'start' and element.tag == 'set' \
                   and element.attrib.get('comment') == 'spin 2':
                    extract_dos_ispin2 = True
                if event == 'end' and element.tag == 'set' \
                   and element.attrib.get('comment') == 'spin 2':
                    extract_dos_ispin2 = False
                if extract_dos_ispin1:
                    if event == 'start' and element.tag == 'r':
                        data.append(element)
                if extract_dos_ispin2:
                    if event == 'start' and element.tag == 'r':
                        data2.append(element)

            if extract_dos_specific:
                if event == 'start' and element.tag == 'i' and \
                   element.attrib.get('name') == 'efermi':
                    data6.append(element)
                if event == 'start' and element.tag == 'set' \
                   and element.attrib.get('comment') == 'spin 1':
                    extract_dos_specific_ispin1 = True
                if event == 'end' and element.tag == 'set' \
                   and element.attrib.get('comment') == 'spin 1':
                    extract_dos_specific_ispin1 = False
                if event == 'start' and element.tag == 'set' \
                   and element.attrib.get('comment') == 'spin 2':
                    extract_dos_specific_ispin2 = True
                if event == 'end' and element.tag == 'set' \
                   and element.attrib.get('comment') == 'spin 2':
                    extract_dos_specific_ispin2 = False
                if extract_dos_specific_ispin1:
                    if event == 'start' and element.tag == 'r':
                        data.append(element)
                if extract_dos_specific_ispin2:
                    if event == 'start' and element.tag == 'r':
                        data2.append(element)

        # now we need to update some elements
        # for static runs, initial is equal to last
        # if cell:
        #     if len(cell) == 1:
        #         cell[2] = cell[1]
        # if pos:
        #     if len(pos) == 1:
        #         pos[2] = pos[1]
        # if force:
        #     if len(force) == 1:
        #         force[2] = force[1]
        # if stress:
        #     if len(stress) == 1:
        #         stress[2] = stress[1]
        # if totens:
        #     if len(totens) == 1:
        #         totens[2] = totens[1]

        # if not extract_all:
        #     # only save initial and last
        #     if cell:
        #         cell = {key: np.asarray(cell[key]) for key in {1, 2}}
        #     if pos:
        #         pos = {key: np.asarray(pos[key]) for key in {1, 2}}
        #     if force:
        #         force = {key: force[key] for key in {1, 2}}
        #     if stress:
        #         stress = {key: stress[key] for key in {1, 2}}
        #     if totens:
        #         totens = {key: totens[key] for key in {1, 2}}

        # if any dict is empty, set to zero
        if not cell:
            cell = None
        if not pos:
            pos = None
        if not force:
            force = None
        if not stress:
            stress = None
        if not totens:
            totens = None

        # store
        self._lattice['unitcell'] = cell
        self._lattice['positions'] = pos
        self._data['forces'] = force
        self._data['stress'] = stress
        self._data['totens'] = totens

        return

    def _fetch_versionw(self, xml):
        """Fetch and set version using etree

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        version : string
            If version is found it is returned.

        Notes
        -----
        Used when detecting the VASP version.

        """

        entry = self._find(
            xml, './/generator/i[@name="version"]')

        if entry is None:
            return None

        return entry.text
    
    def _fetch_symprecw(self, xml):
        """Fetch and set symprec using etree

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        symprec : float
            If SYMPREC is found it is returned.

        Notes
        -----
        Used when detecting symmetry.

        """

        entry = self._find(
            xml, './/parameters/separator[@name="symmetry"]/'
            'i[@name="SYMPREC"]')

        if entry is None:
            return None

        symprec = self._convert_f(entry)

        return symprec

    def _fetch_sigmaw(self, xml):
        """Fetch and set sigma using etree

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        sigma : float
            If SIGMA is found it is returned.

        Notes
        -----
        Determines the smearing used etc.

        """

        entry = self._find(
            xml, './/parameters/separator[@name="electronic"]/'
            'separator[@name="electronic smearing"]/'
            'i[@name="SIGMA"]')

        if entry is None:
            return None

        sigma = self._convert_f(entry)

        return sigma

    def _fetch_nelmw(self, xml):
        """Fetch and set nelm using etree

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        nelm : float
            If NELM is found it is returned.

        Notes
        -----
        Maximum number of eletronic steps. This is needed for checking if the eletronic
        structure is converged. 

        """

        entry = self._find(
            xml, './/parameters/separator[@name="electronic"]/'
            'separator[@name="electronic convergence"]/'
            'i[@name="NELM"]')

        if entry is None:
            return None

        nelm = self._convert_i(entry)

        return nelm

    def _fetch_nsww(self, xml):
        """Fetch and set nsw using etree

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        nsw : int
            If NSW is found it is returned.

        Notes
        -----
        Maximum number of eletronic steps. This is needed for checking ionic convergence. 
        """

        entry = self._find(
            xml, './/parameters/separator[@name="ionic"]/'
            'i[@name="NSW"]')

        if entry is None:
            return None

        nsw = self._convert_i(entry)

        return nsw

    def _fetch_ispinw(self, xml):
        """Fetch and set ispin using etree

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        ispin : int
            If ISPIN is found it is returned.

        Notes
        -----
        Determines if spin is included. ISPIN=2 separates the spins.

        """

        entry = self._find(
            xml, './/parameters/separator[@name="electronic"]/'
            'separator[@name="electronic spin"]/'
            'i[@name="ISPIN"]')
        if entry is None:
            return None

        ispin = self._convert_i(entry)

        return ispin

    def _fetch_ismearw(self, xml):
        """Fetch and set ismear using etree

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        ismear : int
            If ISMEAR is found it is returned.

        Notes
        -----
        Determines which smearing factor is used on the electrons.

        """

        entry = self._find(
            xml, './/parameters/separator[@name="electronic"]/'
            'separator[@name="electronic smearing"]/'
            'i[@name="ISMEAR"]')

        if entry is None:
            return None

        ismear = self._convert_i(entry)

        return ismear

    def _fetch_nbandsw(self, xml):
        """Fetch and set nbands using etree

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        nbands : int
            If NBANDS is found it is returned.

        Notes
        -----
        The number of bands used in the calculation.

        """

        entry = self._find(
            xml, './/parameters/separator[@name="electronic"]/'
            'i[@name="NBANDS"]')
        if entry is None:
            return None

        nbands = self._convert_i(entry)

        return nbands

    def _fetch_nelectw(self, xml):
        """Fetch and set nelect using etree.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        nelect : float
            If NELECT is found it is returned.

        Notes
        -----
        The number of electrons used in the calculation.

        """

        entry = self._find(
            xml, './/parameters/separator[@name="electronic"]/'
            'i[@name="NELECT"]')

        if entry is None:
            return None

        nelect = self._convert_f(entry)

        return nelect

    def _fetch_systemw(self, xml):
        """Fetch and set system using etree.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        system : string
            If SYSTEM is found it is returned.

        Notes
        -----
        A comment that can be specified in the INCAR file.

        """

        entry = self._find(
            xml, './/parameters/separator[@name="general"]/'
            'i[@name="SYSTEM"]')

        if entry is None:
            return None

        system = entry.text

        return system

    def _fetch_bornw(self, xml):
        """Fetch the Born effetive charges.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        born : ndarray
            A ndarray containing the born effective charge
            tensor for each atom.

        """

        entry = self._findall(
            xml, './/calculation/array[@name="born_charges"]/'
            'set/v')

        if entry is None:
            return None

        num_atoms = 0
        species = self._lattice['species']
        if species is None:
            # Try to fetch species again, if still none, we cannot parse it
            # and thus we cannot parse the entries in here either.
            species = self._fetch_speciesw(xml)
            if species is None:
                return None, None, None, None
        num_atoms = species.shape[0]

        born = self._convert_array2D_f(entry, 3)

        born = np.asarray(np.split(born, num_atoms))

        return born

    def _fetch_upfsw(self, xml, extract_all=False):
        """Fetch the unitcell, atomic positions, force and stress.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        extract_all : bool
            Determines which unitcell and positions to get.
            Defaults to the initial and last. If True, extract all.

        Returns
        -------
        cell, pos, force, stress : dict
            An dictionary containing ndarrays of the:
            | unitcells with the vectors as rows in AA.

            | positions with each position as a row in direct coordinates.

            | forces where each row is the force in eV/AA on each atom.

            | stress where each row is the stress matrix for the unitcell in kB.

        """

        cell = {}
        pos = {}
        force = {}
        stress = {}
        num_atoms = 0
        species = self._lattice['species']
        if species is None:
            # Try to fetch species again, if still none, we cannot parse it
            # and thus we cannot parse the entries in here either.
            species = self._fetch_speciesw(xml)
            if species is None:
                return None, None, None, None
        num_atoms = species.shape[0]

        if not extract_all:
            entry = self._findall(
                xml,
                './/structure[@name="finalpos"]/crystal/varray[@name="basis"]/v'
            )
            if entry is not None:
                cell[2] = self._convert_array2D_f(entry, 3)
            else:
                cell[2] = None
            entry = self._findall(
                xml,
                './/structure[@name="initialpos"]/crystal/varray[@name="basis"]/v'
            )
            if entry is not None:
                cell[1] = self._convert_array2D_f(entry, 3)
            else:
                cell[1] = None
            entry = self._findall(
                xml,
                './/structure[@name="finalpos"]/varray[@name="positions"]/v')
            if entry is not None:
                pos[2] = self._convert_array2D_f(entry, 3)
            else:
                pos[2] = None
            entry = self._findall(
                xml,
                './/structure[@name="initialpos"]/varray[@name="positions"]/v')
            if entry is not None:
                pos[1] = self._convert_array2D_f(entry, 3)
            else:
                pos[1] = None

            entry = self._findall(xml,
                                  './/calculation/varray[@name="stress"]/v')

            if entry is not None:
                stress[1] = self._convert_array2D_f(entry[0:3], 3)
                stress[2] = self._convert_array2D_f(entry[-3:], 3)
            else:
                stress[1] = None
                stress[2] = None

            entry = self._findall(xml,
                                  './/calculation/varray[@name="forces"]/v')
            if entry is not None:
                force[1] = self._convert_array2D_f(entry[0:num_atoms], 3)
                force[2] = self._convert_array2D_f(entry[-num_atoms:], 3)
            else:
                force[1] = None
                force[2] = None
        else:
            structures = self._findall(xml, './/calculation/structure')
            entrycell = self._findall(
                xml,
                './/calculation/structure/crystal/varray[@name="basis"]/v')
            entrypos = self._findall(
                xml, './/calculation/structure/varray[@name="positions"]/v')
            entryforce = self._findall(
                xml, './/calculation/varray[@name="forces"]/v')
            entrystress = self._findall(
                xml, './/calculation/varray[@name="stress"]/v')

            if structures is not None:
                num_calcs = len(structures)
            else:
                return None

            num_entrycell = 0
            num_entrypos = 0
            num_entryforce = 0
            num_entrystress = 0

            if entrycell is not None:
                num_entrycell = len(entrycell)
                cell[1] = self._convert_array2D_f(entrycell[0:3], 3)
                if num_entrycell > 3:
                    cell[2] = self._convert_array2D_f(entrycell[-3:], 3)
                else:
                    cell[2] = None
            else:
                cell[1] = None
                cell[2] = None

            if entrypos is not None:
                num_entrypos = len(entrypos)
                pos[1] = self._convert_array2D_f(entrypos[0:num_atoms], 3)
                if num_entrypos > 3:
                    pos[2] = self._convert_array2D_f(entrypos[-num_atoms:], 3)
                else:
                    pos[2] = None
            else:
                pos[1] = None
                pos[2] = None

            if entryforce is not None:
                num_entryforce = len(entryforce)
                force[1] = self._convert_array2D_f(entryforce[0:num_atoms], 3)
                if num_entryforce > 3:
                    force[2] = self._convert_array2D_f(entryforce[-num_atoms:],
                                                       3)
                else:
                    force[2] = None
            else:
                force[1] = None
                force[2] = None

            if entrystress is not None:
                num_entrystress = len(entrystress)
                stress[1] = self._convert_array2D_f(entrystress[0:3], 3)
                if num_entrystress > 3:
                    stress[2] = self._convert_array2D_f(entrystress[-3:], 3)
                else:
                    stress[2] = None
            else:
                stress[1] = None
                stress[2] = None

            max_entries = max(
                num_entrystress,
                max(num_entryforce, max(num_entrycell, num_entrypos)))
            if max_entries > 6:
                for calc in range(1, num_calcs):
                    basecell = calc * 3
                    basepos = calc * num_atoms
                    if entrycell is not None:
                        cell[calc + 1] = self._convert_array2D_f(
                            entrycell[basecell:basecell + 3], 3)
                    if entrypos is not None:
                        pos[calc + 1] = self._convert_array2D_f(
                            entrypos[basepos:basepos + num_atoms], 3)
                    if entryforce is not None:
                        force[calc + 1] = self._convert_array2D_f(
                            entryforce[basepos:basepos + num_atoms], 3)
                    if entrystress is not None:
                        stress[calc + 1] = self._convert_array2D_f(
                            entrystress[basecell:basecell + 3], 3)

        # If we still only have one entry, or number two is None, last and initial should
        # be the same, force them to be similar. We could do this earlier, but this is done
        # to keep it tidy, and if one in the future would like to utilize the returned None.
        if cell:
            if len(cell) == 1 or cell[2] is None:
                cell[2] = cell[1]
        if pos:
            if len(pos) == 1 or pos[2] is None:
                pos[2] = pos[1]
        if force:
            if len(force) == 1 or force[2] is None:
                force[2] = force[1]
        if stress:
            if len(stress) == 1 or stress[2] is None:
                stress[2] = stress[1]

        return cell, pos, force, stress

    def _fetch_speciesw(self, xml):
        """Fetch the atomic species

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        spec : ndarray
            An array containing the atomic species as a number.
            Organized in the same order as the atomic positions.

        """

        entry = self._findall(xml, './/atominfo/'
                              'array[@name="atoms"]/set/rc/c')

        if entry is None:
            return None

        spec = self._convert_species(entry[::2])

        return spec

    def _fetch_hessian(self, xml):
        """Fetch the hessian.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        hessian : ndarray
            An array containing the Hessian matrix.

        """

        entry = self._findall(
            xml, './/calculation/dynmat/'
            'varray[@name="hessian"]/v')

        if entry is None:
            return None

        species = self._lattice['species']
        if species is None:
            # Try to fetch species again, if still none, we cannot parse it
            # and thus we cannot parse the entries in here either.
            species = self._fetch_speciesw(xml)
            if species is None:
                return None
        num_atoms = species.shape[0]

        hessian = self._convert_array2D_f(entry, num_atoms * 3)

        return hessian

    def _fetch_dynmatw(self, xml):
        """Fetch the dynamical matrix data.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        dynmat : dict
            An dict containing the eigenvalues and eigenvectors.

        """

        entry = self._find(xml, './/calculation/dynmat/'
                           'v[@name="eigenvalues"]')

        if entry is None:
            return None

        species = self._lattice['species']
        if species is None:
            # Try to fetch species again, if still none, we cannot parse it
            # and thus we cannot parse the entries in here either.
            species = self._fetch_speciesw(xml)
            if species is None:
                return None
        num_atoms = species.shape[0]

        eigenvalues = self._convert_array_f(entry)

        entry = self._find(
            xml, './/calculation/dynmat/'
            'varray[@name="eigenvectors"]')

        if entry is None:
            return None

        eigenvectors = self._convert_array2D_f(entry, num_atoms * 3)

        dynmat = {'eigenvalues': eigenvalues, 'eigenvectors': eigenvectors}

        return dynmat

    def _fetch_kpointsw(self, xml,
                        path='kpoints/varray[@name="kpointlist"]/v'):
        """Fetch the kpoints.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        path : string
            The path in the XML file containing the k-point set to be extracted.

        Returns
        -------
        kpoints : ndarray
            An array containing the kpoints used in the calculation
            in direct coordinates.

        """

        entry = self._findall(xml, path)

        if entry is None:
            return None

        kpoints = self._convert_array2D_f(entry, 3)

        return kpoints

    def _fetch_kpointsww(self, xml, path='kpoints/varray[@name="weights"]/v'):
        """Fetch the kpoint weights.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        path : string
            The path in the XML file containing the k-point weights to be extracted.

        Returns
        -------
        kpointw : ndarray
            An array containing the kpoint weights used in the
            calculation.

        """

        entry = self._findall(xml, path)

        if entry is None:
            return None

        kpointsw = self._convert_array1D_f(entry)

        return kpointsw

    def _fetch_kpointdivw(self, xml):
        """Fetch the number of kpoints in each direction.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        kpointdiv : list
            An list containing the kpoint divisions used in the
            calculation for the full BZ.

        """

        entry = self._find(xml, 'kpoints/generation/v[@name="divisions"]')

        if entry is None:
            return None

        kpointdiv = self._convert_array_i(entry)

        return kpointdiv

    def _fetch_eigenvaluesw(self, xml):
        """Fetch the eigenvalues.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        eigenvalues, occupancies : tupple
            An tupple of dicts containing ndarrays containing the
            eigenvalues and occupancies for each spin, band and
            kpoint index.

        """

        # spin 1
        entry_ispin1 = self._findall(
            xml, './/calculation/eigenvalues/array/set/'
            'set[@comment="spin 1"]/set/r')

        # spin 2
        entry_ispin2 = self._findall(
            xml, './/calculation/eigenvalues/array/set/'
            'set[@comment="spin 2"]/set/r')

        # if we do not find spin 1 entries return right away
        if entry_ispin1 is None:
            return None, None

        eigenvalues, occupancies = self._extract_eigenvalues(
            entry_ispin1, entry_ispin2, len(self._lattice['kpoints']))

        return eigenvalues, occupancies

    def _fetch_eigenvalues_specificw(self, xml):
        """Fetch the eigenvalues at specific k-point grids.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        eigenvalues : ndarray
            An ndarray containing the
            eigenvalues for each spin, band and
            kpoint index.

        """

        # spin 1
        entry_ispin1 = self._findall(
            xml, './/calculation/eigenvalues/'
            'eigenvalues/array/set/'
            'set[@comment="spin 1"]/set/r')
        # spin 2
        entry_ispin2 = self._findall(
            xml, './/calculation/eigenvalues/'
            'eigenvalues/array/set/'
            'set[@comment="spin 2"]/set/r')
        if entry_ispin1 is not None:
            # Also extract the k-point grids
            self._data['kpoints'] = self._fetch_kpointsw(
                xml,
                path='.//calculation/eigenvalues/'
                'kpoints/varray[@name="kpointlist"]/v')

            self._data['kpointsw'] = self._fetch_kpointsww(
                xml,
                path='//calculation/eigenvalues/'
                'kpoints/varray[@name="weights"]/v')

        # if we do not find spin 1 entries return right away
        if entry_ispin1 is None:
            return None

        eigenvalues, _ = self._extract_eigenvalues(entry_ispin1, entry_ispin2,
                                                   len(self._data['kpoints']))

        return eigenvalues

    def _fetch_eigenvelocitiesw(self, xml):
        """Fetch the eigenvelocities and eigenvalues..

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        eigenvelocities : dict
            A dict with ndarrays containing the
            eigenvalues and eigenvelocities for each spin, band and
            kpoint index.
        kpoints : dict
            A dict with a ndarray containing the k-point in the full BZ
            on which the eigenvelocities were extracted.

        """

        # spin 1
        entry_ispin1 = self._findall(
            xml, './/calculation/eigenvelocities/'
            'eigenvalues/array/set/'
            'set[@comment="spin 1"]/set/r')

        # spin 2
        entry_ispin2 = self._findall(
            xml, './/calculation/eigenvelocities/'
            'eigenvalues/array/set/'
            'set[@comment="spin 2"]/set/r')

        # if we do not find spin 1 entries return right away
        if entry_ispin1 is None:
            return None

        self._data['kpoints'] = self._fetch_kpointsw(
            xml,
            path='.//calculation/eigenvelocities/'
            'kpoints/varray[@name="kpointlist"]/v')

        self._data['kpointsw'] = self._fetch_kpointsww(
            xml,
            path='//calculation/eigenvelocities/'
            'kpoints/varray[@name="weights"]/v')

        eigenvelocities = self._extract_eigenvelocities(
            entry_ispin1, entry_ispin2, len(self._data['kpoints']))

        return eigenvelocities

    def _fetch_projectorsw(self, xml):
        """Fetch the projectors.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        projectors : dict
            An dict containing ndarrays of the projectors
            for each atomic, spin, band and kpoint index.

        """

        # projectors spin 1
        entry_ispin1 = self._findall(
            xml, './/calculation/projected/array/set/'
            'set[@comment="spin1"]/set/set/r')

        # projectors spin 2
        entry_ispin2 = self._findall(
            xml, './/calculation/projected/array/set/'
            'set[@comment="spin2"]/set/set/r')

        # if we do not find spin 1 entries return right away
        if entry_ispin1 is None:
            return None

        projectors = self._extract_projectors(entry_ispin1, entry_ispin2)
        return projectors

    def _fetch_totensw(self, xml):
        """Fetch the total energies

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        eigenvalues, occupancies : tupple
            An tupple of two ndarrays containing the eigenvalues
            for each spin, band and kpoint index.

        """

        # fetch the energies for all electronic
        # steps, due to the fact that the number of steps is not the
        # same between each calculation we need to look at all the
        # children
        #
        # TODO: check in the future if it is faster to fetch all scstep
        # elements and then only how many scstep there is pr. calc
        # and sort from there
        #

        entries = self._findall(xml, './/calculation')

        if entries is None:
            return None

        # this most likely takes too long for very long fpmd calculations,
        # so consider putting in a flag that only extract the
        # energies from each step in the calculation and not the scsteps as
        # well
        energies = {}

        for index, calc in enumerate(entries):
            energies_pr_calc = {}
            for supported_energy, supported_key in _SUPPORTED_TOTAL_ENERGIES.items():
                data = self._findall(calc,
                                     './/scstep/energy/i[@name="' + supported_key + '"]')
                if data is None:
                    return None
                data = self._convert_array1D_f(data)
                energies_pr_calc[supported_energy] = data

                # now fetch the final entry outside the sc steps for each calculation
                # this term might have been corrected and scaled compared to the final sc energy
                # extrapolated energy
                data = self._findall(calc,
                                     './energy/i[@name="' + supported_key + '"]')
                if data is None:
                    return None
                data = self._convert_f(data[0])
                energies_pr_calc[supported_energy+'_final'] = data
            energies[index + 1] = energies_pr_calc

        return energies

    def _fetch_dosw(self, xml):
        """Fetch the density of states.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        dos, dos_specific : dicts
            Dictionaries with ndarrays containing the energies, total and
            integrated density of states for the regular and specific k-point grid, respectively.

        """

        # fetch the Fermi level
        entry = self._find(xml, './/calculation/dos/i[@name="efermi"]')

        if entry is not None:
            fermi_level = self._convert_f(entry)
        else:
            fermi_level = None

        # spin 1
        entry_total_ispin1 = self._findall(
            xml, './/calculation/dos/total/array/set/set[@comment="spin 1"]/r')

        # spin 2
        entry_total_ispin2 = self._findall(
            xml, './/calculation/dos/total/array/set/set[@comment="spin 2"]/r')

        # partial spin 1
        entry_partial_ispin1 = self._findall(
            xml,
            './/calculation/dos/partial/array/set/set/set[@comment="spin 1"]/r'
        )

        # partial spin 2
        entry_partial_ispin2 = self._findall(
            xml,
            './/calculation/dos/partial/array/set/set/set[@comment="spin 2"]/r'
        )

        # if no entries for spin 1, eject right away
        if entry_total_ispin1 is None:
            return None, None

        num_atoms = 0
        species = self._lattice['species']
        if species is None:
            # Try to fetch species again, if still none, we cannot parse it
            # and thus we cannot parse the entries in here either.
            species = self._fetch_speciesw(xml)
            if species is None:
                return None, None
        num_atoms = species.shape[0]

        dos = self._extract_dos(entry_total_ispin1, entry_total_ispin2,
                                entry_partial_ispin1, entry_partial_ispin2,
                                fermi_level, num_atoms)

        # Now extract the density of states for specific k-point grids
        # fetch the Fermi level again as this can be different
        entry = self._find(
            xml,
            './/calculation/dos[@comment="interpolated"]/i[@name="efermi"]')

        if entry is not None:
            fermi_level = self._convert_f(entry)
        else:
            fermi_level = None

        # spin 1
        entry_total_ispin1 = self._findall(
            xml,
            './/calculation/dos[@comment="interpolated"]/total/array/set/set[@comment="spin 1"]/r'
        )

        # spin 2
        entry_total_ispin2 = self._findall(
            xml,
            './/calculation/dos[@comment="interpolated"]/total/array/set/set[@comment="spin 2"]/r'
        )

        # partial spin 1
        entry_partial_ispin1 = self._findall(
            xml,
            './/calculation/dos[@comment="interpolated"]/partial/array/set/set/set[@comment="spin 1"]/r'
        )

        # partial spin 2
        entry_partial_ispin2 = self._findall(
            xml,
            './/calculation/dos[@comment="interpolated"]/partial/array/set/set/set[@comment="spin 2"]/r'
        )

        # if no entries for spin 1, eject right away
        if entry_total_ispin1 is None:
            return dos, None

        dos_specific = self._extract_dos(entry_total_ispin1,
                                         entry_total_ispin2,
                                         entry_partial_ispin1,
                                         entry_partial_ispin2, fermi_level,
                                         num_atoms)

        return dos, dos_specific

    def _extract_dos(self, entry_total_ispin1, entry_total_ispin2,
                     entry_partial_ispin1, entry_partial_ispin2, fermi_level,
                     num_atoms):
        """Extract the density of states.

        Parameters
        ----------
        entry_total_spin1 : list
            A list containing ElementTree objects of the total density of states entry for
            spin channel 1.
        entry_total_spin2 : list
            A list containing ElementTree objects of the total density of states entry for
            spin channel 2.
        entry_total_spin1 : list
            A list containing ElementTree objects of the partial density of states entry for
            spin channel 1.
        entry_total_spin1 : list
            A list containing ElementTree objects of the partial density of states entry for
            spin_channel 2.
        fermi_level : float
            The Fermi level in eV.
        num_atoms : int
            The number of atoms.

        Returns
        -------
        dos : dict
            A dict of ndarrays containing the energies, total and
            integrated density of states.

        """
        if entry_total_ispin2:
            dos = {}
            dos = {'up': None, 'down': None}
            dos_ispin = self._convert_array2D_f(entry_total_ispin1, 3)
            _dosup = {}
            _dosdown = {}
            enrgy = dos_ispin[:, 0]
            _dosup['total'] = dos_ispin[:, 1]
            _dosup['integrated'] = dos_ispin[:, 2]
            # check if partial exists
            if entry_partial_ispin1:
                dos_ispin = self._convert_array2D_f(entry_partial_ispin1, 10)
                # do not need the energy term (similar to total)
                _dosup['partial'] = np.asarray(
                    np.split(dos_ispin[:, 1:10], num_atoms))
            else:
                _dosup['partial'] = None
            dos['up'] = _dosup
            dos_ispin = self._convert_array2D_f(entry_total_ispin2, 3)
            _dosdown['total'] = dos_ispin[:, 1]
            _dosdown['integrated'] = dos_ispin[:, 2]
            if entry_partial_ispin2:
                dos_ispin = self._convert_array2D_f(entry_partial_ispin2, 10)
                # do not need the energy term (similar to total)
                _dosdown['partial'] = np.asarray(
                    np.split(dos_ispin[:, 1:10], num_atoms))
            else:
                _dosdown['partial'] = None
            dos['down'] = _dosdown
            dos['total'] = {'fermi_level': fermi_level, 'energy': enrgy}
        else:
            dos = {}
            dos = {'total': None}
            dos_ispin = self._convert_array2D_f(entry_total_ispin1, 3)
            _dos = {}
            _dos['energy'] = dos_ispin[:, 0]
            _dos['total'] = dos_ispin[:, 1]
            _dos['integrated'] = dos_ispin[:, 2]
            # check if partial exists
            if entry_partial_ispin1:
                dos_ispin = self._convert_array2D_f(entry_partial_ispin1, 10)
                # do not need the energy term (similar to total)
                _dos['partial'] = np.asarray(
                    np.split(dos_ispin[:, 1:10], num_atoms))
            else:
                _dos['partial'] = None
            _dos['fermi_level'] = fermi_level
            dos['total'] = _dos

        return dos

    def _fetch_dielectricsw(self, xml, method='dft', transfer=None):
        """ Fetch the dielectric function from the VASP XML file

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        method : {'dft', 'qp', 'bse'}, optional
            What method was used to obtain the dielectric function. VASP
            uses different output between DFT, QP and BSE calculations
            (defaults to 'dft').
        transfer : {'density', 'current'}, optional
            Which dielectric function do you want? Density-density or
            current-current? Defaults to the density-density.

        Returns
        -------
        diel_imag : (N,6) list of list of float
            If `method` is 'dft'.
            The imaginary dielectric function for N energies for the
            xx, yy, zz, xy, yz and zx component, respectively.
        diel_real : (N,6) list of float
            If `method` is 'dft'.
            The real dielectric function for N energies for the
            xx, yy, zz, xy, yz and zx component, respectively.
        diel_imag_mac : (N,6) list of list of float
            If `method` is 'qp'.
            The imaginary part of the macroscopic dielectric function.
            See `diel_imag` for layout.
        diel_real_mac : (N,6) list of list of float
            If `method` is 'qp'.
            The real part of the polarized dielectric function.
            See `diel_imag` for layout.
        diel_imag_pol : (N,6) list of list of float
            If `method` is 'qp'.
            The imaginary part of the polarized dielectric function.
            See `diel_imag` for layout.
        diel_real_pol : (N,6) list of list of floa
            If `method` is 'qp'.
            The real part of the polarized dielectric function.
            See `diel_imag` for layout.
        diel_imag_invlfrpa : (N,6) list of list of float
            If `method` is 'qp'.
            The imaginary part of the inverse dielectric function with
            local field effects on the RPA level.
            See `diel_imag` for layout.
        diel_real_invlfrpa : (N,6) list of list of float
            If `method` is 'qp'.
            The real part of the inverse dielectric function with
            local field effects on the RPA level.
            See `diel_imag` for layout.
        diel_imag : (N,6) list of list of float
            If `method` is 'bse'.
            The imaginary part of the BSE dielectric function.
            See `diel_imag` above for layout.
            diel_real : (N,6) list of list of float
            If `method` is 'bse'.
            The real part of the BSE dielectric function.
            See `diel_imag` at the top for layout.
        epsilon : (3,3) list of list of float

        """

        if method == 'dft':
            diel = {}
            if transfer == 'density':
                tag = 'dielectricfunction[@comment="density-density"]'
            elif transfer == 'current':
                tag = 'dielectricfunction[@comment="current-current"]'
            else:
                tag = 'dielectricfunction'

            # imaginary part
            entry = self._findall(
                xml, './/calculation/' + tag + '/imag/array/set/r')
            if entry is None:
                diel['imag'] = None
                diel['energy'] = None
            else:
                data = self._convert_array2D_f(entry, 7)
                diel['energy'] = data[:, 0]
                diel['imag'] = data[:, 1:7]

            # real part
            entry = self._findall(
                xml, './/calculation/' + tag + '/real/array/set/r')
            if entry is None:
                diel['real'] = None
            else:
                data = self._convert_array2D_f(entry, 7)
                diel['real'] = data[:, 1:7]

            # epsilon part
            entry = self._findall(xml,
                                  './/calculation/varray[@name="epsilon"]/v')
            if entry is not None:
                diel['epsilon'] = self._convert_array2D_f(entry, 3)
            else:
                diel['epsilon'] = None

            # ionic epsilon part
            entry = self._findall(
                xml, './/calculation/varray[@name="epsilon_ion"]/v')
            if entry is not None:
                diel['epsilon_ion'] = self._convert_array2D_f(entry, 3)
            else:
                diel['epsilon_ion'] = None

            return diel

        # if method == "qp":
        #     try:
        #         dielectric_xml = root.findall('dielectricfunction')
        #     except AttributeError:
        #         logger.error(
        #             "Did not find <dielectricfunction> tag in the current XML."
        #             "Exiting.")
        #         sys.exit(1)

        #     # first head of macroscopic
        #     diel_imag_xml = dielectric_xml[0].find(
        #         'imag').find('array').find('set')
        #     diel_imag_mac = []
        #     # first imag part
        #     for energy in diel_imag_xml.iter('r'):
        #         diel_imag_mac.append([float(x) for x in energy.text.split()])
        #     diel_real_xml = dielectric_xml[0].find(
        #         'real').find('array').find('set')
        #     diel_real_mac = []
        #     # then real part
        #     for energy in diel_real_xml.iter('r'):
        #         diel_real_mac.append([float(x) for x in energy.text.split()])

        #     # then polarized
        #     diel_imag_xml = dielectric_xml[1].find(
        #         'imag').find('array').find('set')
        #     diel_imag_pol = []
        #     # first imag part
        #     for energy in diel_imag_xml.iter('r'):
        #         diel_imag_pol.append([float(x) for x in energy.text.split()])
        #     diel_real_xml = dielectric_xml[1].find(
        #         'real').find('array').find('set')
        #     diel_real_pol = []
        #     # then real part
        #     for energy in diel_real_xml.iter('r'):
        #         diel_real_pol.append([float(x) for x in energy.text.split()])

        #     # then inverse macroscopic (including local field)
        #     diel_imag_xml = dielectric_xml[2].find(
        #         'imag').find('array').find('set')
        #     diel_imag_invlfrpa = []
        #     # first imag part
        #     for energy in diel_imag_xml.iter('r'):
        #         diel_imag_invlfrpa.append([float(x) for x in energy.text.split()])
        #     diel_real_xml = dielectric_xml[2].find(
        #             'real').find('array').find('set')
        #     diel_real_invlfrpa = []
        #     # then real part
        #     for energy in diel_real_xml.iter('r'):
        #         diel_real_invlfrpa.append([float(x) for x in energy.text.split()])
        #     return diel_imag_mac, diel_real_mac, diel_imag_pol, diel_real_pol, \
        #         diel_imag_invlfrpa, diel_real_invlfrpa

        # if method == "bse":
        #     try:
        #         dielectric_xml = root.find('dielectricfunction')
        #     except AttributeError:
        #         logger.error(
        #             "Did not find <dielectricfunction> tag in the current XML."
        #             "Exiting.")
        #         sys.exit(1)
        #     diel_imag_xml = dielectric_xml.find('imag').find('array').find('set')
        #     diel_imag = []
        #     # first imag part
        #     for energy in diel_imag_xml.iter('r'):
        #         diel_imag.append([float(x) for x in energy.text.split()])
        #     diel_real_xml = dielectric_xml.find('real').find('array').find('set')
        #     diel_real = []
        #     # then real part
        #     for energy in diel_real_xml.iter('r'):
        #         diel_real.append([float(x) for x in energy.text.split()])
        #     return diel_imag, diel_real

    def _extract_eigenvalues(self, spin1, spin2, num_kpoints):
        """Extract the eigenvalues.

        Parameters
        ----------
        spin1 : list
            A list of ElementTree object to be used for parsing of the
            ispin=1 entries.
        spin2 : list
            A list of ElementTree object to be used for parsing of the
            ispin=2 entries.
        num_kpoints : int
            The number of k-points to extract

        Returns
        -------
        eigenvalues, occupancies : tupple of dicts
            An tupple of two dicts containing ndarrays with the eigenvalues
            and occupancies for each band and kpoint index.

        """

        # then check if we have asigned ispin
        if self._parameters['ispin'] is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_ISPIN])
            sys.exit(self.ERROR_NO_ISPIN)

        # then check if we have asigned nbands
        if self._parameters['nbands'] is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_NBANDS])
            sys.exit(self.ERROR_NO_NBADS)

        # ispin
        ispin = self._parameters['ispin']

        # number of bands
        num_bands = self._parameters['nbands']

        # set dicts
        eigenvalues = {}
        occupancies = {}

        data = []

        if len(spin1) != num_bands * num_kpoints:
            self._logger.error(
                self.ERROR_MESSAGES[self.ERROR_MISMATCH_KPOINTS_NBANDS])
            sys.exit(self.ERROR_MISMATCH_KPOINTS_NBANDS)

        # check number of elements in first entry of spin1 (we assume all are equal)
        entries = len(spin1[0].text.split())
        if entries > 1:
            data.append(self._convert_array2D_f(spin1, entries))
        else:
            data.append(self._convert_array1D_f(spin1))
        data[0] = np.asarray(np.split(data[0], num_kpoints))
        if spin2 is not None:
            if len(spin2) != num_bands * num_kpoints:
                self._logger.error(
                    self.ERROR_MESSAGES[self.ERROR_MISMATCH_KPOINTS_NBANDS])
                sys.exit(self.ERROR_MISMATCH_KPOINTS_NBANDS)
            if entries > 1:
                data.append(self._convert_array2D_f(spin2, entries))
            else:
                data.append(self._convert_array1D_f(spin2))
            data[1] = np.asarray(np.split(data[1], num_kpoints))

        # convert to numpy arrays
        data = np.asarray(data)
        # swap axis if the band index should be before the kpoint index
        if not self._k_before_band:
            data = np.swapaxes(data, 1, 2)
        if spin2 is not None:
            if entries > 1:
                eigenvalues['up'] = np.ascontiguousarray(data[0, :, :, 0])
                eigenvalues['down'] = np.ascontiguousarray(data[1, :, :, 0])
                occupancies['up'] = np.ascontiguousarray(data[0, :, :, 1])
                occupancies['down'] = np.ascontiguousarray(data[1, :, :, 1])
            else:
                eigenvalues['up'] = np.ascontiguousarray(data[0, :, :])
                eigenvalues['down'] = np.ascontiguousarray(data[1, :, :])
        else:
            if entries > 1:
                eigenvalues['total'] = np.ascontiguousarray(data[0, :, :, 0])
                occupancies['total'] = np.ascontiguousarray(data[0, :, :, 1])
            else:
                eigenvalues['total'] = np.ascontiguousarray(data[0, :, :])
        if entries > 1:
            return eigenvalues, occupancies
        else:
            return eigenvalues, None

    def _extract_eigenvelocities(self, spin1, spin2, num_kpoints):
        """Extract the eigenvalues and eigenvelocities.

        Parameters
        ----------
        spin1 : list
            A list of ElementTree object to be used for parsing of the
            ispin=1 entries.
        spin2 : list
            A list of ElementTree object to be used for parsing of the
            ispin=2 entries.
        num_kpoints : dict
            The number of k-point in the full BZ
            on which the eigenvelocities were extracted.

        Returns
        -------
        eigenvelocities : dict
            A dict containing ndarrays with the eigenvalues
            and eigenvelocities for each band and kpoint index.

        """

        # check if we have asigned ispin
        if self._parameters['ispin'] is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_ISPIN])
            sys.exit(self.ERROR_NO_ISPIN)

        # then check if we have asigned nbands
        if self._parameters['nbands'] is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_NBANDS])
            sys.exit(self.ERROR_NO_NBANDS)

        # ispin
        ispin = self._parameters['ispin']

        # number of bands
        num_bands = self._parameters['nbands']

        # set dicts
        eigenvelocities = {}

        data = []
        if len(spin1) != num_bands * num_kpoints:
            self._logger.error(
                self.ERROR_MESSAGES[self.ERROR_MISMATCH_KPOINTS_NBANDS])
            sys.exit(self.ERROR_MISMATCH_KPOINTS_NBANDS)
        data.append(self._convert_array2D_f(spin1, 4))
        data[0] = np.asarray(np.split(data[0], num_kpoints))
        if spin2 is not None:
            if len(spin2) != num_bands * num_kpoints:
                self._logger.error(
                    self.ERROR_MESSAGES[self.ERROR_MISMATCH_KPOINTS_NBANDS])
                sys.exit(self.ERROR_MISMATCH_KPOINTS_NBANDS)
            data.append(self._convert_array2D_f(spin2, 4))
            data[1] = np.asarray(np.split(data[1], num_kpoints))

        # convert to numpy arrays
        data = np.asarray(data)
        # swap axis if the band index should be before the kpoint index
        if not self._k_before_band:
            data = np.swapaxes(data, 1, 2)
        if spin2 is not None:
            eigenvelocities['up'] = np.ascontiguousarray(data[0])
            eigenvelocities['down'] = np.ascontiguousarray(data[1])
        else:
            eigenvelocities['total'] = np.ascontiguousarray(data[0])

        return eigenvelocities

    def _extract_projectors(self, spin1, spin2):
        """Extract the projectors.

        Parameters
        ----------
        spin1 : list
            A list of ElementTree object to be used for parsing of the
            ispin=1 entries. Contains the projectors.
        spin2 : list
            A list of ElementTree object to be used for parsing of the
            ispin=2 entries. Contains the projectors.


        Returns
        -------
        projectors : dict
            A dict containing ndarrays with the projectors for each atom,
            band and kpoint index.

        """

        # first check if we have extracted the kpoints
        if self._lattice['kpoints'] is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_KPOINTS])
            sys.exit(self.ERROR_NO_KPOINTS)

        # then check if we have asigned ispin
        if self._parameters['ispin'] is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_ISPIN])
            sys.exit(self.ERROR_NO_ISPIN)

        # then check if we have asigned nbands
        if self._parameters['nbands'] is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_NBANDS])
            sys.exit(self.ERROR_NO_NBANDS)

        num_atoms = 0
        # also need the number of atoms if the projected values are supplied
        if self._lattice['species'] is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_NBANDS])
            sys.exit(self.ERROR_NO_NBANDS)
        else:
            num_atoms = self._lattice['species'].shape[0]

        # number of kpoints to disect the eigenvalue sets later
        num_kpoints = self._lattice['kpoints'].shape[0]

        # ispin
        ispin = self._parameters['ispin']

        # number of bands
        num_bands = self._parameters['nbands']

        # set dicts
        projectors = {}

        pdata = []
        if len(spin1) != num_bands * num_kpoints * num_atoms:
            self._logger.error(
                self.ERROR_MESSAGES[self.ERROR_MISMATCH_KPOINTS_NBANDS])
            sys.exit(self.ERROR_MISMATCH_KPOINTS_NBANDS)
        pdata.append(self._convert_array2D_f(spin1, 9))
        pdata[0] = np.asarray(np.split(pdata[0], num_kpoints))
        pdata[0] = np.asarray(np.split(pdata[0], num_bands, axis=1))
        if spin2 is not None:
            if len(spin2) != num_bands * num_kpoints * num_atoms:
                self._logger.error(
                    self.ERROR_MESSAGES[self.ERROR_MISMATCH_KPOINTS_NBANDS])
                sys.exit(self.ERROR_MISMATCH_KPOINTS_NBANDS)
            pdata.append(self._convert_array2D_f(spin2, 9))
            pdata[1] = np.asarray(np.split(pdata[1], num_kpoints))
            pdata[1] = np.asarray(np.split(pdata[1], num_bands, axis=1))

        # convert to numpy arrays
        pdata = np.asarray(pdata)
        # swap axis if the band index should be before the kpoint index
        # make sure atomic index is first
        pdata = np.swapaxes(pdata, 0, 3)
        pdata = np.swapaxes(pdata, 1, 3)
        if not self._k_before_band:
            pdata = np.swapaxes(pdata, 2, 3)

        if spin2 is not None:
            projectors['up'] = np.ascontiguousarray(pdata[:, 0, :, :])
            projectors['down'] = np.ascontiguousarray(pdata[:, 1, :, :])
        else:
            projectors['total'] = np.ascontiguousarray(pdata[:, 0, :, :])

        return projectors

    def _convert_array_i(self, entry):
        """Convert the input entry to numpy array

        Parameters
        ----------
        entry : string
            A string containing N integer elements separated by
            blank spaces.

        Returns
        -------
        data : ndarray
            | Dimension: (N)
            An array containing N integers.

        """

        data = None
        if entry is not None:
            try:
                data = np.fromstring(entry.text, sep=' ', dtype='intc')
            except ValueError as e:
                if str(e) == 'setting an array element with a sequence.':
                    self._logger.error(
                        self.ERROR_MESSAGES[self.ERROR_OVERFLOW])
                    sys.exit(self.ERROR_OVERFLOW)

        return data

    def _convert_array_f(self, entry):
        """Convert the input entry to numpy array

        Parameters
        ----------
        entry : string
            A string containing N float elements separated by
            blank spaces.

        Returns
        -------
        data : ndarray
            | Dimension: (N)
            An array containing N floats.

        """

        data = None
        if entry is not None:
            try:
                data = np.fromstring(entry.text, sep=' ', dtype='double')
            except ValueError as e:
                if str(e) == 'setting an array element with a sequence.':
                    self._logger.error(
                        self.ERROR_MESSAGES[self.ERROR_OVERFLOW])
                    sys.exit(self.ERROR_OVERFLOW)

        return data

    def _convert_array1D_i(self, entry):
        """Convert the input entry to numpy array

        Parameters
        ----------
        entry : list
            A list containing Element objects where each
            element is an integer

        Returns
        -------
        data : ndarray
            | Dimension: (N)
            An array containing N integers.

        """

        data = None
        if entry is not None:
            data = np.zeros(len(entry), dtype='intc')
        for index, element in enumerate(entry):
            try:
                data[index] = np.fromstring(element.text, sep=' ')
            except ValueError as e:
                if str(e) == 'setting an array element with a sequence.':
                    self._logger.error(
                        self.ERROR_MESSAGES[self.ERROR_OVERFLOW])
                    sys.exit(self.ERROR_OVERFLOW)

        return data

    def _convert_array1D_f(self, entry):
        """Convert the input entry to numpy array

        Parameters
        ----------
        entry : list
            A list containing Element objects where each
            element is a float

        Returns
        -------
        data : ndarray
            | Dimension: (N)
            An array containing N double elements.

        """

        data = None

        if entry is not None:
            data = np.zeros(len(entry), dtype='double')
        for index, element in enumerate(entry):
            try:
                data[index] = np.fromstring(element.text, sep=' ')
            except ValueError as e:
                if str(e) == 'setting an array element with a sequence.':
                    self._logger.error(
                        self.ERROR_MESSAGES[self.ERROR_OVERFLOW])
                    sys.exit(self.ERROR_OVERFLOW)

        return data

    def _convert_array2D_f(self, entry, dim):
        """Convert the input entry to numpy array

        Parameters
        ----------
        entry : list
            A list containing Element objects where each
            element is a float
        dim : int
            The dimension of the second index.

        Returns
        -------
        data : ndarray
            | Dimension: (N,M)
            An array containing N elements with M float
            elements.

        """

        data = None
        if entry is not None:
            data = np.zeros((len(entry), dim), dtype='double')

        for index, element in enumerate(entry):
            try:
                data[index] = np.fromstring(element.text, sep=' ')
            except ValueError as e:
                if str(e) == 'setting an array element with a sequence.':
                    self._logger.error(
                        self.ERROR_MESSAGES[self.ERROR_OVERFLOW])
                    sys.exit(self.ERROR_OVERFLOW)

        return data

    def _convert_f(self, entry):
        """Convert the input entry to a float.

        Parameters
        ----------
        entry : object
            An Element object containing an integer value.

        Returns
        -------
        data : float
            The float value.

        """

        data = None
        if entry.text is not None:
            if '****' in entry.text:
                self._logger.error(self.ERROR_MESSAGES[self.ERROR_OVERFLOW])
                sys.exit(self.ERROR_OVERFLOW)
            data = float(entry.text)

        return data

    def _convert_i(self, entry):
        """Convert the input entry to an integer.

        Parameters
        ----------
        entry : object
            An Element object containing an integer value.

        Returns
        -------
        data : int
            The integer value.

        """

        data = None
        if entry.text is not None:
            if '****' in entry.text:
                self._logger.error(self.ERROR_MESSAGES[self.ERROR_OVERFLOW])
                sys.exit(self.ERROR_OVERFLOW)
            data = int(entry.text)

        return data

    def _convert_species(self, entry):
        """Set the atomic species to the correct value

        Parameters
        ----------
        entry : list
            A list containing Element objects, where each
            element is one atomic species.

        Returns
        -------
        unitcell : ndarray
            | Dimension: (N,3)
            An array containing the positions of N atoms in
            direct units.

        """

        species = None
        if entry is not None:
            species = np.zeros(len(entry), dtype='intc')
        for index, spec in enumerate(entry):
            try:
                species[index] = constants.elements[entry[index].text.split()
                                                    [0].lower()]
            except KeyError:
                self._logger.warning(
                    self.ERROR_MESSAGES[self.ERROR_UNKNOWN_ELEMENT])
                sys.exit(self.ERROR_UNKNOWN_ELEMENT)

        return species

    def get_forces(self, status):

        forces = self._data['forces']
        if forces is None:
            return None
        self._check_calc_status(status)
        if status == 'initial':
            return forces[1]
        elif status == 'last':
            largest_key = max(self._data['forces'].keys())
            return forces[largest_key]
        elif status == 'all':
            return forces

    def get_stress(self, status):

        stress = self._data['stress']
        if stress is None:
            return None
        self._check_calc_status(status)
        if status == 'initial':
            return stress[1]
        elif status == 'last':
            largest_key = max(self._data['stress'].keys())
            return stress[largest_key]
        elif status == 'all':
            return stress

    def get_hessian(self):

        hessian = self._data['hessian']
        return hessian

    def get_dynmat(self):

        dynmat = self._data['dynmat']
        return dynmat

    def get_dielectrics(self):

        dielectrics = self._data['dielectrics']
        return dielectrics

    def get_epsilon(self):

        epsilon = self._data['dielectrics']['epsilon']
        return epsilon

    def get_epsilon_ion(self):

        epsilon_ion = self._data['dielectrics']['epsilon_ion']
        return epsilon_ion

    def get_fermi_level(self):

        fermi_level = self._data['dos']['total']['fermi_level']
        return fermi_level

    def get_fermi_level_specific(self):

        fermi_level_specific = self._data['dos_specific']['total'][
            'fermi_level']
        return fermi_level_specific

    def get_born(self):

        born = self._data['born']
        return born

    def get_unitcell(self, status):

        unitcell = self._lattice['unitcell']
        if unitcell is None:
            return None
        self._check_calc_status(status)
        if status == 'initial':
            return unitcell[1]
        elif status == 'last':
            largest_key = max(self._lattice['unitcell'].keys())
            return unitcell[largest_key]
        elif status == 'all':
            return unitcell

    def get_positions(self, status):

        positions = self._lattice['positions']
        if positions is None:
            return None
        self._check_calc_status(status)
        if status == 'initial':
            return positions[1]
        elif status == 'last':
            largest_key = max(self._lattice['positions'].keys())
            return positions[largest_key]
        elif status == 'all':
            return positions

    def get_species(self):

        species = self._lattice['species']
        return species

    def get_lattice(self, status):

        species = self.get_species()
        unitcell = self.get_unitcell(status)
        positions = self.get_positions(status)
        return {
            'unitcell': unitcell,
            'positions': positions,
            'species': species
        }

    def get_kpoints(self):

        kpoints = self._lattice['kpoints']
        return kpoints

    def get_kpointsw(self):

        kpointsw = self._lattice['kpointsw']
        return kpointsw

    def get_energies(self, status, etype=None, nosc=True):

        if etype is None:
            etype = ['energy_extrapolated']
        # Check if the supplied etype is in the support list
        for item in etype:
            if item not in _SUPPORTED_TOTAL_ENERGIES.keys():
                raise ValueError(f'The supplied total energy type: {item} is not supported.')

        return self._get_energies(status, etype, nosc)

    def _get_energies(self, status, etype, nosc):

        enrgies = self._data['totens']
        if enrgies is None:
            return None
        self._check_calc_status(status)
        energies = {}
        # We can have different number of electronic steps for each ionic step, so
        # the array would be staggered. In order to save on storage and still utilize
        # regular NumPy array functions, we flatten the array and store a separate array
        # that keeps track of the number of electronic steps per ionic step.
        for item in etype:
            # For the energies inside the electronic step sections.
            energies_per_etype = np.array([])
            # For the final energy available inside the calculations (ionic steps) sections after closure
            # of the electronic steps and applying corrections.
            energy_per_etype = np.array([])
            electronic_steps = np.array([], dtype=int)
            steps = 1
            if status == 'initial':
                # Initial ionic step
                energy_per_etype = np.append(energy_per_etype, enrgies[1][item+'_final'])
                if nosc:
                    e = enrgies[1][item][-1]
                else:
                    e = enrgies[1][item]
                    steps = len(e)
                energies_per_etype = np.append(energies_per_etype, e)
                electronic_steps = np.append(electronic_steps, steps)
            elif status == 'last':
                # Last ionic step
                largest_key = max(enrgies.keys())
                energy_per_etype = np.append(energy_per_etype, enrgies[largest_key][item+'_final'])
                if nosc:
                    e = enrgies[largest_key][item][-1]
                else:
                    e = enrgies[largest_key][item]
                    steps = len(e)
                energies_per_etype = np.append(energies_per_etype, e)
                electronic_steps = np.append(electronic_steps, steps)
            elif status == 'all':
                # For all the ionic steps
                _energies_per_etype = sorted(enrgies.items())
                for index, element in _energies_per_etype:
                    energy_per_etype = np.append(energy_per_etype, element[item+'_final'])
                    if nosc:
                        e = element[item][-1]
                    else:
                        e = element[item]
                        steps = len(e)
                    energies_per_etype = np.append(energies_per_etype, e)
                    electronic_steps = np.append(electronic_steps, steps)
            energies[item+'_final'] = energy_per_etype
            energies[item] = energies_per_etype

        energies['electronic_steps'] = electronic_steps

        return energies

    def get_dos(self):

        dos = self._data['dos']
        return dos

    def get_dos_specific(self):

        dos_specific = self._data['dos_specific']
        return dos_specific

    def get_eigenvalues(self):

        eigenvalues = self._data['eigenvalues']
        return eigenvalues

    def get_eigenvalues_specific(self):

        eigenvalues_specific = self._data['eigenvalues_specific']
        return eigenvalues_specific

    def get_eigenvelocities(self):

        eigenvelocities = self._data['eigenvelocities']
        return eigenvelocities

    def get_kpoints_specific(self):

        kpoints_specific = self._data['kpoints']
        return kpoints_specific

    def get_kpointsw_specific(self):

        kpointsw_specific = self._data['kpointsw']
        return kpointsw_specific

    def get_occupancies(self):

        occupancies = self._data['occupancies']
        return occupancies

    def get_projectors(self):

        projectors = self._data['projectors']
        return projectors

    def get_dict(self):

        dictionary = {
            'parameters': self._parameters,
            'lattice': self._lattice,
            'data': self._data
        }

        return dictionary

    def get_parameters(self):

        parameters = self._parameters
        return parameters

    def get_version(self):
        from re import compile

        version = self._version.strip()
        # The version entry in the xml file is typically of the form
        # X.Y.Z or X.Y.Z.something so only keep three integers and two dots
        pattern = r'(\d+)\.(\d+)\.(\d+)'
        match = compile(pattern)

        return match.search(version).group(0)
    
    def _check_calc_status(self, status):
        allowed_entries = ['initial', 'last', 'all']
        if status not in allowed_entries:
            self._logger.error(
                self.ERROR_MESSAGES[self.ERROR_UNSUPPORTED_STATUS] +
                ' Please use any of the following values ' +
                str(allowed_entries))
            sys.exit(self.ERROR_UNSUPPORTED_STATUS)

    def _find(self, xml, locator):
        """Wrapper to check if the request returns something.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        locator : string
            The locator string to try.

        Returns
        -------
        entry : object
            An Element object if something is found, otherwise it
            returns None.

        """

        entry = xml.find(locator)

        if entry is None:
            return None
        else:
            return entry

    def _findall(self, xml, locator):
        """Wrapper to check if the request returns something.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        locator : string
            The locator string to try.

        Returns
        -------
        entry : object
            An Element object if something is found, otherwise it
            returns None.

        """

        entry = xml.findall(locator)

        if not entry:
            return None
        else:
            return entry

    def _file_size(self):
        """Returns the file size of a file.

        Returns
        -------
        The file size in megabytes.

        """

        if self._file_path is None and self._file_handler is None:
            self._logger.error(
                self.ERROR_MESSAGES[self.ERROR_ONLY_ONE_ARGUMENT])
            return None
        if self._file_handler is None:
            # check if file exists
            if not utils.file_exists(self._file_path, logger=self._logger):
                self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_SIZE])
                return None

            file_size = os.stat(self._file_path).st_size
        else:
            file_size = os.fstat(self._file_handler.fileno()).st_size

        return file_size / 1048576.0

    def _check_xml(self):
        """Do a primitive check of XML file to see if it is
        truncated.

        Returns
        -------
        bool
            True if xml is truncated, False otherwise.

        """

        if self._file_handler is not None:
            handler = self._file_handler
            mapping = mmap.mmap(handler.fileno(), 0, prot=mmap.PROT_READ)
        else:
            handler = open(self._file_path)
            with handler as source:
                mapping = mmap.mmap(source.fileno(), 0, prot=mmap.PROT_READ)

        last_line = mapping[mapping.rfind(b'\n', 0, -1) + 1:]
        if last_line == '</modeling>\n':
            return False
        else:
            return True
