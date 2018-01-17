#!/usr/bin/python
import sys
import os
import numpy as np
import logging

import constants

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
                    logging.error(
                        "Failed to import ElementTree. Exiting.")
                    sys.exit(1)


class XmlParser(object):

    def __init__(self, logger, file_path):
        """Initialize the XmlParser by first trying the lxml and
        fall back to the standard ElementTree if that is not present.

        Parameters
        ----------
        logger : object
            The logger to be used for outputting messages
        file_path : string
            The path of the XML file that is to be opened

        Notes
        -----
        lxml should be used and is required for large files
        """

        self._file_path = file_path
        self._sizecutoff = 500

        self._is_forces = False
        self._is_stress = False
        self._is_positions = False
        self._is_symbols = False
        self._is_basis = False
        self._is_energy = False
        self._is_k_weights = False
        self._is_eigenvalues = False
        self._is_epsilon = False
        self._is_born = False
        self._is_efermi = False

        self._is_v = False
        self._is_i = False
        self._is_rc = False
        self._is_c = False
        self._is_set = False
        self._is_r = False

        self._is_scstep = False
        self._is_structure = False
        self._is_projected = False
        self._is_proj_eig = False

        self._all_forces = []
        self._all_stress = []
        self._all_points = []
        self._all_lattice = []
        self._symbols = []
        self._all_energies = []
        self._born = []
        self._forces = None
        self._stress = None
        self._points = None
        self._lattice = None
        self._energies = None
        self._epsilon = None
        self._born_atom = None
        self._efermi = None
        self._k_weights = None
        self._eigenvalues = None
        self._eig_state = [0, 0]
        self._projectors = None
        self._proj_state = [0, 0, 0]

        # dictionaries
        self._parameters = {"symprec": None}
        self._lattice = {"unitcell": None,
                         "species": None,
                         "positions:": None}

        # parse
        self._parse()

    def _parse(self):
        """Perform the actual parsing

        """

        # Check size of the XML file. For large files we need to
        # perform event driven parsing. For smaller files this is
        # not necessary and is too slow.
        if self._file_size(self._file_path) < self._sizecutoff:
            self._parsew()
            # self._parsee()
        else:
            self._parsee()

    def _parsew(self):
        """Performs parsing on the whole XML files. For smaller files

        """

        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)
        logger.debug("Running parsew.")

        # now open the complete file
        self._check_file(self._file_path)
        vaspxml = etree.parse(self._file_path)

        # let us start to parse the content
        self._fetch_symprecw(vaspxml)
        self._fetch_latticew(vaspxml)

    def _parsee(self):
        """Performs parsing in an event driven fashion on the XML file.
        For bigger files.

        """

        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)
        logger.debug("Running parsee.")

        # bool to control extraction of content
        extract_parameters = False

        for event, element in etree.iterparse(self._file_path, events=("start", "end")):
            # first we detect sections that we want to read

            # parameter section
            if event == "start" and element.tag == "parameters":
                extract_parameters = True
            if event == "end" and element.tag == "parameters":
                extract_parameters = False

            if extract_parameters:
                try:
                    if element.attrib["name"] == "SYMPREC":
                        self._parameters[
                            "symprec"] = self._convert_symprec(element.text)
                except KeyError:
                    pass

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

        entry = xml.find('.//parameters/separator[@name="symmetry"]/'
                         'i[@name="SYMPREC"]')
        symprec = self._convert_symprec(entry)

        return symprec

    def _fetch_latticew(self, xml, status="final"):
        """Fetch and set the lattice

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        status : {"initial", "final"}
            Determines which lattice to get. Defaults to the
            final positions.

        Returns
        -------
        lattice : dict
            A dictionary containing the lattice requested

        """

        # first fetch and set the unitcell
        status = status + "pos"
        print './/structure[@name="' + status + '"] / crystal/varray[@name="basis"]/v'
        entry = xml.findall(
            './/structure[@name="' + status + '"] / crystal/'
            'varray[@name="basis"]/v')

        self._lattice["unitcell"] = self._convert_unitcell(entry)

        # then the atomic positions and species
        entry = xml.findall(
            './/structure[@name="finalpos"]/'
            'varray[@name="positions"]/v')
        entryspec = xml.findall('.//atominfo/'
                                'array[@name="atoms"]/set/rc')

        pos, spec = self._convert_positions_species(entry, entryspec)
        self._lattice["positions"] = pos
        self._lattice["species"] = spec

    def _convert_symprec(self, entry):
        """Set the symprec to correct value

        """

        symprec = None
        if entry.text is not None:
            symprec = float(entry.text)

        return symprec

    def _convert_unitcell(self, entry):
        """Set the unitcell to the correct value

        """

        unitcell = None
        if entry is not None:
            unitcell = np.zeros((3, 3))
            for index, unitcell_vector in enumerate(entry):
                unitcell[index] = np.fromstring(entry.text, sep=' ')

        return unitcell

    def _convert_positions_species(self, entrypos, entryspec):
        """Set the atomic positions to the correct value

        """

        positions = None
        species = None
        if entrypos is not None:
            positions = np.zeros((len(entrypos), 3),
                                 dtype='double')
        if entryspec is not None:
            species = np.zeros(len(entrypos), dtype='intc', order='C')
        for index, position in enumerate(entrypos):
            positions[index] = np.fromstring(position.text, sep=' ')
            species[index] = constants.elements[
                entryspec[index][0].text.split()[0].lower()]

        return positions

    def get_forces(self):
        return np.array(self._all_forces)

    def get_stress(self):
        return np.array(self._all_stress)

    def get_epsilon(self):
        return np.array(self._epsilon)

    def get_efermi(self):
        return self._efermi

    def get_born(self):
        return np.array(self._born)

    def get_points(self):
        return np.array(self._all_points)

    def get_lattice(self):
        return np.array(self._all_lattice)

    def get_symbols(self):
        return self._symbols

    def get_energies(self):
        return np.array(self._all_energies)

    def get_k_weights(self):
        return self._k_weights

    def get_eigenvalues(self):
        return self._eigenvalues

    def get_projectors(self):
        return self._projectors

    def _check_file(self, file_path):
        """
        Check if a file exists

        Parameters
        ----------
        file_path : string
            The path of the file to be checked.

        Returns
        -------
        None

        """

        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)
        logger.debug("Running check_file.")

        if not os.path.isfile(file_path):
            logger.error(file_path + "was not found. Exiting.")
            sys.exit(1)

    def _file_size(self, file_path):
        """Returns the file size of a file.

        Parameters
        ----------
        filepath : string
            The file path to the file to be checked.

        Returns
        -------
        The file size in megabytes.

        """
        file_info = os.stat(file_path)
        file_size = file_info.st_size
        return file_size / 1048576.0
