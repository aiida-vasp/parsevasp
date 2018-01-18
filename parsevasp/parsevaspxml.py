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

        # set logger
        self._logger = logger

        # dictionaries that contain the output of the parsing
        self._parameters = {"symprec": None}
        self._lattice = {"unitcell": None,
                         "species": None,
                         "positions": None,
                         "kpoints": None,
                         "kpointsw": None,
                         "kpointdiv": None}

        # parse parse parse
        self._parse()

    def _parse(self):
        """Perform the actual parsing

        """

        # Check size of the XML file. For large files we need to
        # perform event driven parsing. For smaller files this is
        # not necessary and is too slow.
        if self._file_size(self._file_path) < self._sizecutoff:
            #self._parsew()
            self._parsee()
        else:
            self._parsee()

        print self._lattice

    def _parsew(self):
        """Performs parsing on the whole XML files. For smaller files

        """

        self._logger.debug("Running parsew.")

        # now open the complete file
        self._check_file(self._file_path)
        vaspxml = etree.parse(self._file_path)

        # let us start to parse the content
        self._parameters["symprec"] = self._fetch_symprecw(vaspxml)
        self._lattice["unitcell"] = self._fetch_unitcellw(vaspxml)
        self._lattice["species"] = self._fetch_speciesw(vaspxml)
        self._lattice["positions"] = self._fetch_positionsw(vaspxml)
        self._lattice["kpoints"] = self._fetch_kpointsw(vaspxml)
        self._lattice["kpointsw"] = self._fetch_kpointsww(vaspxml)
        self._lattice["kpointdiv"] = self._fetch_kpointdiv(vaspxml)
        
    def _parsee(self):
        """Performs parsing in an event driven fashion on the XML file.
        Slower, but suitable for bigger files.

        """

        # set logger
        self._logger.debug("Running parsee.")

        # helper list
        data = []

        # bool to control extraction of content
        extract_parameters = False
        extract_latticedata = False
        extract_unitcell = False
        extract_positions = False
        extract_species = False
        extract_kpointdata = False
        extract_kpoints = False
        extract_kpointsw = False
        extract_kpointdiv = False

        # need to feed this in later
        status = "finalpos"

        for event, element in etree.iterparse(self._file_path, events=("start", "end")):
            # set extraction points (what to read and when to read it)
            # here we also set the relevant data elements when the tags
            # close when they contain more than one element
            if event == "start" and element.tag == "parameters":
                extract_parameters = True
            if event == "end" and element.tag == "parameters":
                extract_parameters = False
            if event == "start" and element.tag == "structure" \
               and element.attrib.get("name") == status:
                extract_latticedata = True
            if event == "end" and element.tag == "structure" \
               and element.attrib.get("name") == status:
                extract_latticedata = False
            if event == "start" and element.tag == "array" \
               and element.attrib.get("name") == "atoms":
                extract_species = True
            if event == "end" and element.tag == "array" \
               and element.attrib.get("name") == "atoms":
                # only need every other element (element, not atomtype)
                self._lattice["species"] = self._convert_species(data[::2])
                data = []
                extract_species = False
            if event == "start" and element.tag == "kpoints":
                extract_kpointdata = True
            if event == "end" and element.tag == "kpoints":
                extract_kpointdata = False

            # now fetch the data
            if extract_parameters:
                try:
                    if event == "start" and element.attrib["name"] == "SYMPREC":
                        self._parameters["symprec"] = self._convert_f(element)
                except KeyError:
                    pass

            if extract_latticedata:
                # print event, element.tag, element.text, element.attrib
                if event == "start" and element.tag == "varray" \
                   and element.attrib.get("name") == "basis":
                    extract_unitcell = True
                if event == "end" and element.tag == "varray" \
                   and element.attrib.get("name") == "basis":
                    self._lattice["unitcell"] = self._convert_array2D_f(data)
                    data = []
                    extract_unitcell = False
                if event == "start" and element.tag == "varray" \
                   and element.attrib.get("name") == "positions":
                    extract_positions = True
                if event == "end" and element.tag == "varray" \
                   and element.attrib.get("name") == "positions":
                    self._lattice["positions"] = self._convert_array2D_f(data)
                    data = []
                    extract_positions = False

                if extract_unitcell:
                    if event == "start" and element.tag == "v":
                        data.append(element)
                if extract_positions:
                    if event == "start" and element.tag == "v":
                        data.append(element)

            if extract_species:
                if event == "start" and element.tag == "c":
                    data.append(element)

            if extract_kpointdata:
                try:
                    if event == "start" and element.tag == "v" and element.attrib["name"] == "divisions":
                        self._lattice["kpointdiv"] = self._convert_array_i(element)
                except KeyError:
                    pass

                if event == "start" and element.tag == "varray" \
                   and element.attrib.get("name") == "kpointlist":
                    extract_kpoints = True
                if event == "end" and element.tag == "varray" \
                   and element.attrib.get("name") == "kpointlist":
                    self._lattice["kpoints"] = self._convert_array2D_f(data)
                    data = []
                    extract_kpoints = False
                if event == "start" and element.tag == "varray" \
                   and element.attrib.get("name") == "weights":
                    extract_kpointsw = True
                if event == "end" and element.tag == "varray" \
                   and element.attrib.get("name") == "weights":
                    self._lattice["kpointsw"] = self._convert_array1D_f(data)
                    data = []
                    extract_kpointsw = False
                if extract_kpoints:
                    if event == "start" and element.tag == "v":
                        data.append(element)
                if extract_kpointsw:
                    if event == "start" and element.tag == "v":
                        data.append(element)
                    

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
        symprec = self._convert_f(entry)

        return symprec

    def _fetch_unitcellw(self, xml, status="final"):
        """Fetch the unitcell

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        status : {"initial", "final"}
            Determines which unitcell to get. Defaults to the
            final unitcell.

        Returns
        -------
        cell : ndarray
            An array containing the unitcell vectors as rows in
            units of AA.

        """

        entry = xml.findall(
            './/structure[@name="' + status + 'pos"]/crystal/'
            'varray[@name="basis"]/v')

        cell = self._convert_array2D_f(entry)

        return cell

    def _fetch_positionsw(self, xml, status="final"):
        """Fetch the atomic positions.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        status : {"initial", "final"}
            Determines which atomic positions to get. Defaults to the
            final positions.

        Returns
        -------
        pos : ndarray
            An array containing the atomic positions in direct coordinates.

        """
        
        entry = xml.findall(
            './/structure[@name="' + status + 'pos"]/'
            'varray[@name="positions"]/v')

        pos = self._convert_array2D_f(entry)

        return pos
    
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

        entry = xml.findall('.//atominfo/'
                            'array[@name="atoms"]/set/rc/c')[::2]

        spec = self._convert_species(entry)

        return spec

    def _fetch_kpointsw(self, xml):
        """Fetch the kpoints.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        kpoints : ndarray
            An array containing the kpoints used in the calculation 
            in direct coordinates.

        """

        entry = xml.findall(
            'kpoints/varray[@name="kpointlist"]/v')
                
        kpoints = self._convert_array2D_f(entry)

        return kpoints

    def _fetch_kpointsww(self, xml):
        """Fetch the kpoint weights.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        kpointw : ndarray
            An array containing the kpoint weights used in the 
            calculation.

        """

        entry = xml.findall(
            'kpoints/varray[@name="weights"]/v')
                
        kpointsw = self._convert_array1D_f(entry)

        return kpointsw

    def _fetch_kpointdiv(self, xml):
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

        entry = xml.find(
            'kpoints/generation/v[@name="divisions"]')

        kpointdiv = self._convert_array_i(entry)

        return kpointdiv

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
            data = np.fromstring(entry.text, sep=' ', dtype='intc')

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
            data = np.fromstring(entry.text, sep=' ', dtype='double')

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
            data = np.zeros(len(entry),dtype='intc')
        for index, element in enumerate(entry):
            data[index] = np.fromstring(element.text, sep=' ')

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
            data = np.zeros(len(entry),dtype='double')
        for index, element in enumerate(entry):
            data[index] = np.fromstring(element.text, sep=' ')

        return data

    def _convert_array2D_f(self, entry):
        """Convert the input entry to numpy array

        Parameters
        ----------
        entry : list
            A list containing Element objects where each
            element is a float

        Returns
        -------
        data : ndarray
            | Dimension: (N,3)
            An array containing N elements with three float
            elements.

        """

        data = None
        if entry is not None:
            data = np.zeros((len(entry), 3),
                                 dtype='double')

        for index, element in enumerate(entry):
            data[index] = np.fromstring(element.text, sep=' ')

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
                species[index] = constants.elements[
                    entry[index].text.split()[0].lower()]
            except KeyError:
                self._logger.error("There is an atomic element present in the "
                                   "XML file that is unknown. Exiting.")
                sys.exit(1)

        return species

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
