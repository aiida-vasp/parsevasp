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

        # extract data from all calculations (e.g. ionic steps)
        self._extract_all = True

        # kpoint index before band index (for instance for the ordering
        # of the eigenvalue data etc.)?
        self._k_before_band = False

        # dictionaries that contain the output of the parsing
        self._parameters = {"symprec": None,
                            "ismear": None,
                            "sigma": None,
                            "ispin": None,
                            "nbands": None}
        self._lattice = {"unitcell": None,
                         "species": None,
                         "positions": None,
                         "kpoints": None,
                         "kpointsw": None,
                         "kpointdiv": None}
        self._data = {"eigenvalues": None,
                      "occupancies": None,
                      "dos": None}
        self._calculations = {}

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

        # Check size of the XML file. For large files we need to
        # perform event driven parsing. For smaller files this is
        # not necessary and is too slow.
        if self._file_size(self._file_path) < self._sizecutoff:
            self._parsew()
            # self._parsee()
        else:
            self._parsee()

        print self._data["dos"]

    def _parsew(self):
        """Performs parsing on the whole XML files. For smaller files

        """

        self._logger.debug("Running parsew.")

        # now open the complete file
        self._check_file(self._file_path)
        vaspxml = etree.parse(self._file_path)

        # do we want to extract data from all calculations (e.g. ionic steps)
        all = self._extract_all

        # let us start to parse the content
        self._parameters["symprec"] = self._fetch_symprecw(vaspxml)
        self._parameters["sigma"] = self._fetch_sigmaw(vaspxml)
        self._parameters["ismear"] = self._fetch_ismearw(vaspxml)
        self._parameters["ispin"] = self._fetch_ispinw(vaspxml)
        self._parameters["nbands"] = self._fetch_nbandsw(vaspxml)
        self._lattice["unitcell"] = self._fetch_unitcellw(vaspxml, all=all)
        self._lattice["species"] = self._fetch_speciesw(vaspxml)
        self._lattice["positions"] = self._fetch_positionsw(vaspxml, all=all)
        self._lattice["kpoints"] = self._fetch_kpointsw(vaspxml)
        self._lattice["kpointsw"] = self._fetch_kpointsww(vaspxml)
        self._lattice["kpointdiv"] = self._fetch_kpointdivw(vaspxml)
        self._data["eigenvalues"], self._data[
            "occupancies"] = self._fetch_eigenvaluesw(vaspxml)
        self._data["dos"] = self._fetch_dosw(vaspxml)

    def _parsee(self):
        """Performs parsing in an event driven fashion on the XML file.
        Slower, but suitable for bigger files.

        """

        # set logger
        self._logger.debug("Running parsee.")

        # helper list
        data = []
        data2 = []

        # helper dicts
        cell = {}
        pos = {}

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
        extract_eigenvalues = False
        extract_eigenvalues_spin1 = False
        extract_eigenvalues_spin2 = False

        # do we want to extract data from all calculations (e.g. ionic steps)
        all = self._extract_all

        # index that control the calculation step (e.g. ionic step)
        calc = 0
        for event, element in etree.iterparse(self._file_path, events=("start", "end")):
            # set extraction points (what to read and when to read it)
            # here we also set the relevant data elements when the tags
            # close when they contain more than one element
            if event == "start" and element.tag == "parameters":
                extract_parameters = True
            if event == "end" and element.tag == "parameters":
                extract_parameters = False
            if all:
                if event == "start" and element.tag == "structure":
                    extract_latticedata = True
                if event == "end" and element.tag == "structure":
                    extract_latticedata = False
            else:
                if event == "start" and element.tag == "structure":
                    extract_latticedata = True
                if event == "end" and element.tag == "structure":
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
            if event == "start" and element.tag == "eigenvalues":
                extract_eigenvalues = True
            if event == "end" and element.tag == "eigenvalues":
                eigenvalues, occupancies = self._extract_eigenvalues(
                    data, data2)
                self._data["eigenvalues"] = eigenvalues
                self._data["occupancies"] = occupancies
                data = []
                data2 = []
                extract_eigenvalues = False

            # now fetch the data
            if extract_parameters:
                try:
                    if event == "start" and element.attrib["name"] == "SYMPREC":
                        self._parameters["symprec"] = self._convert_f(element)
                except KeyError:
                    pass
                try:
                    if event == "start" and element.attrib["name"] == "ISPIN":
                        self._parameters["ispin"] = self._convert_i(element)
                except KeyError:
                    pass
                try:
                    if event == "start" and element.attrib["name"] == "ISMEAR":
                        self._parameters["ismear"] = self._convert_i(element)
                except KeyError:
                    pass
                try:
                    if event == "start" and element.attrib["name"] == "SIGMA":
                        self._parameters["sigma"] = self._convert_f(element)
                except KeyError:
                    pass
                try:
                    if event == "start" and element.attrib["name"] == "NBANDS":
                        self._parameters["nbands"] = self._convert_i(element)
                except KeyError:
                    pass

            if extract_latticedata:
                # print event, element.tag, element.text, element.attrib
                if event == "start" and element.tag == "varray" \
                   and element.attrib.get("name") == "basis":
                    extract_unitcell = True
                if event == "end" and element.tag == "varray" \
                   and element.attrib.get("name") == "basis":
                    if calc == 0:
                        attribute = "initial"
                    else:
                        attribute = "step_" + str(calc)
                    cell[attribute] = self._convert_array2D3_f(data)
                    data = []
                    extract_unitcell = False
                if event == "start" and element.tag == "varray" \
                   and element.attrib.get("name") == "positions":
                    extract_positions = True
                if event == "end" and element.tag == "varray" \
                   and element.attrib.get("name") == "positions":
                    if calc == 0:
                        attribute = "initial"
                    else:
                        attribute = "step_" + str(calc)
                    pos[attribute] = self._convert_array2D3_f(data)
                    data = []
                    # if we do multiple calculations (e.g. ionic) the
                    # position data is always updated? we assume this for now
                    # otherwise this needs to be more fancy
                    calc = calc + 1
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
                        self._lattice[
                            "kpointdiv"] = self._convert_array_i(element)
                except KeyError:
                    pass

                if event == "start" and element.tag == "varray" \
                   and element.attrib.get("name") == "kpointlist":
                    extract_kpoints = True
                if event == "end" and element.tag == "varray" \
                   and element.attrib.get("name") == "kpointlist":
                    self._lattice["kpoints"] = self._convert_array2D3_f(data)
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

            if extract_eigenvalues:
                if event == "start" and element.tag == "set" \
                   and element.attrib.get("comment") == "spin 1":
                    extract_eigenvalues_spin1 = True
                if event == "end" and element.tag == "set" \
                   and element.attrib.get("comment") == "spin 1":
                    extract_eigenvalues_spin1 = False
                if event == "start" and element.tag == "set" \
                   and element.attrib.get("comment") == "spin 2":
                    extract_eigenvalues_spin2 = True
                if event == "end" and element.tag == "set" \
                   and element.attrib.get("comment") == "spin 2":
                    extract_eigenvalues_spin2 = False

                if extract_eigenvalues_spin1:
                    if event == "start" and element.tag == "r":
                        data.append(element)

                if extract_eigenvalues_spin2:
                    if event == "start" and element.tag == "r":
                        data2.append(element)

        # now we need to update some elements
        last_element = len(pos) - 1
        # the two last and two first elements should be the same,
        # so remove them
        del cell["step_1"]
        del pos["step_1"]
        if last_element > 2:
            # in cases where a static run is done, we will have
            # initial, step_1 and final so we only have to delete
            # step_1, which is done above, here we delete the next
            # last item for the other cases
            del cell["step_" + str(last_element - 1)]
            del pos["step_" + str(last_element - 1)]
        cell["final"] = cell.pop("step_" + str(last_element))
        pos["final"] = pos.pop("step_" + str(last_element))
        if not all:
            # only save initial and final
            self._lattice["unitcell"] = {key: cell[key]
                                         for key in {"initial", "final"}}
            self._lattice["positions"] = {key: cell[key]
                                          for key in {"initial", "final"}}
        else:
            # save all
            self._lattice["unitcell"] = cell
            self._lattice["positions"] = pos

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

        symprec = None
        entry = xml.find('.//parameters/separator[@name="symmetry"]/'
                         'i[@name="SYMPREC"]')
        if entry is not None:
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

        sigma = None
        entry = xml.find('.//parameters/separator[@name="electronic"]/'
                         'separator[@name="electronic smearing"]/'
                         'i[@name="SIGMA"]')

        if entry is not None:
            sigma = self._convert_f(entry)

        return sigma

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

        ispin = None
        entry = xml.find('.//parameters/separator[@name="electronic"]/'
                         'separator[@name="electronic spin"]/'
                         'i[@name="ISPIN"]')
        if entry is not None:
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

        ismear = None
        entry = xml.find('.//parameters/separator[@name="electronic"]/'
                         'separator[@name="electronic smearing"]/'
                         'i[@name="ISMEAR"]')
        if entry is not None:
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

        nbands = None
        entry = xml.find('.//parameters/separator[@name="electronic"]/'
                         'i[@name="NBANDS"]')
        if entry is not None:
            nbands = self._convert_i(entry)

        return nbands

    def _fetch_unitcellw(self, xml, all=False):
        """Fetch the unitcell

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        all : bool 
            Determines which unitcell to get. Defaults to the initial
            and final. If True, extract all.

        Returns
        -------
        cell : dict
            An dictionary containing ndarrays of the unitcell with 
            vectors as rows in units of AA for each calculation.

        """

        cell = {}
        if not all:
            entry = xml.findall(
                './/structure[@name="finalpos"]/crystal/varray[@name="basis"]/v')
            cell["final"] = self._convert_array2D3_f(entry)
            entry = xml.findall(
                './/structure[@name="initialpos"]/crystal/varray[@name="basis"]/v')
            cell["initial"] = self._convert_array2D3_f(entry)
        else:
            entry = xml.findall(
                './/calculation/structure/crystal/varray[@name="basis"]/v')
            entries = len(entry)
            num_calcs = entries / 3
            cell["initial"] = self._convert_array2D3_f(entry[0:3])
            cell["final"] = self._convert_array2D3_f(entry[-3:])
            for calc in range(1, num_calcs - 1):
                base = calc * 3
                cell[
                    "step_" + str(calc + 1)] = self._convert_array2D3_f(entry[base:base + 3])

        return cell

    def _fetch_positionsw(self, xml, all=False):
        """Fetch the atomic positions.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        all : bool
            Determines which atomic positions to get. Defaults to the
            initial and final. If True, extract all.

        Returns
        -------
        pos : dict
            An dictionary containing ndarrays of the positions with 
            each position as rows in direct coordinates for each calculation.

        """

        pos = {}
        if not all:
            entry = xml.findall(
                './/structure[@name="finalpos"]/varray[@name="positions"]/v')
            pos["final"] = self._convert_array2D3_f(entry)
            entry = xml.findall(
                './/structure[@name="initialpos"]/varray[@name="positions"]/v')
            pos["initial"] = self._convert_array2D3_f(entry)
        else:
            entry = xml.findall(
                './/calculation/structure/varray[@name="positions"]/v')
            entries = len(entry)
            num_atoms = self._lattice["species"].shape[0]
            num_calcs = entries / num_atoms
            pos["initial"] = self._convert_array2D3_f(entry[0:num_atoms])
            pos["final"] = self._convert_array2D3_f(entry[-num_atoms:])
            for calc in range(1, num_calcs - 1):
                base = calc * num_atoms
                pos["step_" + str(calc + 1)
                    ] = self._convert_array2D3_f(entry[base:base + num_atoms])

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

        kpoints = self._convert_array2D3_f(entry)

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

        entry = xml.find(
            'kpoints/generation/v[@name="divisions"]')

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
            An tupple of two ndarrays containing the eigenvalues
            for each spin, band and kpoint index.

        """

        # spin 1
        entry_ispin1 = xml.findall(
            './/calculation/eigenvalues/array/set/set[@comment="spin 1"]/set/r')

        # spin 2
        entry_ispin2 = xml.findall(
            './/calculation/eigenvalues/array/set/set[@comment="spin 2"]/set/r')

        eigenvalues, occupancies = self._extract_eigenvalues(entry_ispin1,
                                                             entry_ispin2)

        return eigenvalues, occupancies

    def _extract_eigenvalues(self, data1, data2):
        """Extract the eigenvalues.

        Parameters
        ----------
        data1 : list
            A list of ElementTree object to be used for parsing of the
            ispin=1 entries.
        data2 : list
            A list of ElementTree object to be used for parsing of the
            ispin=2 entries.

        Returns
        -------
        eigenvalues, occupancies : tupple
            An tupple of two ndarrays containing the eigenvalues
            for each spin, band and kpoint index.

        """

        # first check if we have extracted the kpoints
        if self._lattice["kpoints"] is None:
            self._logger.error("Before extracting the eigenvalues, please"
                               "extract the kpoints. Exiting.")
            sys.exit(1)

        # then check if we have asigned ispin
        if self._parameters["ispin"] is None:
            self._logger.error("Before extracting the eigenvalues, please"
                               "extract ISPIN. Exiting.")
            sys.exit(1)

        # then check if we have asigned nbands
        if self._parameters["nbands"] is None:
            self._logger.error("Before extracting the eigenvalues, please"
                               "extract NBANDS. Exiting.")
            sys.exit(1)

        # number of kpoints to disect the eigenvalue sets later
        num_kpoints = self._lattice["kpoints"].shape[0]

        # ispin
        ispin = self._parameters["ispin"]

        # number of bands
        num_bands = self._parameters["nbands"]

        data = []

        if len(data1) != num_bands * num_kpoints:
            self._logger.error("The number of eigenvalues found does not match "
                               "the number of located kpoints and NBANDS. "
                               "Exiting.")
            sys.exit(1)

        data.append([])

        if data2:
            if len(data2) != num_bands * num_kpoints:
                self._logger.error("The number of eigenvalues found does not match "
                                   "the number of located kpoints and NBANDS. "
                                   "Exiting.")
                sys.exit(1)
            data.append([])

        # a bit of a nasty and explicit loop, improve performance later
        for kpoint in range(num_kpoints):
            base = kpoint * num_bands
            data[0].append(
                self._convert_array2D2_f(data1[base:base + num_bands]))
            if data2:
                data[1].append(
                    self._convert_array2D2_f(data2[base:base + num_bands]))

        # convert to numpy array
        data = np.asarray(data)

        # swap axis if the band index should be before the kpoint index
        if not self._k_before_band:
            data = np.swapaxes(data, 1, 2)

        # set data and make sure it is continues in memory
        eigenvalues = np.ascontiguousarray(data[:, :, :, 0])
        occupancies = np.ascontiguousarray(data[:, :, :, 1])

        return eigenvalues, occupancies

    def _fetch_dosw(self, xml):
        """Fetch the density of states.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.

        Returns
        -------
        dos : dict
            A dict of ndarrays containing the energies, total and
            integrated density of states.

        """

        # spin 1
        entry_total_ispin1 = xml.findall(
            './/calculation/dos/total/array/set/set[@comment="spin 1"]/r')

        # spin 2
        entry_total_ispin2 = xml.findall(
            './/calculation/dos/total/array/set/set[@comment="spin 2"]/r')

        dos_ispin1 = self._convert_array2D3_f(entry_total_ispin1)

        if entry_total_ispin2:
            dos_ispin2 = self._convert_array2D3_f(entry_total_ispin2)
            _dos = np.concatenate((dos_ispin1, dos_ispin2), axis=0)
        else:
            _dos = dos_ispin1

        _dos = np.ascontiguousarray(_dos)

        # now make a dict
        dos = {}
        dos["energy"] = _dos[:, 0]
        dos["total"] = _dos[:, 1]
        dos["integrated"] = _dos[:, 2]

        return dos

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
            data = np.zeros(len(entry), dtype='intc')
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
            data = np.zeros(len(entry), dtype='double')
        for index, element in enumerate(entry):
            data[index] = np.fromstring(element.text, sep=' ')

        return data

    def _convert_array2D3_f(self, entry):
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

    def _convert_array2D2_f(self, entry):
        """Convert the input entry to numpy array

        Parameters
        ----------
        entry : list
            A list containing Element objects where each
            element is a float

        Returns
        -------
        data : ndarray
            | Dimension: (N,2)
            An array containing N elements with three float
            elements.

        """

        data = None
        if entry is not None:
            data = np.zeros((len(entry), 2),
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
