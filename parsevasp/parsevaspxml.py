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
                            "nbands": None,
                            "nelect": None,
                            "system": None}
        self._lattice = {"unitcell": None,
                         "species": None,
                         "positions": None,
                         "kpoints": None,
                         "kpointsw": None,
                         "kpointdiv": None}
        self._data = {"eigenvalues": None,
                      "occupancies": None,
                      "dos": None,
                      "totens": None,
                      "forces": None,
                      "stress": None,
                      "dielectrics": None}

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
            #self._parsew()
            self._parsee()
        else:
            self._parsee()
        #print self._data["dos"]
        print self._parameters
        print self._data["dielectrics"]
        
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
        self._parameters["nelect"] = self._fetch_nelectw(vaspxml)
        self._parameters["system"] = self._fetch_systemw(vaspxml)
        self._lattice["species"] = self._fetch_speciesw(vaspxml)
        self._lattice["unitcell"], self._lattice["positions"], \
            self._data["forces"], self._data["stress"] = \
                self._fetch_upfsw(vaspxml, all=all)
        self._lattice["kpoints"] = self._fetch_kpointsw(vaspxml)
        self._lattice["kpointsw"] = self._fetch_kpointsww(vaspxml)
        self._lattice["kpointdiv"] = self._fetch_kpointdivw(vaspxml)
        self._data["eigenvalues"], self._data[
            "occupancies"] = self._fetch_eigenvaluesw(vaspxml)
        self._data["dos"] = self._fetch_dosw(vaspxml)
        self._data["totens"] = self._fetch_totensw(vaspxml)
        self._data["dielectrics"] = self._fetch_dielectricsw(vaspxml)

    def _parsee(self):
        """Performs parsing in an event driven fashion on the XML file.
        Slower, but suitable for bigger files.

        """

        # set logger
        self._logger.debug("Running parsee.")

        # helper list
        data = []
        data2 = []
        data3 = []
        data4 = []

        # dicts
        cell = {}
        pos = {}
        force = {}
        stress = {}
        dos = {}
        totens = {}
        _dos = {}
        _dos2 = {}

        # bool to control extraction of content
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
        extract_eigenvalues = False
        extract_eigenvalues_spin1 = False
        extract_eigenvalues_spin2 = False
        extract_dos = False
        extract_total_dos = False
        extract_partial_dos = False
        extract_dos_ispin1 = False
        extract_dos_ispin2 = False
        extract_projected = False
        extract_forces = False
        extract_force = False
        extract_stress = False
        extract_stres = False
        extract_ewoe = False
        extract_scstep = False
        extract_dielectrics = False

        # do we want to extract data from all calculations (e.g. ionic steps)
        all = self._extract_all

        # index that control the calculation step (e.g. ionic step)
        calc = 1
        for event, element in etree.iterparse(self._file_path, events=("start", "end")):
            # set extraction points (what to read and when to read it)
            # here we also set the relevant data elements when the tags
            # close when they contain more than one element
            if event == "start" and element.tag == "parameters":
                extract_parameters = True
            if event == "end" and element.tag == "parameters":
                extract_parameters = False
            if event == "start" and element.tag == "calculation":
                extract_calculation = True
            if event == "end" and element.tag == "calculation":
                data3 = self._convert_array1D_f(data3)
                totens["step_"+str(calc)] = {"energy_wo_entropy":
                                             [data3[data3.shape[0]-1], data3]}
                data3 = []
                # update index for the calculation
                calc = calc + 1
                extract_calculation = False
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
            if event == "start" and element.tag == "projected":
                extract_projected = True
            if event == "end" and element.tag == "projected":
                extract_projected = False
            if event == "start" and element.tag == "dos":
                extract_dos = True
            if event == "end" and element.tag == "total":
                if data2:
                    dos_ispin = self._convert_array2D3_f(data)
                    _dos["energy"] = dos_ispin[:,0]
                    _dos["total"] = dos_ispin[:,1]
                    _dos["integrated"] = dos_ispin[:,2]
                    dos_ispin = self._convert_array2D3_f(data2)
                    _dos2["energy"] = dos_ispin[:,0]
                    _dos2["total"] = dos_ispin[:,1]
                    _dos2["integrated"] = dos_ispin[:,2]
                else:
                    dos_ispin = self._convert_array2D3_f(data)
                    _dos["energy"] = dos_ispin[:,0]
                    _dos["total"] = dos_ispin[:,1]
                    _dos["integrated"] = dos_ispin[:,2]
                data = []
                data2 = []
            if event == "end" and element.tag == "partial":
                num_atoms = 0
                if self._lattice["species"] is not None:
                    num_atoms = self._lattice["species"].shape[0]
                else:
                    self._logger.error("Before extracting the density of "
                                       "states, please extract the species. "
                                       "Exiting.")
                    sys.exit(1)
                if data2:
                    dos_ispin = self._convert_array2D10_f(data)
                    # do not need the energy term (similar to total)
                    _dos["partial"] = np.asarray(
                        np.split(dos_ispin[:,1:10], num_atoms))
                    dos_ispin = self._convert_array2D10_f(data2)
                    # do not need the energy term (similar to total)
                    _dos2["partial"] = np.asarray(
                        np.split(dos_ispin[:,1:10], num_atoms))
                else:
                    dos_ispin = self._convert_array2D10_f(data)
                    # do not need the energy term (similar to total)
                    _dos["partial"] = np.asarray(
                        np.split(dos_ispin[:,1:10], num_atoms))
                data = []
                data2 = []
            if event == "end" and element.tag == "dos":
                # check the Fermi level
                if len(data4) == 1:
                    fermi_level = self._convert_f(data4[0])
                elif len(data4) > 1:
                    self._logger.error("Multiple entries of efermi was located. "
                                       "Exiting.")
                    sys.exit(1)
                else:
                    fermi_level = None
                    
                if data2:
                    dos["up"] = _dos
                    dos["down"] = _dos2
                    dos["total"] = {"fermi_level": fermi_level}
                else:
                    _dos["fermi_level"] = fermi_level
                    dos["total"] = _dos
                self._data["dos"] = dos
                data = []
                data2 = []
                data4 = []
                extract_dos = False
                 
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
                try:
                    if event == "start" and element.attrib["name"] == "NELECT":
                        self._parameters["nelect"] = self._convert_f(element)
                except KeyError:
                    pass
                try:
                    if event == "start" and element.attrib["name"] == "SYSTEM":
                        self._parameters["system"] = element.text
                except KeyError:
                    pass
                

            if extract_calculation:
                attribute = "step_" + str(calc)
                # it would be very tempting just to fill the data and disect
                # it later, would be faster, but it is not so easy since
                # we do not know how many calculations have been performed
                # or how many scteps there are per calculation
                if event == "start" and element.tag == "structure":
                    extract_latticedata = True
                if event == "end" and element.tag == "structure":
                    extract_latticedata = False
                if event == "start" and element.tag == "varray" and \
                   element.attrib["name"] == "forces":
                    extract_forces = True
                if event == "end" and element.tag == "varray" and \
                   element.attrib["name"] == "forces":
                    extract_forces = False
                if event == "start" and element.tag == "varray" and \
                   element.attrib["name"] == "stress":
                    extract_stress = True
                if event == "end" and element.tag == "varray" and \
                   element.attrib["name"] == "stress":
                    extract_stress = False
                if event == "start" and element.tag == "scstep":
                    extract_scstep = True
                if event == "end" and element.tag == "scstep":
                    extract_scstep = False
                if event == "start" and element.tag == "eigenvalues":
                    # do not start regular eigenvalue extraction
                    # if projected dataset is found
                    # do not think we need tests like this on the end
                    # statements as the xml is always serial?
                    if not extract_projected:
                        extract_eigenvalues = True
                if event == "end" and element.tag == "eigenvalues" and not \
                   extract_projected:
                    # we do not do this for the projected eigenvalues,
                    # that needs special threatment
                    eigenvalues, occupancies = self._extract_eigenvalues(
                        data, data2)
                    self._data["eigenvalues"] = eigenvalues
                    self._data["occupancies"] = occupancies
                    data = []
                    data2 = []
                    extract_eigenvalues = False
                if event == "start" and element.tag == "dielectricfunction":
                    extract_dielectrics = True
                if event == "end" and element.tag == "dielectricfunction":
                    _diel = {}
                    diel = np.split(self._convert_array2D7_f(data), 2)
                    _diel["energies"] = diel[0][:,0]
                    _diel["imag"] = diel[0][:,1-7]
                    _diel["real"] = diel[1][:,1-7]
                    self._data["dielectrics"] = _diel
                    data = []
                    extract_dielectrics = False
                    
                # now extract data
                if extract_scstep:
                    #print element.tag, element.attrib, element.text
                    # energy without entropy
                    if event == "start" and element.tag == "i" and \
                       element.attrib["name"] == "e_wo_entrp":
                        extract_ewoe = True
                    if event == "end" and element.tag == "i" and \
                       element.attrib["name"] == "e_wo_entrp":
                        extract_ewoe = False
                    if extract_ewoe:
                        data3.append(element)
                if extract_latticedata:
                    if event == "start" and element.tag == "varray" \
                       and element.attrib.get("name") == "basis":
                        extract_unitcell = True
                    if event == "end" and element.tag == "varray" \
                       and element.attrib.get("name") == "basis":
                        cell[attribute] = self._convert_array2D3_f(data)
                        data = []
                        extract_unitcell = False

                    if event == "start" and element.tag == "varray" \
                       and element.attrib.get("name") == "positions":
                        extract_positions = True
                    if event == "end" and element.tag == "varray" \
                       and element.attrib.get("name") == "positions":
                        pos[attribute] = self._convert_array2D3_f(data)
                        data = []
                        extract_positions = False

                    if extract_unitcell:
                        if event == "start" and element.tag == "v":
                            data.append(element)
                    if extract_positions:
                        if event == "start" and element.tag == "v":
                            data.append(element)
                if extract_forces:
                    if event == "start" and element.tag == "v":
                        extract_force = True
                    if event == "end" and element.tag == "v":
                        force[attribute] = self._convert_array2D3_f(data)
                        data = []
                        extract_force = False
                    if extract_force:
                        data.append(element)

                if extract_stress:
                    if event == "start" and element.tag == "v":
                        extract_stres = True
                    if event == "end" and element.tag == "v":
                        stress[attribute] = self._convert_array2D3_f(data)
                        data = []
                        extract_stres = False
                    if extract_stres:
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

                if extract_dielectrics:
                    if event == "start" and element.tag == "r":
                        data.append(element)
                        
            if extract_species:
                if event == "start" and element.tag == "c":
                    data.append(element)

            if extract_kpointdata:
                try:
                    if event == "start" and element.tag == "v" and \
                       element.attrib["name"] == "divisions":
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

            if extract_dos:
                try:
                    if event == "start" and element.tag == "i" and \
                       element.attrib["name"] == "efermi":
                        data4.append(element)
                except KeyError:
                    pass

                if event == "start" and element.tag == "set" \
                   and element.attrib.get("comment") == "spin 1":
                    extract_dos_ispin1 = True
                if event == "end" and element.tag == "set" \
                   and element.attrib.get("comment") == "spin 1":
                    extract_dos_ispin1 = False
                if event == "start" and element.tag == "set" \
                   and element.attrib.get("comment") == "spin 2":
                    extract_dos_ispin2 = True
                if event == "end" and element.tag == "set" \
                   and element.attrib.get("comment") == "spin 2":
                    extract_dos_ispin2 = False
                if extract_dos_ispin1:
                    if event == "start" and element.tag == "r":
                        data.append(element)
                if extract_dos_ispin2:
                    if event == "start" and element.tag == "r":
                        data2.append(element)
                        
            if extract_force:
                if event == "start" and element.tag == "v":
                    data.append(element)

            if extract_stress:
                if event == "start" and element.tag == "v":
                    data.append(element)

        # now we need to update some elements
        # first element should be initial
        cell["initial"] = cell.pop("step_1")
        pos["initial"] = pos.pop("step_1")
        force["initial"] = force.pop("step_1")
        stress["initial"] = stress.pop("step_1")
        totens["initial"] = totens.pop("step_1")
        if len(cell) == 1:
            # for static runs, initial is equal to final
            cell["final"] = cell["initial"]
            pos["final"] = pos["initial"]
            force["final"] = force["initial"]
            stress["final"] = stress["initial"]
            totens["final"] = totens["initial"]
        else:
            last_element = "step_"+str(len(cell)-1)
            cell["final"] = cell.pop(last_element)
            pos["final"] = pos.pop(last_element)
            force["final"] = force.pop(last_element)
            stress["final"] = stress.pop(last_element)
            totens["final"] = totens.pop(last_element)

        if not all:
            # only save initial and final
            self._lattice["unitcell"] = {key: cell[key]
                                         for key in {"initial", "final"}}
            self._lattice["positions"] = {key: cell[key]
                                          for key in {"initial", "final"}}
            self._data["forces"] = {key: force[key]
                                         for key in {"initial", "final"}}
            self._data["stress"] = {key: stress[key]
                                         for key in {"initial", "final"}}
            self._data["totens"] = {key: totens[key]
                                         for key in {"initial", "final"}}
                        
        else:
            # save all
            self._lattice["unitcell"] = cell
            self._lattice["positions"] = pos
            self._data["forces"] = force
            self._data["stress"] = stress
            self._data["totens"] = totens

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

        nelect = None
        entry = xml.find('.//parameters/separator[@name="electronic"]/'
                         'i[@name="NELECT"]')
        if entry is not None:
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

        system = None
        entry = xml.find('.//parameters/separator[@name="general"]/'
                         'i[@name="SYSTEM"]')
        if entry is not None:
            system = entry.text
            
        return system
    
    def _fetch_upfsw(self, xml, all=False):
        """Fetch the unitcell, atomic positions, force and stress.

        Parameters
        ----------
        xml : object
            An ElementTree object to be used for parsing.
        all : bool 
            Determines which unitcell and positions to get. 
            Defaults to the initial and final. If True, extract all.

        Returns
        -------
        cell, pos, force, stress : dict
            An dictionary containing ndarrays of the:
            | unitcells with the vectors as rows in AA.

            | positions with each position as a row in direct coordinates.

            | forces where each row is the force in eVAA on each atom.

            | stress where each row is the stress vector for the unitcell.

        """

        cell = {}
        pos = {}
        force = {}
        stress = {}
        if not all:
            entry = xml.findall(
                './/structure[@name="finalpos"]/crystal/varray[@name="basis"]/v')
            cell["final"] = self._convert_array2D3_f(entry)
            entry = xml.findall(
                './/structure[@name="initialpos"]/crystal/varray[@name="basis"]/v')
            cell["initial"] = self._convert_array2D3_f(entry)
            entry = xml.findall(
                './/structure[@name="finalpos"]/varray[@name="positions"]/v')
            pos["final"] = self._convert_array2D3_f(entry)
            entry = xml.findall(
                './/structure[@name="initialpos"]/varray[@name="positions"]/v')
            pos["initial"] = self._convert_array2D3_f(entry)
        else:
            num_atoms = 0
            if self._lattice["species"] is not None:
                num_atoms = self._lattice["species"].shape[0]
            else:
                self._logger.error("Before extracting the unitcell and"
                                   "positions please extract the species. "
                                   "Exiting.")
                sys.exit(1)
        
            entrycell = xml.findall(
                './/calculation/structure/crystal/varray[@name="basis"]/v')
            entrypos = xml.findall(
                './/calculation/structure/varray[@name="positions"]/v')
            entryforce = xml.findall(
                './/calculation/varray[@name="forces"]/v')
            entrystress = xml.findall(
                './/calculation/varray[@name="stress"]/v')        
            entries = len(entrycell)
            num_calcs = entries / 3
            cell["initial"] = self._convert_array2D3_f(entrycell[0:3])
            cell["final"] = self._convert_array2D3_f(entrycell[-3:])
            pos["initial"] = self._convert_array2D3_f(entrypos[0:num_atoms])
            pos["final"] = self._convert_array2D3_f(entrypos[-num_atoms:])
            force["initial"] = self._convert_array2D3_f(entryforce[0:num_atoms])
            force["final"] = self._convert_array2D3_f(entryforce[-num_atoms:])
            stress["initial"] = self._convert_array2D3_f(entrystress[0:3])
            stress["final"] = self._convert_array2D3_f(entrystress[-3:])
            for calc in range(1, num_calcs - 1):
                basecell = calc * 3
                basepos = calc * num_atoms
                cell[
                    "step_" + str(calc + 1)] = self._convert_array2D3_f(
                        entrycell[basecell:basecell + 3])
                pos[
                    "step_" + str(calc + 1)] = self._convert_array2D3_f(
                        entrypos[basepos:basepos + num_atoms])
                force[
                    "step_" + str(calc + 1)] = self._convert_array2D3_f(
                        entryforce[basepos:basepos + num_atoms])
                stress[
                    "step_" + str(calc + 1)] = self._convert_array2D3_f(
                        entrystress[basecell:basecell + 3])
                
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
        entries = xml.findall(
            './/calculation')

        # this most likely takes too long for very long fpmd calculations,
        # so consider putting in a flag that only extract the
        # energies from each step in the calculation and not the scsteps as
        # well
        energies = {}
        for index, calc in enumerate(entries):
            # energy without entropy
            data = calc.findall('.//scstep/energy/i[@name="e_wo_entrp"]')
            data = self._convert_array1D_f(data)
            if index == 0:
                energies["initial"] = {"energy_no_entropy":
                                       [data[data.shape[0]-1], data]}
            else:
                energies["step_"+str(index+1)] = {"energy_no_entropy":
                                                  [data[data.shape[0]-1], data]}
                
        if len(energies) == 1:
            energies["final"] = energies["initial"]
        else:
            # replace key on the last element to final
            energies["final"] = energies.pop("step_"+str(len(entries)))

        return energies        

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

        # fetch the Fermi level
        entry = xml.find(
            './/calculation/dos/i[@name="efermi"]')

        if entry is not None:
            fermi_level = self._convert_f(entry)
        else:
            fermi_level = None

        # spin 1
        entry_total_ispin1 = xml.findall(
            './/calculation/dos/total/array/set/set[@comment="spin 1"]/r')

        # spin 2
        entry_total_ispin2 = xml.findall(
            './/calculation/dos/total/array/set/set[@comment="spin 2"]/r')

        # partial spin 1
        entry_partial_ispin1 = xml.findall(
            './/calculation/dos/partial/array/set/set/set[@comment="spin 1"]/r')

        # partial spin 2
        entry_partial_ispin2 = xml.findall(
            './/calculation/dos/partial/array/set/set/set[@comment="spin 2"]/r')

        # check if we have extracted the species (number of atoms)
        if self._lattice["species"] is None:
            self._logger.error("Before extracting the density of states, please"
                               "extract the species. Exiting.")
            sys.exit(1)

        # number of atoms
        num_atoms = self._lattice["species"].shape[0]
        
        if entry_total_ispin2:
            dos = {}
            dos = {"up": None, "down": None}
            dos_ispin = self._convert_array2D3_f(entry_total_ispin1)
            _dos = {}
            _dos["energy"] = dos_ispin[:,0]
            _dos["total"] = dos_ispin[:,1]
            _dos["integrated"] = dos_ispin[:,2]
            # check if partial exists
            if entry_partial_ispin1:
                dos_ispin = self._convert_array2D10_f(entry_partial_ispin1)
                # do not need the energy term (similar to total)
                _dos["partial"] = np.asarray(np.split(dos_ispin[:,1:10],num_atoms))
            else:
                _dos["partial"] = None
            dos["up"] = _dos
            dos_ispin = self._convert_array2D3_f(entry_total_ispin2)
            _dos["energy"] = dos_ispin[:,0]
            _dos["total"] = dos_ispin[:,1]
            _dos["integrated"] = dos_ispin[:,2]
            if entry_partial_ispin2:
                dos_ispin = self._convert_array2D10_f(entry_partial_ispin2)
                # do not need the energy term (similar to total)
                _dos["partial"] = np.asarray(np.split(dos_ispin[:,1:10],num_atoms))
            dos["down"] = _dos
            dos["total"] = {"fermi_level": fermi_level}
        else:
            dos = {}
            dos = {"total": None}
            dos_ispin = self._convert_array2D3_f(entry_total_ispin1)
            _dos = {}
            _dos["energy"] = dos_ispin[:,0]
            _dos["total"] = dos_ispin[:,1]
            _dos["integrated"] = dos_ispin[:,2]
            # check if partial exists
            if entry_partial_ispin1:
                dos_ispin = self._convert_array2D10_f(entry_partial_ispin1)
                # do not need the energy term (similar to total)
                _dos["partial"] = np.asarray(np.split(dos_ispin[:,1:10],num_atoms))
            else:
                _dos["partial"] = None
            _dos["fermi_level"] = fermi_level
            dos["total"] = _dos
        
        return dos

    def _fetch_dielectricsw(self, xml, method="dft", transfer=None):
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

        """
    
        if method == "dft":
            diel = {}
            if transfer == "density":
                tag = 'dielectricfunction[@comment="density-density"]'
            elif transfer == "current":
                tag = 'dielectricfunction[@comment="current-current"]'
            else:
                tag = 'dielectricfunction'

            # imaginary part
            entry = xml.findall(
                './/calculation/'+tag+'/imag/array/set/r')
            if not entry:
                self._logger.error("Did not find <dielectricfunction> in "
                                   "the XML file. Exiting.")
                sys.exit(1)
            data = self._convert_array2D7_f(entry)
            diel["energy"] = data[:,0]
            diel["imag"] = data[:,1:7]

            # real part
            entry = xml.findall(
                './/calculation/'+tag+'/real/array/set/r')
            data = self._convert_array2D7_f(entry)
            diel["real"] = data[:,1:7]

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
        eigenvalues, occupancies : tupple of dicts
            An tupple of two dicts containing ndarrays with the eigenvalues
            and occupancies for each band and kpoint index.

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

        # set dicts
        eigenvalues = {}
        occupancies = {}
        if data2:
            eigenvalues["up"] = np.ascontiguousarray(data[0, :, :, 0])
            occupancies["up"] = np.ascontiguousarray(data[0, :, :, 1])
            eigenvalues["down"] = np.ascontiguousarray(data[1, :, :, 0])
            occupancies["down"] = np.ascontiguousarray(data[1, :, :, 1])
        else:
            eigenvalues["total"] = np.ascontiguousarray(data[0, :, :, 0])
            occupancies["total"] = np.ascontiguousarray(data[0, :, :, 1])

        return eigenvalues, occupancies
        
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
    
    def _convert_array2D10_f(self, entry):
        """Convert the input entry to numpy array

        Parameters
        ----------
        entry : list
            A list containing Element objects where each
            element is a float

        Returns
        -------
        data : ndarray
            | Dimension: (N,10)
            An array containing N elements with ten float
            elements.

        """

        data = None
        if entry is not None:
            data = np.zeros((len(entry), 10),
                            dtype='double')

        for index, element in enumerate(entry):
            data[index] = np.fromstring(element.text, sep=' ')

        return data

    def _convert_array2D7_f(self, entry):
        """Convert the input entry to numpy array

        Parameters
        ----------
        entry : list
            A list containing Element objects where each
            element is a float

        Returns
        -------
        data : ndarray
            | Dimension: (N,7)
            An array containing N elements with ten float
            elements.

        """

        data = None
        if entry is not None:
            data = np.zeros((len(entry), 7),
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
