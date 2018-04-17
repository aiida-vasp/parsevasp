#!/usr/bin/python
import sys
import logging
import numpy as np
from collections import Counter
import StringIO

import utils


class Poscar(object):

    def __init__(self, poscar_string=None, poscar_dict=None,
                 file_path=None, logger=None):
        """Initialize a POSCAR object and set content as a dictionary.

        Parameters
        ----------
        poscar_string : string
            A string containing POSCAR entries. Must contain line
            breaks if multiline, otherwise the POSCAR will be mangled.
        poscar_dict : dict
            A dictionary containing the POSCAR entries.
        file_path : string
            The file path in which the POSCAR is read.
        logger : object, optional
            A standard Python logger object.


        """

        self.file_path = file_path
        self.poscar_dict = poscar_dict
        self.poscar_string = poscar_string

        # set logger
        if logger is not None:
            self._logger = logger
        else:
            logging.basicConfig(level=logging.DEBUG)
            self._logger = logging.getLogger('PoscarParser')

        # check that only one argument is supplied
        if (poscar_string is not None and poscar_dict is not None) \
           or (poscar_string is not None and file_path is not None) \
           or (poscar_dict is not None and file_path is not None):
            self._logger.error("Please only supply one argument when "
                              "initializing Poscar. Exiting.")
            sys.exit(1)
        # check that at least one is suplpied
        if (poscar_string is None and poscar_dict is None
                and file_path is None):
            self._logger.error("Please supply one argument when "
                              "initializing Poscar. Exiting.")
            sys.exit(1)

        if file_path is not None:
            # create dictionary from a file
            poscar_dict = self._from_file()

        if poscar_string is not None:
            # create dictionary from a string
            poscar_dict = self._from_string()

        # store entries
        self.entries = poscar_dict

        # validate dictionary
        self._validate()


    def _from_file(self):
        """Create rudimentary dictionary of entries from a
        file.

        """

        poscar = utils.readlines_from_file(self.file_path)
        poscar_dict = self._from_list(poscar)
        return poscar_dict

    def _from_string(self):
        """Create rudimentary dictionary of entries from a
        string.

        """

        poscar = self.poscar_string.splitlines()
        poscar_dict = self._from_list(poscar)
        return poscar_dict

    def _from_list(self, poscar):
        """Go through the list and analyze for = and ; in order to
        deentangle grouped entries etc.

        Parameters
        ----------
        poscar : list
            A list of strings containing each line in the POSCAR file.

        Returns
        -------
        poscar_dict : dictionary
            A dictionary containing each POSCAR tag as a key with the
            associated element.

        Notes
        -----
        No checking for consistency is done here. We do this at a later step
        in order to be able to keep the input methods as clean as posible.

        """

        comment = poscar[0]
        vasp5 = True
        # check for VASP 5 POSCAR
        if (utils.is_number(poscar[5])):
            vasp5 = False
        # set selective, test is done later
        selective = False
        # check scaling factor
        scaling = float(poscar[1].split()[0])
        if (scaling < 0.0):
            self._logger.error("Currently negative scaling values in POSCAR is "
                              "not supported. Exiting.")
            sys.exit(1)
        unitcell = [[0.0 for y in range(3)] for x in range(3)]
        nions = 0
        spec = None
        loopmax = 8
        if (vasp5):
            unitcell[0] = [float(x) for x in poscar[2].split()]
            unitcell[1] = [float(x) for x in poscar[3].split()]
            unitcell[2] = [float(x) for x in poscar[4].split()]
            spec = poscar[5].split()
            atoms = [int(x) for x in poscar[6].split()]
            for i in range(len(atoms)):
                nions = nions + atoms[i]
            if poscar[7][0].lower() == "s":
                selective = True
                loopmax = 9
                if not poscar[8][0].lower() == "d":
                    self._logger.error(
                        "Please supply a POSCAR in direct coordinates. "
                        "Exiting.")
                    sys.exit(1)
            if not selective:
                if not poscar[7][0].lower() == "d":
                    self._logger.error(
                        "Please supply a POSCAR in direct coordinates. "
                        "Exiting.")
                    sys.exit(1)
        else:
            self._logger.error(
                "VASP 4 POSCAR is not supported. User, please modernize. "
                "Exiting.")
            sys.exit(1)

        # create site objects
        specie_slot = 0
        index = 1
        sites_temp = []
        velocities = False
        predictor = False
        # loop positions
        for i in range(nions):
            # fetch specie
            if index > atoms[specie_slot]:
                specie_slot = specie_slot + 1
                index = 1
            specie = spec[specie_slot]
            # fetch positions
            line = poscar[i + loopmax].split()
            position = np.zeros(3)
            position[0] = float(line[0])
            position[1] = float(line[1])
            position[2] = float(line[2])
            # fetch selective flags
            flags = [True, True, True]
            if selective:
                if "f" in line[3].lower():
                    flags[0] = False
                if "f" in line[4].lower():
                    flags[1] = False
                if "f" in line[5].lower():
                    flags[2] = False
            # create a site object and add to sites list
            index = index + 1
            velo = None
            pred = None
            sites_temp.append([specie, position, flags, velo, pred])
        # now check if there is more in the POSCAR
        loopmax_pos = nions+loopmax
        if len(poscar) > loopmax_pos:
            first_char = poscar[loopmax_pos][0].lower()
            if first_char == 'd' \
               or first_char == 'c':
                velocities = True
            elif poscar[loopmax_pos] == '\n':
                predictor = True
        # now check that the next line is in fact a coordinate
        loopmax_pos = loopmax_pos + 1
        # allow for blank lines at the end of the positions
        if len(poscar) > loopmax_pos:
            if not utils.is_number(poscar[loopmax_pos].split()[0]):
                self._logger.error("A velocity or predictor-corrector "
                                  "coordinate was not detected. Exiting.")
                sys.exit(1)
        else:
            # but make sure the predictor is set back to False
            # if we only have a blank line and nothing else following
            # the coordinates
            predictor = False
        if velocities:
            for i in range(nions):
                # fetch velocities
                line = poscar[i + loopmax_pos].split()
                vel = np.zeros(3)
                vel[0] = float(line[0])
                vel[1] = float(line[1])
                vel[2] = float(line[2])
                sites_temp[i][3] = vel
        # now check if there is predictor-corrector coordinates following
        # the velocities
        loopmax_pos = nions + loopmax_pos
        if len(poscar) > loopmax_pos:
            if poscar[loopmax_pos] == '\n':
                loopmax_pos = loopmax_pos + 1
                if utils.is_number(poscar[loopmax_pos].split()[0]):
                    for i in range(nions):
                        line = poscar[i + loopmax_pos].split()
                        pre = np.zeros(3)
                        pre[0] = float(line[0])
                        pre[1] = float(line[1])
                        pre[2] = float(line[2])
                        sites_temp[i][4] = pre
        # do one final loop to create the objects and read
        # predictors if they exist
        sites = []
        loopmax_pos = nions + loopmax + 1
        for i in range(nions):
            pre = np.zeros(3)
            if predictor:
                line = poscar[i + loopmax_pos].split()
                pre[0] = float(line[0])
                pre[1] = float(line[1])
                pre[2] = float(line[2])
                sites_temp[i][4] = pre
            site = Site(sites_temp[i][0], sites_temp[i][1],
                        selective = sites_temp[i][2],
                        velocities = sites_temp[i][3],
                        predictors = sites_temp[i][4])
            sites.append(site)
            
        # apply scaling factor
        unitcell = [[x * scaling for x in y] for y in unitcell]
        # build dictionary and convert to NumPy
        poscar_dict = {}
        poscar_dict["comment"] = comment.strip()
        poscar_dict["unitcell"] = np.asarray(unitcell)
        poscar_dict["sites"] = sites
        return poscar_dict

    def modify(self, entry, value, site_number=None):
        """Modify an entry tag in the Poscar dictionary.
        If it is not found add it.

        Parameters
        ----------
        tag : string
            The entry tag of the POSCAR entry.
        value
            The entry value of the POSCAR entry.
            Can be either a string for the comment,
            an ndarray for the unitcell or a Site object
            for a position.
        site_number : int, optional
            The site to be modified. If not supplied
            the value have to be a list of Site objects.

        """

        # check allowed entries
        self._check_allowed_entries(entry)
        # check that entries exists
        self._check_dict()
        if site_number is not None:
            # check that a Site() object is supplied
            self._check_site(value)
            # check site number
            self._check_site_number(site_number)
            # check that position is an integer
            if not utils.is_number(site_number):
                self._logger.error("The supplied 'site_number' is "
                                  "not a number (i.e. the index) "
                                  "starting from 1 for the site "
                                  "position to be modified. Exiting.")
                sys.exit(1)
            self.entries["sites"][site_number] = value
        else:
            if entry == "sites":
                self._check_sites(sites = value)
            if entry == "comment":
                self._check_comment(comment = value)
            if entry == "unitcell":
                self._check_unitcell(unitcell = value)
            
            self.entries[entry] = value

    def delete_site(self, site_number):
        """Delete the site with the supplied
        number.

        Parameters
        ----------
        site_number : int
            The site number to be deleted, starting
            from 0.

        """

        self._check_sites()
        self._check_site_number(site_number)
        del self.entries["sites"][site_number]

    def add_site(self, site_number):
        """Add a site with the supplied
        number. If not supplied, add at the end
        of the last element of the specie group

        Parameters
        ----------
        site_number : int
            The site number to be deleted, starting
            from 0.

        """

        # EFL: ADD LATER

    def _check_dict(self):
        """Check that entries is present.

        """

        try:
            entries = self.entries
        except AttributeError:
            self._logger.error("There is no 'entries'. Exiting.")
            sys.exit(1)
        

    def _check_allowed_entries(self, entry):
        """Check the allowed values of entry.

        Parameters
        ----------
        entry : string
            Contains the entry to be checked.
        
        """

        if not (("comment" in entry) or ("unitcell" in entry) or
                ("sites" in entry)):
            self._logger.error("Only 'comment', 'unitcell' or "
                              "'sites' is allowed as input for "
                              "entry. Exiting.")
            sys.exit(1)

    def _check_unitcell(self, unitcell = None):
        """Check that the unitcell entries are present and
        are of a 3x3 ndarray type.

        Parameters
        ----------
        unitcell, optional
            The unitcell to be checked. If not supplied the
            'unitcell' key in the 'entries' is checked.

        """

        if unitcell is None:
            try:
                unitcell = self.entries["unitcell"]
            except KeyError:
                self._logger.error("There is no key 'unitcell' in "
                                  "'entries'. Exiting.")
                sys.exit(1)
            
        if (not isinstance(unitcell, np.ndarray)) \
           or (unitcell.shape != (3, 3)):
            self._logger.error("The value of 'unitcell' is "
                              "not an 3x3 ndarray. Exiting.")
            sys.exit(1)

    def _check_comment(self, comment = None):
        """Check that the comment entry is present and
        is a string.
        
        Parameters
        ----------
        comment, optional
            The comment to be checked. If not supplied the
            'comment' key in the 'entries' is checked.

        """
        if comment is None:
            try:
                comment = self.entries["comment"]
            except KeyError:
                self._logger.error("There is no key 'comment' in "
                                  "'entries'. Exiting.")
                sys.exit(1)
        # allow None for comment
        if self.entries["comment"] is not None:
            if not isinstance(comment, str):
                self._logger.error("The 'comment' is not a string. Exiting.")
                sys.exit(1)

    def _check_sites(self, sites = None):
        """Check that the sites entries are present.

        Parameters
        ----------
        sites, optional
            The sites to be checked. If not supplied the
            'sites' key in the 'entries' is checked.

        """
        
        if sites is None:
            try:
                sites = self.entries["sites"]
            except KeyError:
                self._logger.error("There is no key 'sites' in "
                                  "'entries'. Exiting.")
                sys.exit(1)

        if not isinstance(sites, list):
            self._logger.error("The 'sites' entry "
                              "have to be a list. Exiting.")
            sys.exit(1)
        
    def _check_site(self, site = None):
        """Check that the site entry is a Site() object.

        Parameters
        ----------
        site, optional
            The site to be checked. If not supplied the entries under
            the 'sites' key in the 'entries' is checked.

        """
        if site is None:
            try:
                sites = self.entries["sites"]
            except KeyError:
                self._logger.error("There is no key 'sites' in "
                                  "'entries'. Exiting.")
                sys.exit(1)

            for site in sites:
                if not isinstance(site, Site):
                    self._logger.error("At least one of the values in 'sites' "
                                      "is not a Site() object. Exiting.")
                    sys.exit(1)
        else:
            if not isinstance(site, Site):
                self._logger.error("The 'site' "
                                  "is not a Site() object. Exiting.")
                sys.exit(1)
                
    def _check_site_number(self, site_number):
        """Check that the site_number is an int and that
        it is not out of bounds.

        Parameters
        ----------
        site_number : int
            The site_number to be checked

        """
        
        if not isinstance(site_number, int):
            self._logger.error("The supplied 'site_number' is "
                              "not an integer. Exiting.")
            sys.exit(1)
        sites = self.entries["sites"]
        if site_number > (len(sites)-1):
            self._logger.error("The supplied site_number is larger "
                              "than the number of sites. Exiting.")
            sys.exit(1)

    def _validate(self):
        """Validate the content of entries

        """

        self._check_dict()
        self._check_comment()
        self._check_unitcell()
        self._check_sites()

    def _sort_and_group_sites(self):
        """Sort and group the positions and species to
        VASP specifications.

        Returns
        -------
        sites : list
            Contains site info for each site. Each site element
            contains a string describing the specie, a ndarray
            of floats describing the position, a ndarray of booleans
            to describe the selective flags and a boolean that 
            contains a flag that is True if the positions are in direct 
            coordinates.
        species : list of strings
            Contains the number of unique species
        num_species : list of ints
            Contains the occurancy of each specie in the same order as
            'species'.
        selective : bool
            True if any selective flags are enabled, False otherwise.

        """

        sites = []
        species = []
        selective = False
        velocities = False
        predictors = False
        for site in self.entries["sites"]:
            specie = site.get_specie()
            select = site.get_selective()
            position = site.get_position()
            direct = site.get_direct()
            vel = site.get_velocities()
            pre = site.get_predictors()
            if direct is False:
                # make sure it is direct
                position = self._to_direct(position)
                self._logger.error("Cartesian is not yet implemented. "
                                  "Exiting.")
                sys.exit(1)
            if False in select:
                selective = True
            if vel is not None:
                velocities = True
            if pre is not None:
                predictors = True
            sites.append([specie, position, select, direct,
                          vel, pre])
            species.append(specie)
            
        # find unique entries and their number
        counter = Counter(species)
        # Counter does not order, so order now with the
        # least occuring element first (typical for compounds)
        sorted_keys = sorted(counter, key=counter.get)        
        species = []
        num_species = []
        for key in sorted_keys:
            species.append(key)
            num_species.append(counter[key])

        # now make sure the sites is on the same order
        ordered_sites = []
        for specie in species:
            ordered_sites.extend([site for site in sites if specie == site[0]])
        # EFL: consider to also sort on coordinate after specie
        
        return ordered_sites, species, num_species, selective, \
            velocities, predictors
            
    def _get_key(self, item):
        """Key fetcher for the sorted function.

        """

        return item[0]

    def _to_direct(self, position):
        """Transforms the position from cartesian to direct
        coordinates.

        Parameters
        ----------
        position : ndarray
            An ndarray containing the position in cartesian coordinates.

        Returns
        -------
        position_cart : ndarray
            An ndarray containing the position in direct coordinates.

        """

        return position

    def _to_cart(self, position):
        """Transforms the position from direct to cartesian
        coordinates.

        Parameters
        ----------
        position : ndarray
            An ndarray containing the position in direct coordinates.

        Returns
        -------
        position_dir : ndarray
            An ndarray containing the position in cartesian coordinates.

        """

        return position
        
    def get(self, tag):
        """Return the value and comment of the entry with tag.

        Parameters
        ----------
        tag : string
            The entry tag of the POSCAR entry.        

        Returns
        -------
        value : string, int, float or list
            The value of the tag entry

        """

        value = None
        try:
            value = self.entries[tag].get_value()
        except KeyError:
            pass

        return value

    def get_dict(self):
        """Get a true dictionary containing the entries in an
        POSCAR compatible fashion.

        Returns
        -------
        poscar_dict : dict
            A dictionary on POSCAR compatible form.

        """
        
        dictionary = {}
        for key, entry in self.entries.iteritems():
            if key == 'sites':
                dictionary[key] = [[element.get_specie().capitalize(),
                                    element.get_position(),
                                    element.get_selective(),
                                    element.get_velocities(),
                                    element.get_predictors(),
                                    element.get_direct()]  for element in entry]
            else:
                dictionary[key] = entry

        return dictionary

    def get_string(self):
        """Get a string containing the entries in a POSCAR
        compatible fashion. Each line is broken by a newline
        character

        Returns
        -------

        poscar_string : str
            A string containing the POSCAR entries of the
            current instance.

        """
        
        string_object = StringIO.StringIO()
        self._write(poscar = string_object)
        poscar_string = string_object.getvalue()
        string_object.close()

        return poscar_string
        
    def write(self, file_path):
        """ Write POSCAR like files

        Parameters
        ----------
        file_path : str
            The location and filename of the POSCAR like file to be
            written.

        """

        poscar = utils.file_handler(file_path, status='w')
        self._write(poscar = poscar)
        utils.file_handler(file_handler=poscar)        
        
    def _write(self, poscar):
        """ Write POSCAR like files to a file or string

        Parameters
        ----------
        poscar : object
            Either a file object or a StringIO object.

        """

        self._validate()
        entries = self.entries
        comment = entries["comment"]
        unitcell = entries["unitcell"]
        # sort and group to VASP specifications
        sites, species, num_species, selective, velocities, predictors = \
            self._sort_and_group_sites()
        # update comment
        compound = ""
        for index, specie in enumerate(species):
            if num_species[index] == 1:
                num_string = ""
            else:
                num_string = str(num_species[index])
            compound = compound + str(specie).capitalize()+num_string
        if comment is None:
            comment = "# Compound: " + compound + "."
        else:
            comment = "# Compound: " + compound + \
                      ". Old comment: " + comment + "."
        poscar.write(comment + "\n")
        # we avoid usage of the scaling factor
        poscar.write("1.0\n")
        # write unitcell
        for i in range(3):
            poscar.write(str(unitcell[i][0]) + " " +
                         str(unitcell[i][1]) + " " +
                         str(unitcell[i][2]) + "\n")
        for specie in species:
            poscar.write(str(specie).capitalize() + " ")
        poscar.write("\n")
        # write number of species
        for number in num_species:
            poscar.write(str(number) + " ")
        poscar.write("\n")
        # write selective if any flags are True
        if selective:
            poscar.write("Selective dynamics\n")
        # always write direct
        poscar.write("Direct\n")
        # write positions
        for site in sites:
            poscar.write(str(site[1][0]) + " " +
                         str(site[1][1]) + " " + str(site[1][2]))
            if selective:
                sel = ["T", "T", "T"]
                flags = site[2]
                for index, flag in enumerate(flags):
                    if not flag:
                        sel[index] = "F"
                    
                poscar.write(" " + sel[0] + " " +
                             sel[1] + " " +
                             sel[2])
            poscar.write("\n")
        # write velocities if they exist (again, always direct)
        if velocities:
            poscar.write("Direct\n")
            for site in sites:
                poscar.write(str(site[4][0]) + " " +
                             str(site[4][1]) + " " + str(site[4][2])+"\n")
        if predictors:
            poscar.write("\n")
            for site in sites:
                poscar.write(str(site[5][0]) + " " +
                             str(site[5][1]) + " " + str(site[5][2])+"\n")


class Site(object):

    def __init__(self, specie, position, selective=[True, True, True],
                 velocities=None, predictors=None,
                 direct = True, logger=None):
        """A site, typically a position in POSCAR.

        Parameters
        ----------
        specie : string
            The specie of this site.
        position : ndarray
            The position of the current site as a ndarray of floats.
        selective : ndarray, optional
            The selective tags as a ndarray of booleans. If not
            supplied, defaults to True.
        velocities : ndarray, optional
            The velocities for each position. Defaults to None if not
            supplied.
        predictors : ndarray, optional
            The predictor-corrector coordinates. Defaults to None if not
            supplied.
        direct : bool, optional
            True if the position is in direct coordinates. This is the
            default.
        logger : object, optional
            A standard Python logger object.

        """

        # make sure specie is lowercase
        self.specie = specie.lower()
        self.position = position
        self.selective = selective
        self.velocities = velocities
        self.predictors = predictors
        self.direct = direct

    def get_specie(self):
        return self.specie

    def get_position(self):
        return self.position

    def get_selective(self):
        return self.selective

    def get_velocities(self):
        return self.velocities

    def get_predictors(self):
        return self.predictors
    
    def get_direct(self):
        return self.direct
