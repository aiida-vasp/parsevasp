"""Handle POSCAR."""
# pylint: disable=C0302, consider-using-f-string
import io
import sys
from collections import Counter

import numpy as np

from parsevasp import utils
from parsevasp.base import BaseParser


class Poscar(BaseParser):
    """Class to handle POSCAR."""

    ERROR_UNSUPPORTED_TAG = 300
    ERROR_VASPFOUR = 301
    ERROR_NO_VEL_OR_PRED = 302
    ERROR_SITE_NUMBER = 303
    ERROR_TOO_LARGE_SITE_INDEX = 304
    ERROR_INVALID_ENTRY = 305
    ERROR_NO_DIRECT = 306
    ERROR_NEGATIVE_SCALING = 307
    BaseParser.ERROR_MESSAGES.update({
        ERROR_NEGATIVE_SCALING: 'Currently negative scaling values in POSCAR is not supported.',
        ERROR_VASPFOUR: 'VASP 4 POSCAR is not supported. User, please modernize. ',
        ERROR_NO_VEL_OR_PRED: 'A velocity or predictor-corrector coordinate was not detected.',
        ERROR_SITE_NUMBER: "The supplied 'site_number' is not a number (i.e. the index) "
        'starting from 1 for the site position to be modified.',
        ERROR_TOO_LARGE_SITE_INDEX: 'The supplied site_number is larger than the number of sites.',
        ERROR_INVALID_ENTRY: "Only 'comment', 'unitcell' or 'sites' is allowed as input for entry.",
        ERROR_NO_DIRECT: 'Coordinate should be direct. Did you hack this?'
    })
    ERROR_MESSAGES = BaseParser.ERROR_MESSAGES

    def __init__(
        self,
        poscar_string=None,
        poscar_dict=None,
        file_path=None,
        file_handler=None,
        logger=None,
        prec=None,
        conserve_order=False,
        write_direct=True
    ):
        """Initialize a POSCAR object and set content as a dictionary.

        Parameters
        ----------
        poscar_string : string
            A string containing POSCAR entries. Must contain line
            breaks if multiline, otherwise the POSCAR will be mangled.
        poscar_dict : dict
            A dictionary containing the POSCAR entries.
        file_path : string
            A string containing the file path to the file that is going to be parsed.
        file_handler : object
            A file like object that acts as a handler for the content to be parsed.
        logger : object
            A logger object if you would like to use an external logger for messages
            ejected inside this parser.
        prec : int, optional
            An integer describing how many decimals the users wants
            when printing files.
        conserve_order : bool
            If True, do keep the ordering of the supplied positions
            and atomic species.
        write_direct : bool (optional)
            If True (the default), all write operations are done using direct coordinates,
            while if False in cartesian coordinates.

        """

        super().__init__(file_path=file_path, file_handler=file_handler, logger=logger)

        self._write_direct = write_direct
        self._poscar_dict = poscar_dict
        self._poscar_string = poscar_string
        self._conserve_order = conserve_order

        # Check that only one argument is supplied
        # pylint: disable=R0916
        if (self._poscar_string is not None and self._poscar_dict is not None) or (
            self._poscar_string is not None and self._file_path is not None
        ) or (self._poscar_dict is not None and self._file_path is not None and self._file_handler is not None):
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)
        # Check that at least one is suplpied
        if (
            self._poscar_string is None and self._poscar_dict is None and self._file_path is None and
            self._file_handler is None
        ):
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        # Set precision
        if prec is None:
            self._prec = 12
        else:
            self._prec = prec
        self._width = self._prec + 4

        if self._file_path is not None or self._file_handler is not None:
            # Create dictionary from a file
            self._poscar_dict = self._from_file()

        if self._poscar_string is not None:
            # Create dictionary from a string
            self._poscar_dict = self._from_string()

        if self._poscar_dict is not None:
            # Update site entries to Site objects
            self._from_dict()

        # Store entries
        self.entries = self._poscar_dict

        # Validate dictionary
        self._validate()

    def _from_file(self):
        """Create rudimentary dictionary of entries from a
        file.

        """

        poscar = utils.read_from_file(self._file_path, self._file_handler, encoding='utf8')
        poscar_dict = self._from_list(poscar)
        return poscar_dict

    def _from_string(self):
        """Create rudimentary dictionary of entries from a
        string.

        """

        poscar = self._poscar_string.splitlines(True)
        poscar_dict = self._from_list(poscar)
        return poscar_dict

    def _from_dict(self):
        """If necessary, generate Site objects for each entry in site and/or
        convert to direct coordinates.

        """

        sites = self._poscar_dict['sites']
        for site in sites:
            if not isinstance(site, Site):
                # Entry is not of a Site type, convert it
                direct = site['direct']
                position = site['position']
                velocities = site['velocities']
                if not direct:
                    # Convert to direct
                    position = self._to_direct(position, self._poscar_dict['unitcell'])
                    if velocities is not None:
                        velocities = self._to_direct(velocities, self._poscar_dict['unitcell'])
                    direct = True
                site = Site(site['specie'], position, site['selective'], velocities, site['predictors'], direct)
            else:
                if not site.get_direct():
                    # Cartesian, so convert.
                    position = site.get_position()
                    velocities = site.get_velocities()
                    position = self._to_direct(position, self._poscar_dict['unitcell'])
                    site.set_position(position)
                    if velocities is not None:
                        velocities = self._to_direct(velocities, self._poscar_dict['unitcell'])
                        site.set_velocities(velocities)
                    site.set_direct(True)

    def _from_list(self, poscar):  # pylint: disable=R0915
        """
        Go through the list and analyze for = and ; in order to
        disentangle grouped entries etc.

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

        comment = poscar[0].replace('#', '').strip()
        vasp5 = True
        # Check for VASP 5 POSCAR
        if utils.is_numbers(poscar[5]):
            vasp5 = False
        # Set direct, test is done later
        direct = True
        # Set selective, test is done later
        selective = False
        # Check scaling factor
        scaling = float(poscar[1].split()[0])
        if scaling < 0.0:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NEGATIVE_SCALING])
            sys.exit(self.ERROR_NEGATIVE_SCALING)
        unitcell = [[0.0 for y in range(3)] for x in range(3)]
        nions = 0
        spec = None
        loopmax = 8
        if vasp5:
            # Could go straight to numpy with fromstring, consider
            # to change in the future
            unitcell[0] = [float(x) for x in poscar[2].split()]
            unitcell[1] = [float(x) for x in poscar[3].split()]
            unitcell[2] = [float(x) for x in poscar[4].split()]
            unitcell = np.asarray(unitcell)
            # Apply scaling factor
            unitcell = scaling * unitcell
            spec = poscar[5].split()
            atoms = [int(x) for x in poscar[6].split()]
            for num_ions in atoms:
                nions = nions + num_ions
            if poscar[7][0].lower() == 's':
                selective = True
                loopmax = 9
                if not poscar[8][0].lower() == 'd':
                    direct = False
            if not selective:
                if not poscar[7][0].lower() == 'd':
                    direct = False
        else:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_VASPFOUR])
            sys.exit(self.ERROR_VASPFOUR)

        # Create site objects
        specie_slot = 0
        index = 1
        sites_temp = []
        velocities = False
        predictor = False
        # Loop positions
        for i in range(nions):
            # Fetch specie
            if index > atoms[specie_slot]:
                specie_slot = specie_slot + 1
                index = 1
            specie = spec[specie_slot]
            # Fetch positions
            line = poscar[i + loopmax].split()
            position = np.zeros(3)
            if direct:
                # When we have direct coordinates, the
                # positions should not be scaled by the
                # scaling factor
                scaling = 1.0
            position[0] = float(line[0]) * scaling
            position[1] = float(line[1]) * scaling
            position[2] = float(line[2]) * scaling
            if not direct:
                # Convert to direct
                position = self._to_direct(position, unitcell)
            # Fetch selective flags
            flags = [True, True, True]
            if selective:
                if 'f' in line[3].lower():
                    flags[0] = False
                if 'f' in line[4].lower():
                    flags[1] = False
                if 'f' in line[5].lower():
                    flags[2] = False
            # Create a site object and add to sites list
            index = index + 1
            velo = None
            pred = None
            sites_temp.append([specie, position, flags, velo, pred])
        # Now check if there is more in the POSCAR
        loopmax_pos = nions + loopmax
        if len(poscar) > loopmax_pos:
            first_char = poscar[loopmax_pos][0].lower()
            if first_char in ('d', 'c'):
                velocities = True
                if first_char == 'c':
                    # Make sure we convert velocities to direct
                    direct = False
            elif poscar[loopmax_pos].replace(' ', '') == '\n':
                predictor = True
        # Now check that the next line is in fact a coordinate
        loopmax_pos = loopmax_pos + 1
        # Allow for blank lines at the end of the positions
        if len(poscar) > loopmax_pos:
            if not utils.is_number(poscar[loopmax_pos].split()[0]):
                self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_VEL_OR_PRED])
                sys.exit(self.ERROR_NO_VEL_OR_PRED)
        else:
            # But make sure the predictor is set back to False
            # if we only have a blank line and nothing else following
            # the coordinates
            predictor = False
        if velocities:
            for i in range(nions):
                # Fetch velocities
                line = poscar[i + loopmax_pos].split()
                vel = np.zeros(3)
                vel[0] = float(line[0])
                vel[1] = float(line[1])
                vel[2] = float(line[2])
                if not direct:
                    # Convert to direct
                    vel = self._to_direct(vel, unitcell)
                sites_temp[i][3] = vel
            # Now check if there is predictor-corrector coordinates following
            # the velocities
            loopmax_pos = nions + loopmax_pos
            if len(poscar) > loopmax_pos:
                if poscar[loopmax_pos].replace(' ', '') == '\n':
                    loopmax_pos = loopmax_pos + 1
                    if utils.is_number(poscar[loopmax_pos].split()[0]):
                        for i in range(nions):
                            line = poscar[i + loopmax_pos].split()
                            pre = np.zeros(3)
                            pre[0] = float(line[0])
                            pre[1] = float(line[1])
                            pre[2] = float(line[2])
                            sites_temp[i][4] = pre
        # Do one final loop to create the objects and read
        # predictors if they exist
        sites = []
        loopmax_pos = nions + loopmax + 1
        for i in range(nions):
            pre = np.zeros(3)
            if predictor and not velocities:
                line = poscar[i + loopmax_pos].split()
                pre[0] = float(line[0])
                pre[1] = float(line[1])
                pre[2] = float(line[2])
                sites_temp[i][4] = pre
            site = Site(
                sites_temp[i][0],
                sites_temp[i][1],
                selective=sites_temp[i][2],
                velocities=sites_temp[i][3],
                predictors=sites_temp[i][4]
            )
            sites.append(site)

        # Build dictionary and convert to NumPy
        poscar_dict = {}
        poscar_dict['comment'] = comment
        poscar_dict['unitcell'] = np.asarray(unitcell)
        poscar_dict['sites'] = sites
        return poscar_dict

    def modify(self, entry, value, site_number=None):
        """
        Modify an entry tag in the Poscar dictionary.
        If it is not found add it.

        Parameters
        ----------
        entry : string
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

        # Check allowed entries
        self._check_allowed_entries(entry)
        # Check that entries exists
        self._check_dict()
        if site_number is not None:
            # Check that a Site() object is supplied
            self._check_site(value)
            # Check site number
            self._check_site_number(site_number)
            # Check that position is an integer
            if not utils.is_number(site_number):
                self._logger.error(self.ERROR_MESSAGES[self.ERROR_SITE_NUMBER])
                sys.exit(self.ERROR_SITE_NUMBER)
            self.entries['sites'][site_number] = value
        else:
            if entry == 'sites':
                self._check_sites(sites=value)
            if entry == 'comment':
                self._check_comment(comment=value)
            if entry == 'unitcell':
                self._check_unitcell(unitcell=value)

            self.entries[entry] = value

    def delete_site(self, site_number):
        """
        Delete the site with the supplied
        number.

        Parameters
        ----------
        site_number : int
            The site number to be deleted, starting
            from 0.

        """

        self._check_sites()
        self._check_site_number(site_number)
        del self.entries['sites'][site_number]

    def add_site(self, site_number):
        """
        Add a site with the supplied
        number. If not supplied, add at the end
        of the last element of the specie group

        Parameters
        ----------
        site_number : int
            The site number to be deleted, starting
            from 0.

        """

        raise NotImplementedError

    def _check_dict(self):
        """Check that entries is present.

        """

        try:
            _ = self.entries
        except AttributeError:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_ENTRIES])
            sys.exit(self.ERROR_NO_ENTRIES)

    def _check_allowed_entries(self, entry):
        """
        Check the allowed values of entry.

        Parameters
        ----------
        entry : string
            Contains the entry to be checked.

        """

        if not (('comment' in entry) or ('unitcell' in entry) or ('sites' in entry)):
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_INVALID_ENTRY])
            sys.exit(self.ERROR_INVALID_ENTRY)

    def _check_unitcell(self, unitcell=None):
        """
        Check that the unitcell entries are present and
        are of a 3x3 ndarray type.

        Parameters
        ----------
        unitcell, optional
            The unitcell to be checked. If not supplied the
            'unitcell' key in the 'entries' is checked.

        """

        if unitcell is None:
            try:
                unitcell = self.entries['unitcell']
            except KeyError:
                self._logger.error(self.ERROR_MESSAGES[f"{self.ERROR_NO_KEY} The key in question is 'unitcell'."])
                sys.exit(self.ERROR_NO_KEY)

        if (not isinstance(unitcell, np.ndarray)) \
           or (unitcell.shape != (3, 3)):
            self._logger.error(
                f"{self.ERROR_MESSAGES[self.ERROR_KEY_INVALID_TYPE]} The value of 'unitcell' is not an 3x3 ndarray."
            )
            sys.exit(self.ERROR_KEY_INVALID_TYPE)

    def _check_comment(self, comment=None):
        """
        Check that the comment entry is present and
        is a string.

        Parameters
        ----------
        comment, optional
            The comment to be checked. If not supplied the
            'comment' key in the 'entries' is checked.

        """
        if comment is None:
            try:
                comment = self.entries['comment']
            except KeyError:
                self._logger.error(self.ERROR_MESSAGES[f"{self.ERROR_NO_KEY} The key in question is 'comment'."])
                sys.exit(self.ERROR_NO_KEY)
        # Allow None for comment
        if self.entries['comment'] is not None:
            if not isinstance(comment, str):
                self._logger.error(
                    f"{self.ERROR_MESSAGES[self.ERROR_KEY_INVALID_TYPE]} The key 'comment' is not a string."
                )
                sys.exit(self.ERROR_KEY_INVALID_TYPE)

    def _check_sites(self, sites=None):
        """
        Check that the sites entries are present.

        Parameters
        ----------
        sites, optional
            The sites to be checked. If not supplied the
            'sites' key in the 'entries' is checked.

        """

        if sites is None:
            try:
                sites = self.entries['sites']
            except KeyError:
                self._logger.error(self.ERROR_MESSAGES[f"{self.ERROR_NO_KEY} The key in question is 'sites'."])
                sys.exit(self.ERROR_NO_KEY)
        if not isinstance(sites, list):
            self._logger.error(f"{self.ERROR_MESSAGES[self.ERROR_KEY_INVALID_TYPE]} The key 'sites' is not a list.")
            sys.exit(self.ERROR_KEY_INVALID_TYPE)

    def _check_site(self, site=None):
        """
        Check that the site entry is a Site() object.

        Parameters
        ----------
        site, optional
            The site to be checked. If not supplied the entries under
            the 'sites' key in the 'entries' is checked.

        """
        if site is None:
            try:
                sites = self.entries['sites']
            except KeyError:
                self._logger.error(self.ERROR_MESSAGES[f"{self.ERROR_NO_KEY} The key in question is 'sites'."])
                sys.exit(self.ERROR_NO_KEY)
            for _site in sites:
                if not isinstance(_site, Site):
                    self._logger.error(
                        f'{self.ERROR_MESSAGES[self.ERROR_KEY_INVALID_TYPE]} '
                        "The elements of the key 'sites' are not Site() objects."
                    )
                    sys.exit(self.ERROR_KEY_INVALID_TYPE)
        else:
            if not isinstance(site, Site):
                self._logger.error(
                    f"{self.ERROR_MESSAGES[self.ERROR_KEY_INVALID_TYPE]} The key 'site' is not a Site() object."
                )
                sys.exit(self.ERROR_KEY_INVALID_TYPE)

    def _check_site_number(self, site_number):
        """
        Check that the site_number is an int and that
        it is not out of bounds.

        Parameters
        ----------
        site_number : int
            The site_number to be checked

        """

        if not isinstance(site_number, int):
            self._logger.error(
                f"{self.ERROR_MESSAGES[self.ERROR_KEY_INVALID_TYPE]} The key 'site_number' is not an integer."
            )
            sys.exit(self.ERROR_KEY_INVALID_TYPE)
        sites = self.entries['sites']
        if site_number > (len(sites) - 1):
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_TOO_LARGE_SITE_INDEX])
            sys.exit(self.ERROR_TOO_LARGE_SITE_INDEX)

    def _validate(self):
        """Validate the content of entries

        """

        self._check_dict()
        self._check_comment()
        self._check_unitcell()
        self._check_sites()

    def _sort_and_group_sites(self):
        """
        Sort and group the positions and species to
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
            Contains the occupancy of each specie in the same order as
            'species'.
        selective : bool
            True if any selective flags are enabled, False otherwise.

        """

        sites = []
        species = []
        selective = False
        velocities = False
        predictors = False
        for site in self.entries['sites']:
            specie = site.get_specie()
            select = site.get_selective()
            position = site.get_position()
            direct = site.get_direct()
            vel = site.get_velocities()
            pre = site.get_predictors()
            if direct is False:
                # Make sure it is direct as the writer only
                # supports this
                self._logger.error(self.ERROR_MESSAGES[self.ERROR_NO_DIRECT])
                sys.exit(self.ERROR_NO_DIRECT)
            if False in select:
                selective = True
            if vel is not None:
                velocities = True
            if pre is not None:
                predictors = True
            sites.append([specie, position, select, direct, vel, pre])
            species.append(specie)

        if not self._conserve_order:
            # Find unique entries and their number
            counter = Counter(species)
            # Counter does not order, so order now with the
            # least occuring element first (typical for compounds)
            sorted_keys = sorted(counter, key=counter.get)
            species = []
            num_species = []
            for key in sorted_keys:
                species.append(key)
                num_species.append(counter[key])

            # Now make sure the sites is on the same order
            ordered_sites = []
            for specie in species:
                ordered_sites.extend([site for site in sites if specie == site[0]])
            # Consider to also sort on coordinate after specie

            return ordered_sites, species, num_species, selective, velocities, predictors

        # Do not order, but we still need to group similar species
        # that follow each other
        num_species = []
        species_concat = []
        for specie in species:
            if species_concat and species_concat[-1] == specie:
                num_species[-1] = num_species[-1] + 1
            else:
                species_concat.append(specie)
                num_species.append(1)
        return sites, species_concat, num_species, selective, velocities, predictors

    def _get_key(self, item):
        """Key fetcher for the sorted function."""

        return item[0]

    def _to_direct(self, position_cart, unitcell):
        """
        Transforms the position from cartesian to direct
        coordinates.

        Parameters
        ----------
        position_cart : ndarray
            | Dimension: (3)

            An ndarray containing the position in cartesian coordinates.
        unitcell : ndarray
            | Dimension: (3,3)

            Contains the unitcell.

        Returns
        -------
        position : ndarray
            An ndarray containing the position in direct coordinates.

        """

        position = utils.cart_to_dir(position_cart, unitcell)

        return position

    def _to_cart(self, position_dir, unitcell):
        """
        Transforms the position from direct to cartesian
        coordinates.

        Parameters
        ----------
        position : ndarray
            | Dimension: (3)

            An ndarray containing the position in direct coordinates.
        unitcell : ndarray
            | Dimension: (3,3)

            Contains the unitcell.

        Returns
        -------
        position : ndarray
            An ndarray containing the position in cartesian coordinates.

        """

        position = utils.dir_to_cart(position_dir, unitcell)

        return position

    def get(self, tag):
        """
        Return the value and comment of the entry with tag.

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

    def get_dict(self, direct=True):
        """
        Get a true dictionary containing the entries in an
        POSCAR compatible fashion.

        Returns
        -------
        dictionary : dict
            A dictionary on POSCAR compatible form.
        direct : bool, optional
            If True, all coordinates are returned in direct, otherwise
            in cartesian.

        """

        dictionary = {}
        for key, entry in self.entries.items():
            if key == 'sites':
                sites_temp = []
                for element in entry:
                    position = element.get_position()
                    velocities = element.get_velocities()
                    temp_direct = element.get_direct()
                    if not direct:
                        # Convert to cartesian
                        position = self._to_cart(position, self.entries['unitcell'])
                        if velocities is not None:
                            velocities = self._to_cart(velocities, self.entries['unitcell'])
                        temp_direct = False
                    sites_temp.append({
                        'specie': element.get_specie().capitalize(),
                        'position': position,
                        'selective': element.get_selective(),
                        'velocities': velocities,
                        'predictors': element.get_predictors(),
                        'direct': temp_direct
                    })
                dictionary[key] = sites_temp

            else:
                dictionary[key] = entry

        return dictionary

    def get_string(self):
        """
        Get a string containing the entries in a POSCAR
        compatible fashion. Each line is broken by a newline
        character

        Returns
        -------
        poscar_string : str
            A string containing the POSCAR entries of the
            current instance.

        """

        string_object = io.StringIO()
        self._write(file_handler=string_object)
        poscar_string = string_object.getvalue()
        string_object.close()

        return poscar_string

    def _write(self, file_handler, **kwargs):  # pylint: disable=R0915
        """
        Write POSCAR like files to a file or string

        Parameters
        ----------
        file_handler : object
            Either a file object or a StringIO object.

        """

        self._validate()
        entries = self.entries
        comment = entries['comment']
        unitcell = entries['unitcell']
        # Sort and group to VASP specifications
        sites, species, num_species, selective, velocities, predictors = \
            self._sort_and_group_sites()
        # Update comment
        compound = ''
        for index, specie in enumerate(species):
            if num_species[index] == 1:
                num_string = ''
            else:
                num_string = str(num_species[index])
            compound = compound + str(specie).capitalize() + num_string
        compound = 'Compound: ' + compound + '.'
        if comment is None:
            comment = '# ' + compound
        elif compound not in comment:
            comment = '# ' + compound + ' Old comment: ' + comment
        else:
            comment = '# ' + comment
        file_handler.write(comment + '\n')
        # We avoid usage of the scaling factor
        file_handler.write('{:{width}.{prec}f}\n'.format(1.0, prec=self._prec, width=self._width))
        # Write unitcell
        for i in range(3):
            file_handler.write(
                '{:{width}.{prec}f} {:{width}.{prec}f} '
                '{:{width}.{prec}f}\n'.format(
                    unitcell[i][0], unitcell[i][1], unitcell[i][2], prec=self._prec, width=self._width
                )
            )
        # Write specie types
        tempostring = ''
        for specie in species:
            tempostring = tempostring + f'{specie.capitalize():5s} '
        file_handler.write(f'{tempostring.rstrip()}\n')
        # Write number of species
        tempostring = ''
        for number in num_species:
            tempostring = tempostring + f'{number:5d} '
        file_handler.write(f'{tempostring.rstrip()}\n')
        # Write selective if any flags are True
        if selective:
            file_handler.write('Selective dynamics\n')
        if not self._write_direct:
            file_handler.write('Cartesian\n')
        else:
            file_handler.write('Direct\n')

        # Write positions
        for site in sites:
            if self._write_direct:
                _site = site[1]
            else:
                _site = self._to_cart(site[1][0:3], unitcell)
            file_handler.write(
                '{:{width}.{prec}f} {:{width}.{prec}f} '
                '{:{width}.{prec}f}'.format(_site[0], _site[1], _site[2], prec=self._prec, width=self._width)
            )
            if selective:
                sel = ['T', 'T', 'T']
                flags = site[2]
                for index, flag in enumerate(flags):
                    if not flag:
                        sel[index] = 'F'
                file_handler.write(f' {sel[0]} {sel[1]} {sel[2]}')
            file_handler.write('\n')

        # Write velocities if they exist
        if velocities:
            if self._write_direct:
                file_handler.write('Direct\n')
            else:
                file_handler.write('Cartesian\n')
            for site in sites:
                if self._write_direct:
                    _site = site[4]
                else:
                    _site = self._to_cart(site[4][0:3], unitcell)
                file_handler.write(
                    '{:{width}.{prec}f} {:{width}.{prec}f} '
                    '{:{width}.{prec}f}\n'.format(_site[0], _site[1], _site[2], prec=self._prec, width=self._width)
                )

        # Write predictors if they exist
        if predictors:
            file_handler.write('\n')
            for site in sites:
                file_handler.write(
                    '{:{width}.{prec}f} {:{width}.{prec}f} '
                    '{:{width}.{prec}f}\n'.format(
                        site[5][0], site[5][1], site[5][2], prec=self._prec, width=self._width
                    )
                )


class Site:
    """Class to represent atomic site."""

    def __init__(self, specie, position, selective=None, velocities=None, predictors=None, direct=True):
        """
        A site, typically a position in POSCAR.

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

        """

        # make sure specie is lowercase
        self.specie = specie.lower()
        self.position = position
        if selective is None:
            self.selective = [True, True, True]
        else:
            self.selective = selective
        self.velocities = velocities
        self.predictors = predictors
        self.direct = direct

    def get_specie(self):
        """
        Return the specie.

        Returns
        -------
        specie : string
            A string containing the atomic capitalized atomic specie.

        """

        specie = self.specie.capitalize()
        return specie

    def set_position(self, position):
        """
        Sets the position.

        Parameters
        ----------
        position : ndarray
            The coordinate that are to be set as position.

        """

        self.position = position

    def get_position(self):
        """
        Return the position.

        Returns
        -------
        position : ndarray

        """

        position = self.position
        return position

    def get_selective(self):
        """
        Return the selective flags.

        Returns
        -------
        selective : list
            A list of three bool, either True or False, depending on
            which directions to perform selective dynamics.

        """

        selective = self.selective
        return selective

    def set_velocities(self, velocities):
        """
        Sets the velocities.

        Parameters
        ----------
        velocities : ndarray
            The three velocities that are to be set as velocities.

        """

        self.velocities = velocities

    def get_velocities(self):
        """
        Return the velocities.

        Returns
        -------
        velocities : ndarray
            An ndarray of three floats containing the velocities along
            each direction.

        """

        velocities = self.velocities
        return velocities

    def get_predictors(self):
        """
        Return the predictors.

        Returns
        -------
        predictors : ndarray
            An ndarray of three floats containing the predictors along
            each direction.

        """

        predictors = self.predictors
        return predictors

    def set_direct(self, direct):
        """
        Sets direct.

        Parameters
        ----------
        direct : bool
            The boolean that determines if the position and velocities are
            in cartesian or direct coordinates.

        """

        self.direct = direct

    def get_direct(self):
        """
        Return the direct status of the coordinate.

        Returns
        -------
        direct : bool
            True if the coordinates are given in direct coordinates,
            otherwise for direct, False.

        """

        direct = self.direct
        return direct
