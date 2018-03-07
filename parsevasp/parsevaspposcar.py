#!/usr/bin/python
import sys
import logging
import numpy as np

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

        if logger is None:
            logger = logging.getLogger(__name__)
        self.logger = logger

        # check that only one argument is supplied
        if (poscar_string is not None and poscar_dict is not None) \
           or (poscar_string is not None and file_path is not None) \
           or (poscar_dict is not None and file_path is not None):
            self.logger.error("Please only supply one argument when "
                              "initializing Poscar. Exiting.")
            sys.exit(1)
        # check that at least one is suplpied
        if (poscar_string is None and poscar_dict is None
                and file_path is None):
            self.logger.error("Please supply one argument when "
                              "initializing Poscar. Exiting.")
            sys.exit(1)

        if file_path is not None:
            # create dictionary from a file
            poscar_dict = self._from_file()

        if poscar_string is not None:
            # create dictionary from a string
            poscar_dict = self._from_string()

        # validate dictionary
        self.validate()

        # store entries
        self.entries = poscar_dict

    def _from_file(self):
        """Create rudimentary dictionary of entries from a
        file.

        """

        poscar = utils.read_file(self.file_path)
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
            self.logger.error("Currently negative scaling values in POSCAR is "
                              "not supported. Exiting.")
            sys.exit(1)
        lattice = [[0.0 for y in range(3)] for x in range(3)]
        nions = 0
        spec = None
        if (vasp5):
            lattice[0] = [float(x) for x in poscar[2].split()]
            lattice[1] = [float(x) for x in poscar[3].split()]
            lattice[2] = [float(x) for x in poscar[4].split()]
            spec = poscar[5].split()
            atoms = [int(x) for x in poscar[6].split()]
            for i in range(len(atoms)):
                nions = nions + atoms[i]
            if (poscar[7][0] == "s") or (poscar[7][0] == "S"):
                self.logger.info("Selective tag found")
                selective = True
                loopmax = 9
                if not (poscar[8][0] == "D" or poscar[8][0] == "d"):
                    self.logger.error(
                        "Please supply a POSCAR in direct coordinates. "
                        "Exiting.")
                    sys.exit(1)
            else:
                loopmax = 8

            if not (poscar[7][0] == "D" or poscar[7][0] == "d") and selective:
                self.logger.error(
                    "Please supply a POSCAR in direct coordinates. Exiting.")
                sys.exit(1)
        else:
            self.logger.error(
                "VASP 4 POSCAR is not supported. User, please modernize. "
                "Exiting.")
            sys.exit(1)

        # create site objects
        position = np.zeros(3)
        specie = ""
        specie_slot = 0
        index = 1
        sites = []
        for i in range(nions):
            # fetch specie
            if index > atoms[specie_slot]:
                specie_slot = specie_slot + 1
                index = 1
            specie = spec[specie_slot]
            # fetch positions
            line = poscar[i + loopmax].split()
            position[0] = float(line[0])
            position[1] = float(line[1])
            position[2] = float(line[2])
            # fetch selective flags
            selective = [False, False, False]
            if (selective):
                if "t" in line[3].lower():
                    selective[0] = True
                if "t" in line[4].lower():
                    selective[1] = True
                if "t" in line[5].lower():
                    selective[2] = True
            # create a site object and add to sites list
            site = Site(specie, position, selective)
            sites.append(site)
            index = index + 1
        # apply scaling factor
        lattice = [[x * scaling for x in y] for y in lattice]
        # build dictionary and convert to NumPy
        poscar_dict = {}
        poscar_dict["comment"] = comment.split('#')[1]
        poscar_dict["lattice"] = np.asarray(lattice)
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
            an ndarray for the lattice or a Site object
            for a position.
        site_number : int, optional
            The site to be modified. If not supplied
            the value have to be a list of Site objects.

        """

        if ("comment" not in entry) \
           or ("lattice" not in entry) \
           or ("sites" not in entry):
            self.logger.error("Only 'comment', 'lattice' or "
                              "'sites' is allowed as input for "
                              "entry. Exiting.")
            sys.exit(1)

        if site_number is not None:
            # check that sites is supplied if a site is supplied
            if "sites" not in entry:
                self.logger.error("If 'site_number' is "
                                  "supplied, please also make "
                                  "sure that 'entry' is set to "
                                  "'sites'. Exiting.")
                sys.exit(1)
            # check that position is an integer
            if not utils.is_number(site_number):
                self.logger.error("The supplied 'site_number' is "
                                  "not a number (i.e. the index) "
                                  "starting from 1 for the site "
                                  "position to be modified. Exiting.")
                sys.exit(1)
            # check that value is of the right instance
            if not isinstance(value, Site):
                self.logger.error("The supplied 'value' is not an "
                                  "instance of Site(). Exiting.")
                sys.exit(1)
            self.poscar_dict["sites"][site_number] = value
        else:
            if "sites" in entry:
                # the value have to be a list
                if not isinstance(value, list):
                    self.logger.error("The supplied sites 'value' "
                                      "have to be a list of Site() "
                                      "objects. Exiting.")
                    sys.exit(1)
                else:
                    for element in value:
                        if not isinstance(element, Site):
                            self.logger.error("The supplied sites 'value' "
                                              "have to be a list of Site() "
                                              "objects. Exiting.")
                            sys.exit(1)

            if entry == "comment":
                # check if the supplied value is a string
                if not isinstance(value, str):
                    self.logger.error("The supplied comment 'value' is "
                                      "not a string. Exiting.")
                    sys.exit(1)

            if entry == "lattice":
                # check that the supplied type is ndaray
                if (not isinstance(value, np.ndarray)) \
                   or (value.shape != (3, 3)):
                    self.logger.error("The supplied lattice 'value' is "
                                      "not an 3x3 ndarray. Exiting.")
                    sys.exit(1)

            self.poscar_dict[entry] = value

    def delete_site(self, site_number):
        """Delete the site with the supplied
        number.

        Parameters
        ----------
        site_number : int
            The site number to be deleted, starting
            from 1.

        """

        if not isinstance(site_number, int):
            self.logger.error("The supplied 'site_number' is "
                              "not an integer. Exiting.")
            sys.exit(1)
        # add test for out of bounds
        del self.poscar_dict["sites"][site_number]

    def add_site(self, site_number):
        """Add a site with the supplied
        number. If not supplied, add at the end
        of the last element of the specie group

        Parameters
        ----------
        site_number : int
            The site number to be deleted, starting
            from 1.

        """

        if not isinstance(site_number, int):
            self.logger.error("The supplied 'site_number' is "
                              "not an integer. Exiting.")
            sys.exit(1)
        # add test for out of bounds
        del self.poscar_dict["sites"][site_number]

    def get(self, tag, comment=False):
        """Return the value and comment of the entry with tag.

        Parameters
        ----------
        tag : string
            The entry tag of the POSCAR entry.        
        comment : bool, optional
            If set to True, the comment is also returned, otherwise
            not.

        Returns
        -------
        value : string, int, float or list
            The value of the tag entry
        com : string, optional
            If comment is set to True, the comment associated with
            the tag entry will be returned here.

        """

        value = None
        com = None
        try:
            value = self.entries[tag].get_value()
            if comment:
                com = self.entries[tag].get_value()
        except KeyError:
            pass

        if comment:
            return value, com
        else:
            return value

    def write(self, file_path, comments=False):
        """Write the content of the current Poscar instance to
        file.

        Parameters
        ----------
        file_path : string
            The location to write POSCAR.
        comments : bool, optional
            If set to true, the comments are also dumped to the file,
            else not.

        """

        poscar = utils.file_handler(file_path, status='w')
        # write in alfabetical order
        keys = sorted(self.entries)
        entries = self.entries
        for key in keys:
            entry = entries[key]
            value = entry.get_value()
            comment = entry.get_comment()
            if comment is None or not comments:
                comment = ""
            else:
                comment = " # " + comment
            value = self._convert_value_to_string(value)
            string = str(key.upper()) + " = " + value + comment + "\n"
            poscar.write(string)
        utils.file_handler(file_handler=poscar)


class Site(object):

    def __init__(self, specie, position, selective=None, logger=None):
        """A site, typically a position in POSCAR.

        Parameters
        ----------
        specie : string
            The specie of this site.
        position : ndarray
            The position of the current site as a ndarray of floats.
        selective : ndarray, optional
            The selective tags as a ndarray of booleans. If not
            supplied, defaults to None.
        logger : object, optional
            A standard Python logger object.

        """

        if logger is None:
            logger = logging.getLogger(__name__)
        self.logger = logger

        # make sure specie is lowercase
        self.specie = specie.lower()
        self.position = position
        self.selective = selective

    def get_specie(self):
        return self.specie

    def get_position(self):
        return self.position

    def get_selective(self):
        return self.selective

    def _clean_entry(self, tag, value, comment):
        """Cleans the tag and value, for instance makes sures
        there are no spaces around the tag, that it is in lower case,
        and that the value is build into a list of several different
        parameters are given.

        Parameters
        ----------
        tag : string
            The entry tag of the POSCAR entry.
        value : string
            The entry value of the POSCAR entry.
        comment : string
            The entry comment of the POSCAR entry.

        Returns
        ------
        clean_tag : string
            The entry tag of the POSCAR entry, cleaned and checked.
        clean_value : string
            The entry value of the POSCAR entry, cleaned and checked.
        clean_comment : string
            The entry comment of the POSCAR entry, cleaned and checked.

        """

        # remove possible spaces on tag
        clean_tag = tag.split()

        # make sure tag is lowerscore
        clean_tag = clean_tag[0].lower()

        # possible solutions for the value are (if POSCAR is read):
        # 1 - integer
        # .TRUE. - bool
        # Something - a string
        # Something anotherthing - set of strings
        # 1 2 3 4 - set of integers
        # 1.0 - float
        # 1.0 2.0 3.0 - set of floats

        # however, let us also open for the fact that users might
        # give a value what they would in POSCAR

        # make sure we keep compatibility between Python 2 and 3
        try:
            basestring
        except NameError:
            basestring = str

        if isinstance(value, basestring):
            # values is a string, so we have read an POSCAR or the user
            # try to assign an entry with just a string
            values = value.split()

            # if we have some kind of set, check if all are ints, floats
            # or strings
            content_type = []
            clean_value = []
            for element in values:
                cnt_type = self._test_string_content(element)
                content_type.append(cnt_type)
                if cnt_type == 'int':
                    cnt = int(element)
                elif cnt_type == 'float':
                    cnt = float(element)
                else:
                    # we have here also have a bool, so check that
                    cnt = self._test_string_for_bool(element)
                clean_value.append(cnt)

            # check if all values are the same type (they should be)
            if not all(x == content_type[0] for x in content_type):
                self.logger.error("All values of the tag " + clean_tag.upper() +
                                  " are not of the same type. Did you "
                                  "possibly forget a comment tag(#)? Exiting.")
                sys.exit(1)

        else:
            # if the user wants to assign an element as a list, int,
            # float or bool, accept this
            if not (isinstance(value, int) or isinstance(value, float)
                    or isinstance(value, bool) or (type(value) is list)):
                self.logger.error("The type: " + str(type(value)) + " of the "
                                  "supplied value for the POSCAR tag is not "
                                  "recognized. Exiting.")
                sys.exit(1)
            clean_value = value

        # finally clean comment by removing any redundant spaces
        if comment is not None:
            clean_comment = comment.split()[0]
        else:
            clean_comment = None

        return clean_tag, clean_value, clean_comment

    def _test_string_for_bool(self, string):
        """Detects if string contains Fortran bool.

        Parameters
        ----------
        string : string
            A string containing what is to be tested.

        Returns
        -------
        string : bool or string
            The return is the detected boolean or the input string,
            depending on what is detected.

        """
        if ".True." in string \
           or ".TRUE." in string \
           or ".true." in string:
            return True
        elif ".False." in string \
             or ".FALSE." in string \
             or ".false." in string:
            return False
        else:
            return string

    def _test_string_content(self, string):
        """Detects if string is integer, float or string.

        Parameters
        ----------
        string : string
            An input string to be tested.

        Returns
        -------
        string
            A string with value 'int' if input is an integer,
            'float' if the input is a float and 'string' if it
            is just a regular string.

        """
        try:
            float(string)
            return 'int' if string.count('.') == 0 else 'float'
        except ValueError:
            return 'string'
