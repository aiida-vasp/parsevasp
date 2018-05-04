#!/usr/bin/python
import sys
import os
import numpy as np
import logging
import mmap
import StringIO
# sys.path.append(os.path.dirname(__file__))
# print sys.path

import utils
import constants


class Incar(object):

    def __init__(self, incar_string=None, incar_dict=None,
                 file_path=None, logger=None, prec = None):
        """Initialize an INCAR object and set content as a dictionary.

        Parameters
        ----------
        incar_string : string, optional
            A string containing INCAR entries. Must contain line
            breaks if multiline, otherwise the INCAR will be mangled.
        incar_dict : dict, optional
            A dictionary containing the INCAR entries.
        file_path : string, optional
            The file path in which the INCAR is read.
        logger : object, optional
            A standard Python logger object.
        prec : int, optional
            An integer describing how many decimals the users wants
            when printing files.

        """
        print(incar_string, incar_dict, file_path)
        self.file_path = file_path
        self.incar_dict = incar_dict
        self.incar_string = incar_string

        # set logger
        if logger is not None:
            self._logger = logger
        else:
            logging.basicConfig(level=logging.DEBUG)
            self._logger = logging.getLogger('IncarParser')

        # set precision
        if prec is None:
            self._prec = 12
        else:
            self._prec = prec
        self._width = self._prec + 4    

        # check that only one argument is supplied
        if (incar_string is not None and incar_dict is not None) \
           or (incar_string is not None and file_path is not None) \
           or (incar_dict is not None and file_path is not None):
            self._logger.error("Please only supply one argument when "
                              "initializing Incar. Exiting.")
            sys.exit(1)
            
        # check that at least one is supplied
        if (incar_string is None and incar_dict is None
                and file_path is None):
            self._logger.error("Please supply one argument when "
                              "initializing Incar. Exiting.")
            sys.exit(1)

        if file_path is not None:
            # create list from a file
            incar_list = self._from_file()

        if incar_string is not None:
            # create list from a string
            incar_list = self._from_string()

        if incar_dict is None:
            incar = self._from_list(incar_list)
        else:
            incar = self._from_dict(incar_dict)
            
        # validate dictionary
        self.validate()

        # store entries
        self.entries = incar

    def _from_file(self):
        """Create rudimentary dictionary of entries from a
        file.

        """

        incar = utils.readlines_from_file(self.file_path)
        return incar

    def _from_string(self):
        """Create rudimentary dictionary of entries from a
        string.

        """

        incar = self.incar_string.splitlines()
        return incar

    def _from_list(self, incar):
        """Go through the list and analyze for = and ; in order to
        deentangle grouped entries etc. Also set up IncarItem elements.

        Parameters
        ----------
        incar : list
            A list of strings containing each line in the INCAR file.

        Returns
        -------
        incar_dict : dictionary
            A dictionary containing each INCAR tag as a key with the
            associated element.

        Notes
        -----
        No checking for consistency is done here. We do this at a later step
        in order to be able to keep the input methods as clean as posible.

        """

        incar_dict = {}
        for line in incar:
            # check for comment at the start of a line, if so skip
            comment = not line.split('#')[0].split()
            if not comment:
                # now check if a line contains a comment, if so, only take
                # what is in front of this, truly...
                splitted = line.split('#', 1)
                if len(splitted) > 1:
                    comment = splitted[1]
                else:
                    comment = None
                # now split on ; as we could have combined entries
                splitted = splitted[0].split(';')
                for ntry in splitted:
                    # then split on = and analyze each entry
                    final_split = ntry.split("=")
                    if len(final_split) > 2:
                        self._logger.error("Detected two equal signs for an entry in "
                                          "the INCAR file. The following line contains "
                                          "the problem:\n\n" +
                                          ntry + "\n\nPlease correct. Exiting.")
                        sys.exit(1)
                    if len(final_split) == 1:
                        self._logger.error("Detected a comment line that does not start "
                                          "with a #. Please correct and be consistent. "
                                          "Exiting.")
                        sys.exit(1)
                    tag = final_split[0]
                    value = final_split[1]
                    # create new instance of entry
                    entry = IncarItem(tag, value, comment)
                    clean_tag = entry.get_tag()
                    if clean_tag in incar_dict:
                        self._logger.info("Tag " + entry.get_tag() + " already "
                                         "found in the INCAR dictionary, "
                                         "overwriting it.")
                    incar_dict[clean_tag] = entry

        return incar_dict

    def _from_dict(self, incar):
        """Go through the dict and setup the IncarItem elements.

        Parameters
        ----------
        incar : dict
            A dict containing each line in the INCAR file.

        Returns
        -------
        incar_dict : dictionary
            A dictionary containing each INCAR tag as a key with the
            associated element.

        Notes
        -----
        No checking for consistency is done here. We do this at a later step
        in order to be able to keep the input methods as clean as posible.

        """

        incar_dict = {}
        for tag, value in incar.iteritems():
            # check for comment in value, if so skip shuffle to comment
            comment = None
            if isinstance(value, str):
                comment = value.split('#')
                if len(comment) > 2:
                    self._logger.info("Multiple comment tags detected for tag " +
                                     str(tag) + ". Exiting.")
                    sys.exit(1)
                if len(comment) == 1:
                    comment = None
            # create new instance of entry
            if comment is not None:
                if len(comment) == 2:
                    entry = IncarItem(tag, comment[0], comment[1])
                else:
                    entry = IncarItem(tag, value, comment)
            else:
                entry = IncarItem(tag, value, comment)
            clean_tag = entry.get_tag()
            if clean_tag in incar_dict:
                self._logger.info("Tag " + entry.get_tag() + " already "
                                 "found in the INCAR dictionary, "
                                 "overwriting it.")
            incar_dict[clean_tag] = entry

        return incar_dict
    
    def _convert_value_to_string(self, value):
        """Converts a value for an INCAR entry to a string that
        is compatible with VASP.


        Parameters
        ----------
        value : string
            The entry value of the INCAR entry.

        Returns
        -------
        string : string
            The entry value of the INCAR entry as a string.

        """

        # possible values are:
        # 1 - integer
        # True - bool
        # Something - a string
        # [Something, anotherthing] - list of strings
        # [1, 2, 3, 4] - list of integers
        # 1.0 - float
        # [1.0, 2.0, 3.0] - list of floats

        
        if type(value) is list:
            # list of values (we know all are either string, int or float)
            string = ' '.join(map(str, value))
        else:
            if isinstance(value, bool):
                if value:
                    return ".TRUE."
                else:
                    return ".FALSE."
            string = str(value)

        return string
    
    def validate(self):
        """Validate the content of the current Incar instance.

        Notes
        -----
        Uses a table that is synced with the VASP developers and thus gives
        an additional consistency check with respect to allowed parameters
        and their combinations.

        """

        return

    def modify(self, tag, value, comment=None):
        """Modify the entry tag in INCAR. If it is not found,
        add it.

        Parameters
        ----------
        tag : string
            The entry tag of the INCAR entry.
        value : string
            The entry value of the INCAR entry.
        comment : string
            The entry comment of the INCAR entry.

        """

        # create a new INCAR item and check it
        entry = IncarItem(tag, value, comment)
        # store or modify
        self.entries[entry.get_tag()] = entry

    def delete(self, tag):
        """Delete the entry with the supplied tag.

        Parameters
        ----------
        tag : string
            The entry tag of the INCAR entry.        

        """

        try:
            del self.entries[tag]
        except KeyError:
            pass

    def get(self, tag, comment=False):
        """Return the value and comment of the entry with tag.

        Parameters
        ----------
        tag : string
            The entry tag of the INCAR entry.        
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

    def get_dict(self):
        """Get a true dictionary containing the entries in an
        INCAR compatible fashion.

        Returns
        -------
        incar_dict : dict
            A dictionary on INCAR compatible form.

        """
        
        dictionary = {}
        for key, entry in self.entries.iteritems():
            dictionary[key] = entry.get_value()

        return dictionary

    def get_string(self):
        """Get a string containing the entries in an
        POSCAR compatible fashion. Each line is terminated by
        a newline character.

        Returns
        -------
        incar_string : string
            A string on POSCAR compatible form.

        """

        string_object = StringIO.StringIO()
        self._write(incar = string_object)
        incar_string = string_object.getvalue()
        string_object.close()

        return incar_string
        
        
    def write(self, file_path, comments=False):
        """Write the content of the current Incar instance to
        file.

        Parameters
        ----------
        file_path : string
            The location to write INCAR.
        comments : bool, optional
            If set to true, the comments are also dumped to the file,
            else not.

        """
        
        incar = utils.file_handler(file_path, status='w')        
        self._write(incar = incar)
        utils.file_handler(file_handler=incar)        
        
    def _write(self, incar, comments=False):
        """Write the content of the current Incar instance to
        file or string.

        Parameters
        ----------
        incar : object
            Either a file object or a StringIO object.
        comments : bool, optional
            If set to true, the comments are also dumped to the file,
            else not.

        """


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
            incar.write(string)

class IncarItem(object):

    def __init__(self, tag, value, comment, logger=None):
        """Initialize an entry in INCAR.

        Parameters
        ----------
        tag : string
            The entry tag of the INCAR entry.
        value : string
            The entry value of the INCAR entry.
        comment : string
            The entry comment of the INCAR entry.
        logger : object, optional
            A standard Python logger object.

        """

        # clean tag and value
        clean_tag, clean_value, clean_comment = self._clean_entry(tag,
                                                                  value,
                                                                  comment)
        self.tag = clean_tag
        self.value = clean_value
        self.comment = clean_comment
        
    def get_tag(self):
        return self.tag

    def get_value(self):
        return self.value

    def get_comment(self):
        return self.comment

    def _clean_entry(self, tag, value, comment):
        """Cleans the tag and value, for instance makes sures
        there are no spaces around the tag, that it is in lower case,
        and that the value is build into a list of several different
        parameters are given.

        Parameters
        ----------
        tag : string
            The entry tag of the INCAR entry.
        value : string
            The entry value of the INCAR entry.
        comment : string
            The entry comment of the INCAR entry.

        Returns
        ------
        clean_tag : string
            The entry tag of the INCAR entry, cleaned and checked.
        clean_value : string
            The entry value of the INCAR entry, cleaned and checked.
        clean_comment : string
            The entry comment of the INCAR entry, cleaned and checked.

        """

        # remove possible spaces on tag
        clean_tag = tag.split()

        # make sure tag is lowerscore
        clean_tag = clean_tag[0].lower()

        # possible solutions for the value are (if INCAR is read):
        # 1 - integer
        # .TRUE. - bool
        # Something - a string
        # Something anotherthing - set of strings
        # 1 2 3 4 - set of integers
        # 1.0 - float
        # 1.0 2.0 3.0 - set of floats

        # however, let us also open for the fact that users might
        # give a value what they would in INCAR

        # make sure we keep compatibility between Python 2 and 3
        try:
            basestring
        except NameError:
            basestring = str

        if isinstance(value, basestring):
            if clean_tag == "system":
                # if value is SYSTEM, treat it a bit special and
                # leave its string intact but remove grub
                clean_value = value.strip()
                return clean_tag, clean_value, None

            # values is a string, so we have read an INCAR or the user
            # tries to assign an entry with just a string
            values = value.split()

            # if we have some kind of set, check if all are ints, floats
            # or strings
            content_type = []
            clean_value = []
            for element in values:
                cnt_type = utils.test_string_content(element)
                content_type.append(cnt_type)
                if cnt_type == 'int':
                    cnt = int(element)
                elif cnt_type == 'float':
                    cnt = float(element)
                else:
                    # we have here also have a bool, so check that
                    cnt = self._test_string_for_bool(element)
                clean_value.append(cnt)
                
            # now if there is only one element in clean_value
            # remove list
            if len(clean_value) == 1:
                clean_value = clean_value[0]

            # check if all values are the same type (they should be)
            if not all(x == content_type[0] for x in content_type):
                self._logger.error("All values of the tag " + clean_tag.upper() +
                                  " are not of the same type. Most likely "
                                  "you forgot to add # as a comment tag. "
                                  "Exiting.")
                sys.exit(1)

        else:
            # if the user wants to assign an element as a list, int,
            # float or bool, accept this
            if not (isinstance(value, int) or isinstance(value, float)
                    or isinstance(value, bool) or (type(value) is list)):
                self._logger.error("The type: " + str(type(value)) + " of the "
                                  "supplied value for the INCAR tag is not "
                                  "recognized. Exiting.")
                sys.exit(1)
            clean_value = value

        # finally clean comment by removing any redundant spaces
        if comment is not None:
            clean_comment = comment.strip()
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
           or ".true." in string \
           or ".t." in string \
           or ".T." in string:
            return True
        elif ".False." in string \
             or ".FALSE." in string \
             or ".false." in string \
             or ".f." in string \
             or ".F." in string:
            return False
        else:
            return string
