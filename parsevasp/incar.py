"""Handle INCAR."""
# pylint: disable=consider-using-f-string
import io
import logging
import sys

from parsevasp import constants, utils
from parsevasp.base import BaseParser


class Incar(BaseParser):
    """Class to handle INCAR."""

    ERROR_UNSUPPORTED_TAG = 100
    ERROR_TWO_EQUALS = 101
    ERROR_INVALID_COMMENT_SIGN = 102
    ERROR_MULTIPLE_COMMENTS = 103
    BaseParser.ERROR_MESSAGES.update({
        ERROR_TWO_EQUALS:
        'Detected two equal signs for an entry in the INCAR file.',
        ERROR_INVALID_COMMENT_SIGN:
        'Detected a comment line that does not start '
        'with a #. Please correct and be consistent.',
        ERROR_MULTIPLE_COMMENTS:
        'Multiple comment tags detected.',
        ERROR_UNSUPPORTED_TAG:
        'The supplied INCAR tag is not '
        'officially supported. Please consult the VASP manual or '
        'set the validate_tags attribute for the Incar class initializer to '
        'False if you want to disable tag checking.'
    })
    ERROR_MESSAGES = BaseParser.ERROR_MESSAGES

    def __init__(
        self,
        incar_string=None,
        incar_dict=None,
        file_path=None,
        file_handler=None,
        logger=None,
        prec=None,
        validate_tags=True
    ):  # pylint: disable=too-many-arguments
        """Initialize an INCAR object and set content as a dictionary.

        Parameters
        ----------
        incar_string : string, optional
            A string containing INCAR entries. Must contain line
            breaks if multiline, otherwise the INCAR will be mangled.
        incar_dict : dict, optional
            A dictionary containing the INCAR entries.
        prec : int, optional
            An integer describing how many decimals the users wants
            when printing files.
        validate_tags : bool, optional
            If True, validate the tags supplied against the VASP documentation.

        """

        super().__init__(file_path=file_path, file_handler=file_handler, logger=logger)

        self._incar_dict = incar_dict
        self._incar_string = incar_string
        self._validate_tags = validate_tags

        # Set precision
        if prec is None:
            self._prec = 12
        else:
            self._prec = prec
        self._width = self._prec + 4

        # Check that only one argument is supplied
        # pylint: disable=R0916
        if (self._incar_string is not None and self._incar_dict is not None) or (
            self._incar_string is not None and self._file_path is not None
        ) or (self._incar_dict is not None and self._file_path is not None) and self._file_handler is not None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        # Check that at least one is supplied
        if (
            self._incar_string is None and self._incar_dict is None and self._file_path is None and
            self._file_handler is None
        ):
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        if self._file_path is not None or self._file_handler is not None:
            # Create list from a file
            incar_list = self._from_file()

        if self._incar_string is not None:
            # Create list from a string
            incar_list = self._from_string()

        if self._incar_dict is None:
            incar = self._from_list(incar_list)
        else:
            incar = self._from_dict(self._incar_dict)

        # Ctore entries
        self.entries = incar

        # Validate dictionary
        self.validate()

    def _from_file(self):
        """Create rudimentary dictionary of entries from a
        file.

        """

        incar = utils.read_from_file(self._file_path, self._file_handler, encoding='utf8')
        return incar

    def _from_string(self):
        """Create rudimentary dictionary of entries from a
        string.

        """

        incar = self._incar_string.splitlines()
        return incar

    def _from_list(self, incar):
        """
        Go through the list and analyze for = and ; in order to
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
            # Check for comment at the start of a line, if so skip
            comment = not line.split('#')[0].split()
            if not comment:
                # Now check if a line contains a comment, if so, only take
                # what is in front of this, truly...
                splitted = line.split('#', 1)
                if len(splitted) > 1:
                    comment = splitted[1]
                else:
                    comment = None
                # Now split on ; as we could have combined entries
                splitted = splitted[0].split(';')
                for ntry in splitted:
                    if ntry == '\n':
                        # Skip if the user used ; at the end of a line
                        continue
                    # Then split on = and analyze each entry
                    final_split = ntry.split('=')
                    if len(final_split) > 2:
                        self._logger.error(
                            '{} The following line contains the problem:\n\n {}'
                            '\n\nPlease correct. Exiting.'.format(self.ERROR_MESSAGES[self.ERROR_TWO_EQUALS], ntry)
                        )
                        sys.exit(self.ERROR_TWO_EQUALS)
                    if len(final_split) == 1:
                        self._logger.error(self.ERROR_MESSAGES[self.ERROR_INVALID_COMMENT_SIGN])
                        sys.exit(self.ERROR_INVALID_COMMENT_SIGN)
                    tag = final_split[0]
                    value = final_split[1]
                    # Create new instance of entry
                    entry = IncarItem(tag, value, comment, logger=self._logger)
                    clean_tag = entry.get_tag()
                    if clean_tag in incar_dict:
                        self._logger.info(
                            f'Tag {entry.get_tag()} already found in the INCAR dictionary, overwriting it.'
                        )
                    incar_dict[clean_tag] = entry

        return incar_dict

    def _from_dict(self, incar):
        """"
        Go through the dict and setup the IncarItem elements.

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
        for tag, value in incar.items():
            # Check for comment in value, if so skip shuffle to comment
            comment = None
            if isinstance(value, str):
                comment = value.split('#')
                if len(comment) > 2:
                    self._logger.info(
                        f'{self.ERROR_MESSAGES[self.ERROR_MULTIPLE_COMMENTS]} The tag {str(tag)} is affected.'
                    )
                    sys.exit(self.ERROR_MULTIPLE_COMMENTS)
                if len(comment) == 1:
                    comment = None
            # Create new instance of entry
            if comment is not None:
                if len(comment) == 2:
                    entry = IncarItem(tag, comment[0], comment[1], logger=self._logger)
                else:
                    entry = IncarItem(tag, value, comment, logger=self._logger)
            else:
                entry = IncarItem(tag, value, comment, logger=self._logger)
            clean_tag = entry.get_tag()
            if clean_tag in incar_dict:
                self._logger.info(f'Tag {entry.get_tag()} already found in the INCAR dictionary, overwriting it.')
            incar_dict[clean_tag] = entry

        return incar_dict

    def _convert_value_to_string(self, value):
        """
        Converts a value for an INCAR entry to a string that
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

        # Possible values are:
        # 1 - integer
        # True - bool
        # Something - a string
        # [Something, anotherthing] - list of strings
        # [1, 2, 3, 4] - list of integers
        # 1.0 - float
        # [1.0, 2.0, 3.0] - list of floats

        if isinstance(value, list):
            # List of values (we know all are either string, int or float)
            string = ' '.join(map(str, value))
        else:
            if isinstance(value, bool):
                if value:
                    return '.TRUE.'
                return '.FALSE.'
            string = str(value)

        return string

    def validate(self):
        """
        Validate the content of the current Incar instance.

        Notes
        -----
        Uses a table that is synced with the VASP developers and thus gives
        an additional consistency check with respect to allowed parameters.
        Currently, only the tag name is checked, not its parameters.

        """
        # If we do not want to validate (e.g. if you supply a set which contains
        # unsupported keys, return now).
        if not self._validate_tags:
            return

        allowed_keys = constants.incar_tags.keys()
        for key, _ in self.entries.items():
            if key not in allowed_keys:
                self._logger.error(
                    f'{self.ERROR_MESSAGES[self.ERROR_UNSUPPORTED_TAG]} The tag in question is {key.upper()}'
                )
                sys.exit(self.ERROR_UNSUPPORTED_TAG)

        return

    def modify(self, tag, value, comment=None):
        """
        Modify the entry tag in INCAR. If it is not found, add it.

        Parameters
        ----------
        tag : string
            The entry tag of the INCAR entry.
        value : string
            The entry value of the INCAR entry.
        comment : string
            The entry comment of the INCAR entry.

        """

        # Create a new INCAR item and check it
        entry = IncarItem(tag, value, comment, logger=self._logger)
        # Store or modify
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
        """
        Return the value and comment of the entry with tag.

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
        return value

    def get_dict(self):
        """
        Get a true dictionary containing the entries in an
        INCAR compatible fashion.

        Returns
        -------
        incar_dict : dict
            A dictionary on INCAR compatible form.

        """

        dictionary = {}
        for key, entry in self.entries.items():
            dictionary[key] = entry.get_value()

        return dictionary

    def get_string(self):
        """
        Get a string containing the entries in an
        POSCAR compatible fashion. Each line is terminated by
        a newline character.

        Returns
        -------
        incar_string : string
            A string on POSCAR compatible form.

        """

        string_object = io.StringIO()
        self._write(file_handler=string_object)
        incar_string = string_object.getvalue()
        string_object.close()

        return incar_string

    def _write(self, file_handler, **kwargs):
        """
        Write the content of the current Incar instance to
        file or string.

        Parameters
        ----------
        file_handler : object
            Either a file object or a StringIO object.
        comments : bool, optional
            If set to true, the comments are also dumped to the file,
            else not.

        """

        comments = kwargs.pop('comments', False)

        # Write in alfabetical order
        keys = sorted(self.entries)
        entries = self.entries
        for key in keys:
            entry = entries[key]
            value = entry.get_value()
            comment = entry.get_comment()
            if comment is None or not comments:
                comment = ''
            else:
                comment = ' # ' + comment
            value = self._convert_value_to_string(value)
            string = str(key.upper()) + ' = ' + value + comment + '\n'
            file_handler.write(string)


class IncarItem:
    """Class to treat each entry in INCAR."""

    ERROR_VALUES_NOT_SAME_TYPE = 104
    ERROR_INVALID_TYPE = 105

    def __init__(self, tag, value, comment, logger=None):
        """
        Initialize an entry in INCAR.

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

        # Set logger
        if logger is not None:
            self._logger = logger
        else:
            logging.basicConfig(level=logging.DEBUG)
            self._logger = logging.getLogger('IncarItem')

        # Clean tag and value
        clean_tag, clean_value, clean_comment = self._clean_entry(tag, value, comment)
        self.tag = clean_tag
        self.value = clean_value
        self.comment = clean_comment

    def get_tag(self):
        """Return tag."""
        return self.tag

    def get_value(self):
        """Return value."""
        return self.value

    def get_comment(self):
        """Return comment."""
        return self.comment

    def _clean_entry(self, tag, value, comment):
        """
        Cleans the tag and value, for instance makes sures
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

        # Remove possible spaces on tag
        clean_tag = tag.split()

        # Make sure tag is lowerscore
        clean_tag = clean_tag[0].lower()

        # Possible solutions for the value are (if INCAR is read):
        # 1 - integer
        # .TRUE. - bool
        # Something - a string
        # Something anotherthing - set of strings
        # 1 2 3 4 - set of integers
        # 1.0 - float
        # 1.0 2.0 3.0 - set of floats
        # 10*1.0 etc., interpreted as strings here

        # However, let us also open for the fact that users might
        # give a value what they would in INCAR

        if isinstance(value, str):
            if clean_tag in ('system', 'magmom', 'm_constr'):
                # If value is SYSTEM or MAGMOM (can contain asterix), treat it a bit special and
                # leave its string intact but remove grub
                clean_value = value.strip()
                return clean_tag, clean_value, None

            # Values is a string, so we have read an INCAR or the user
            # tries to assign an entry with just a string
            values = value.split()

            # If we have some kind of set, check if all are ints, floats
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
                    # We have here also have a bool, so check that
                    cnt = self._test_string_for_bool(element)
                clean_value.append(cnt)

            # Now if there is only one element in clean_value
            # remove list
            if len(clean_value) == 1:
                clean_value = clean_value[0]

            # Check if all values are the same type (they should be)
            if not all(x == content_type[0] for x in content_type):
                self._logger.error(
                    'All values of an INCAR tag are not of the same type. '
                    'Maybe you forgot to add # as a comment tag?'
                    f' The tag in question is: {clean_tag.upper()}'
                )
                sys.exit(self.ERROR_VALUES_NOT_SAME_TYPE)

        else:
            # If the user wants to assign an element as a list, int,
            # float or bool, accept this, including unicode
            if not isinstance(value, (bool, float, int, list)):
                self._logger.error(
                    'The type one of the supplied values for the INCAR tag '
                    f'is not recognized. The tag in question is: {clean_tag.upper()}'
                )
                sys.exit(self.ERROR_INVALID_TYPE)
            clean_value = value

        # Finally clean comment by removing any redundant spaces
        if comment is not None:
            clean_comment = comment.strip()
        else:
            clean_comment = None

        return clean_tag, clean_value, clean_comment

    def _test_string_for_bool(self, string):
        """
        Detects if string contains Fortran bool.

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
        if '.True.' in string \
           or '.TRUE.' in string \
           or '.true.' in string \
           or '.t.' in string \
           or '.T.' in string:
            return True
        if '.False.' in string \
             or '.FALSE.' in string \
             or '.false.' in string \
             or '.f.' in string \
             or '.F.' in string:
            return False
        return string
