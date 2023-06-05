"""
Standard stream parser for VASP.

--------------------------------
Contains parsers for the standard stream that originates from VASP. It fills a set with streams (e.g.
errors and warnings as defined in the streams.yml file).
"""

import re
import sys
from pathlib import Path

import yaml

from parsevasp import utils
from parsevasp.base import BaseParser


class Stream(BaseParser):  # pylint: disable=R0902
    """Class to handle standard stream."""

    def __init__(self, file_path=None, file_handler=None, logger=None, history=False, config=None):
        """Initialize an object that contain a standard stream composed of e.g. the standard output and error.

        Parameters
        ----------
        file_path : string, optional
            The file path that contains the standard stream.
        file_handler : object
            The file handler object.
        logger: object
            A logger object.
        stream : string, optional
            A string determining if a stdout, stderr or a combined stream is supplied. Not implemented.
        history : bool, optional
            If True, keep track of all the stream elements in appearing order.
        config : dict, optional
            A dictionary containing the override configuration of the recognized errors and warnings.
            Setting this will override the supplied error and warning configuration, or add new entries.

        """

        super().__init__(file_path=file_path, file_handler=file_handler, logger=logger)

        self._file_path = file_path
        self._file_handler = file_handler
        self._history = history
        self._streams = []
        self._inverse_streams = []
        self._config = config

        # Check that at least file path or file handler is supplied
        if self._file_path is None and self._file_handler is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        if self._file_path is None and self._file_handler is None:
            self._logger.error(self.ERROR_MESSAGES[self.ERROR_USE_ONE_ARGUMENT])
            sys.exit(self.ERROR_USE_ONE_ARGUMENT)

        # Load stream configuration from the supplied config or the standard config file
        self._stream_config = self._load_config()

        # Now investigate the kinds of streams present in the config
        self._stream_kinds = self._set_streams()

        # Build list of error and warning objects and store the these as stream triggers on which
        # we will react if detected in the stream.
        self._stream_triggers, self._inverse_stream_triggers = self._build_stream_triggers()

        # Parse parse parse
        self._parse()

        # Add inverse triggers not detected during parsing
        self._add_inverse_triggers()

    def __repr__(self):
        """Define representation to list number of streams found."""
        return f'StreamScanner found {len(self._streams)} streams'

    def __str__(self):
        """Define a string representation for the class which can be used for reporting purposes."""
        return f"We detected {self.number_of_entries} unique {', '.join(self._stream_kinds)}'s"

    def _write(self, *args, **kwargs):
        """Write not supported for standard streams."""
        raise NotImplementedError('Writing VASP standard streams files is not supported.')

    @property
    def configured_streams(self):
        """Return the configured streams."""
        return self._stream_triggers

    @property
    def kinds(self):
        """Return a list containing the different kinds."""
        return self._stream_kinds

    @property
    def entries(self):
        """Return the found streams after parsing as a list."""
        return self._streams

    @property
    def has_entries(self):
        """True if there are streams present after parsing."""
        return bool(self._streams)

    @property
    def number_of_entries(self):
        """Return a dict containing the number of unique streams detected."""
        return len(self._streams)

    def _load_config(self):
        """Load the configuration of the stream."""

        # First load the standard entries from file
        stream_config = self._load_config_from_file()
        # Then override or add new entries with the supplied entries.
        if self._config is not None:
            stream_config.update(self._config)

        return stream_config

    def _load_config_from_file(self):
        """Read the configuration of the errors and warnings from a yaml file and save it as the class method."""

        stream_config = None
        fname = Path(__file__).parent / 'stream.yml'
        # Read the config file
        with open(fname, 'r', encoding='utf8') as file_handler:
            stream_config = yaml.safe_load(file_handler)

        return stream_config

    def _set_streams(self):
        """Check the kinds of streams present in the config files."""
        stream_kinds = []
        for _, value in self._stream_config.items():
            kind = value['kind']
            if isinstance(kind, str):
                if kind.upper() in VaspStream._ALLOWED_STREAMS:  # pylint: disable=W0212
                    if not kind in stream_kinds:
                        stream_kinds.append(kind)
            else:
                raise ValueError('One of the kind entries is not a string.')
        return stream_kinds

    def _build_stream_triggers(self):
        """Here we use the stream configs to initialize the triggers"""

        # Define container for the triggers
        triggers = {}
        inverse_triggers = {}
        for stream in self._stream_kinds:
            triggers[''.join([stream.lower(), 's'])] = []
            inverse_triggers[''.join([stream.lower(), 's'])] = []

        for stream in self._stream_kinds:
            for shortname, config in self._stream_config.items():
                if config['kind'] == stream:
                    try:
                        inverse = config['inverse']
                    except KeyError:
                        inverse = False
                    if not inverse:
                        triggers[''.join([stream.lower(), 's'])].append(VaspStream(shortname=shortname, **config))
                    else:
                        inverse_triggers[''.join([stream.lower(),
                                                  's'])].append(VaspStream(shortname=shortname, **config))

        return triggers, inverse_triggers

    def _parse(self):
        """Perform the actual parsing."""

        if self._file_path is None and self._file_handler is None:
            return

        # Create dictionary from a file
        self._from_file()

    def _from_file(self):
        """Create a dictionary of entries from a
        file and store them in the this instance's data dictionary.

        """

        stream = utils.read_from_file(self._file_path, self._file_handler, encoding='utf8')
        self._from_list(stream)

    def _from_list(self, stream):
        """Go through the list and extract any recognized entries.

        Parameters
        ----------
        stream : list
            A list of strings containing each line in the standard stream.

        """
        # Make sure we access all triggers
        stream_triggers = {}
        for kind in self._stream_kinds:
            key = ''.join([kind.lower(), 's'])
            stream_triggers[key] = self._stream_triggers[key]
            stream_triggers[key].extend(self._inverse_stream_triggers[key])

        # Ignore list to avoid multiple occurrences if not requested (default)
        ignore = {}
        for kind, _ in stream_triggers.items():
            ignore[kind] = []
        for _, line in enumerate(stream):  # pylint: disable=too-many-nested-blocks
            # Go though all entries in the stream triggers
            for kind, triggers in stream_triggers.items():
                # Not check all the triggers of the given kind
                for index, trigger in enumerate(triggers):
                    if index not in ignore[kind]:
                        trigger_record = trigger.check_line(line)
                        if trigger_record is not None:
                            if not trigger_record.inverse:
                                # Store the streams that we trigger on and want to store.
                                self._streams.append(trigger_record)
                            else:
                                # Keep track of streams detected that we trigger on, but do not want stored.
                                self._inverse_streams.append(trigger_record.shortname)
                            if not self._history:
                                # Add index to avoid storing same trigger multiple times if we do not want the
                                # full history of streams (e.g. multiple stream occurrences recorded)
                                ignore[kind].append(index)

    def _add_inverse_triggers(self):
        """Adds the triggers that are marked as inverse and are not detected, meaning they should be
        included in the returned streams.

        """

        for _, triggers in self._inverse_stream_triggers.items():
            for trigger in triggers:
                if trigger.shortname not in self._inverse_streams:
                    # Only add inverse triggers if they was not detected during parsing.
                    self._streams.append(
                        VaspStream(
                            trigger.shortname, trigger.kind, trigger.regex, trigger.message, trigger.suggestion,
                            trigger.location, trigger.recover, trigger.inverse
                        )
                    )


class VaspStream:  # pylint: disable=R0902
    """Class representing stream elements given by VASP that we want to trigger on."""

    _ALLOWED_STREAMS = ['ERROR', 'WARNING']
    _ALLOWED_LOCATIONS = ['STDOUT', 'STDERR']

    def __init__(
        self,
        shortname,
        kind,
        regex,
        message,
        suggestion=None,
        location='STDOUT',
        recover=False,
        inverse=False
    ):  # pylint: disable=too-many-arguments
        """
        Initialise a VaspStream object.

        Paramters
        ---------
        shortname : string
            A short and unique string that identifies the stream
        kind : string
            The type of stream regex.
        regex : string
            Regex used for scanning.
        message : string
            Message to the user.
        suggestion : string, optional
            String containing a suggestion on how to address the stream message. Defaults to None.
        location : string, optional
            The location of the stream (typically STDOUT or STDERR). Defaults to STDOUT.
        recover : bool, optional
            True if the stream indicates that we are able to recover using some measures. Defaults to False.
        inverse : bool, optional
            If True, the stream should only be triggered on if we do not detect it.

        """
        self.shortname = shortname
        self.kind = kind
        if isinstance(regex, str):
            self.regex = re.compile(regex)
        else:
            self.regex = regex
        self.message = message
        self.suggestion = suggestion
        self.location = location
        self.recover = recover
        self.inverse = inverse

    def __repr__(self):
        """Set the representation."""
        return (
            f'VaspStream(kind={self.kind}, re={self.regex}, message={self.message}, '
            f'recover={self.recover}, inverse={self.inverse})'
        )

    def __str__(self):
        """Set string representation of the stream entry that can be used in a human readable report."""
        return f'({self.kind}) {self.shortname}: {self.message}'

    @property
    def shortname(self):
        """Return the shortname of the stream."""
        return self._shortname

    @shortname.setter
    def shortname(self, shrt):
        """Setter for the shortname that validates that it is a string."""

        if not isinstance(shrt, str):
            raise ValueError('The supplied shortname is not of type string.')
        self._shortname = shrt

    @property
    def kind(self):
        """Return the kind."""
        return self._kind

    @kind.setter
    def kind(self, knd):
        """Setter for kind that validates if the entries are supported."""
        if isinstance(knd, str):
            if knd.upper() in self._ALLOWED_STREAMS:
                self._kind = knd
            else:
                raise ValueError(
                    f'The type of kind for {self._shortname} is not supported. '
                    f'Currently we support {self._ALLOWED_STREAMS}.'
                )
        else:
            raise ValueError(f'The kind for {self._shortname} is not of type string.')

    @property
    def regex(self):
        """Return the regex."""
        return self._regex

    @regex.setter
    def regex(self, reg):
        """Setter for regex that validates and compiles if necessary."""
        if isinstance(reg, str):
            self._regex = re.compile(reg)
        else:
            self._regex = reg

    @property
    def message(self):
        """Return the message."""
        return self._message

    @message.setter
    def message(self, mes):
        """Setter for message that validates if it is a string."""
        if not isinstance(mes, str):
            raise ValueError('The message needs to be a string.')
        self._message = mes

    @property
    def suggestion(self):
        """Return the suggestion."""
        return self._suggestion

    @suggestion.setter
    def suggestion(self, sug):
        """Setter for the suggestion which validated if it is a string."""
        if sug is not None:
            # Allow None
            if not isinstance(sug, str):
                raise ValueError(f'The suggestion entry {sug} is not of type string.')
        self._suggestion = sug

    @property
    def location(self):
        """Return the location of the stream."""
        return self._location

    @location.setter
    def location(self, loc):
        """Setter for the location that validates if it is an allowed value."""
        if not loc in self._ALLOWED_LOCATIONS:
            raise ValueError(f'The location entry {loc} is not one of the allowed values {self._ALLOWED_LOCATIONS}')
        self._location = loc

    @property
    def recover(self):
        """Return the recover status."""
        return self._recover

    @recover.setter
    def recover(self, rec):
        """Setter for the recover that validates if it is a boolean."""
        if not isinstance(rec, bool):
            raise ValueError(f'The recover entry {rec} is not of type bool.')
        self._recover = rec

    @property
    def recoverable(self):
        """True if the stream is marked as recoverable."""
        return self._recover

    @property
    def inverse(self):
        """Return the inverse status."""
        return self._inverse

    @inverse.setter
    def inverse(self, inv):
        """Setter for the inverse that validates if it is a boolean."""
        if not isinstance(inv, bool):
            raise ValueError(f'The recover entry {inv} is not of type bool.')
        self._inverse = inv

    def check_line(self, line):
        """Check the stream in a line, return True the stream is found"""
        mch = self.regex.search(line)
        if mch:
            # Make a new instance for this particular error (in case we want
            # to save each and every error)
            return VaspStream(
                self.shortname, self.kind, self.regex, self.message, self.suggestion, self.location, self.recover,
                self.inverse
            )

        return None
