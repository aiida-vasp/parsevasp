"""Handle POTCAR

This parser takes parts of the
`pymatgen parser<https://github.com/materialsproject/pymatgen/blob/v2023.3.23/pymatgen/io/vasp/inputs.py#L1616-L2211>`_
"""

import re

from parsevasp import utils
from parsevasp.base import BaseParser


class Potcar(BaseParser):
    """Class to handle the POTCAR"""

    _functional_tags = {
        'pe': {
            'name': 'PBE',
            'class': 'GGA'
        },
        '91': {
            'name': 'PW91',
            'class': 'GGA'
        },
        'rp': {
            'name': 'revPBE',
            'class': 'GGA'
        },
        'am': {
            'name': 'AM05',
            'class': 'GGA'
        },
        'ps': {
            'name': 'PBEsol',
            'class': 'GGA'
        },
        'pw': {
            'name': 'PW86',
            'class': 'GGA'
        },
        'lm': {
            'name': 'Langreth-Mehl-Hu',
            'class': 'GGA'
        },
        'pb': {
            'name': 'Perdew-Becke',
            'class': 'GGA'
        },
        'ca': {
            'name': 'Perdew-Zunger81',
            'class': 'LDA'
        },
        'hl': {
            'name': 'Hedin-Lundquist',
            'class': 'LDA'
        },
        'wi': {
            'name': 'Wigner Interpoloation',
            'class': 'LDA'
        },
    }

    def __init__(self, file_path=None, file_handler=None, logger=None):
        super().__init__(file_path=file_path, file_handler=file_handler, logger=logger)

        self.metadata = None
        self._symbol = None
        self._element = None

        if self._file_path is not None or self._file_handler is not None:
            self._from_file()
        if self._file_path is None and self._file_handler is None:
            raise ValueError('Either "file_path" or "file_handler" should be passed')

    def _from_file(self):
        """Create rudimentary dictionary of entries from a
        file.

        """

        potcar = utils.read_from_file(self._file_path, self._file_handler, encoding='utf8', lines=False)

        return self._generate_metadata(potcar)

    def _generate_metadata(self, potcar_contents: str):
        """Get the metadata from a POTCAR file

        Parameters
        ----------
        potcar_contents: string
            The contents of the POTCAR file as a string

        Returns
        -------
        metadata: dictionary
            A dictionary containing the metadata associated with the POTCAR
        """
        _parameters_to_parse = {
            'VRHFIN': lambda val: val.strip(),
            'LEXCH': lambda val: val.strip(),
            'TITEL': lambda val: val.strip(),
            'LULTRA': lambda val: re.match(r'^\.?([TFtf])[A-Za-z]*\.?', val).group(1).lower() in ['t'],
            'LCOR': lambda val: re.match(r'^\.?([TFtf])[A-Za-z]*\.?', val).group(1).lower() in ['t'],
            'LPAW': lambda val: re.match(r'^\.?([TFtf])[A-Za-z]*\.?', val).group(1).lower() in ['t'],
            'IUNSCR': lambda val: int(re.match(r'^-?[0-9]+', val).group(0)),
            'NDATA': lambda val: int(re.match(r'^-?[0-9]+', val).group(0)),
            'ICORE': lambda val: int(re.match(r'^-?[0-9]+', val).group(0)),
            'EATOM': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'RPACOR': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'POMASS': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'ZVAL': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'RCORE': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'RWIGS': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'ENMAX': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'ENMIN': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'RCLOC': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'EAUG': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'DEXC': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'RMAX': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'RAUG': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'RDEP': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'RDEPT': lambda val: float(re.search(r'^-?\d*\.?\d*[eE]?-?\d*', val).group(0)),
            'STEP': lambda val: [float(y) for y in re.split(r'\s+', val.strip()) if not y.isalpha()],
        }
        search_lines = re.search(
            r'(?s)(parameters from PSCTR are:'
            r'.*?END of PSCTR-controll parameters)',
            potcar_contents,
        ).group(1)

        self.metadata = {}
        for key, val in re.findall(r'(\S+)\s*=\s*(.*?)(?=;|$)', search_lines, flags=re.MULTILINE):
            if key in _parameters_to_parse:
                self.metadata[key] = _parameters_to_parse[key](val)

        try:
            self._symbol = self.metadata['TITEL'].split(' ')[1].strip()
        except IndexError:
            self._symbol = self.metadata['TITEL'].strip()

        self._element = self._symbol.split('_')[0]

    @property
    def symbol(self):
        """
        Get the symbol associated with this POTCAR

        Returns
        -------
        symbol: string
            The POTCAR symbol, e.g. W_pv
        """
        return self._symbol

    @property
    def element(self):
        """
        Get the symbol of the element associated with this POTCAR

        Returns
        -------
        element: string
            The POTCAR element, e.g. W
        """
        return self._element

    @property
    def functional(self):
        """
        Get the functional associated with this POTCAR

        Returns
        -------
        functional: string
            Functional associated with the current POTCAR file.
        """
        return self._functional_tags.get(self.metadata.get('LEXCH').lower(), {}).get('name')

    @property
    def functional_class(self):
        """
        Get the functional class associated with this POTCAR

        Returns
        --------
        functional_class: string
            Functional class associated with the current POTCAR file.
        """
        return self._functional_tags.get(self.metadata.get('LEXCH').lower(), {}).get('class')

    def __getattr__(self, attribute):
        """
        Delegates attributes to keywords. For example, you can use
        potcar.enmax to get the ENMAX of the POTCAR.
        For float type properties, they are converted to the correct float. By
        default, all energies in eV and all length scales are in Angstroms.
        """
        try:
            return self.metadata[attribute.upper()]
        except Exception as exc:
            raise AttributeError(attribute) from exc

    def _write(self, file_handler, **kwargs):
        pass
