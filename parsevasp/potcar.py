"""Handle POTCAR"""

import re

from parsevasp import utils
from parsevasp.base import BaseParser


class Potcar(BaseParser):

    def __init__(self, file_path=None, file_handler=None, logger=None):
        super().__init__(file_path, file_handler, logger)

    if self._file_path is not None or self._file_handler is not None:
        self.metadata = self._from_file()



    def _from_file(self):
        """Create rudimentary dictionary of entries from a
        file.

        """

        potcar = utils.read_from_file(self._file_path, self._file_handler, encoding='utf8')

        return self._generate_metadata(potcar)

    def _generate_metadata(potcar_contents:str) -> dict:
        _parameters_to_parse = {
            'VRHFIN': str,
            'LEXCH': str,
            'EATOM': float,
            'TITEL': str,
            'LULTRA': bool,
            'IUNSCR': int,
            'RPACOR': float,
            'POMASS': float, 'ZVAL': float,
            'RCORE': float,
            'RWIGS': float,
            'ENMAX': float, 'ENMIN': float,
            'RCLOC': float,
            'LCOR': bool,
            'LPAW': bool,
            'EAUG': float,
            'DEXC': float,
            'RMAX': float,
            'RAUG': float,
            'RDEP': float,
            'RDEPT': float,
            'NDATA': int,
            'STEP': list,
        }
        search_lines = re.search(
            r"(?s)(parameters from PSCTR are:"
            r".*?END of PSCTR-controll parameters)",
            potcar_contents,
        ).group(1)

        metadata = {}
        for key, val in re.findall(r"(\S+)\s*=\s*(.*?)(?=;|$)", search_lines, flags=re.MULTILINE):
            if key in _parameters_to_parse:
                if isinstance(_parameters_to_parse[key], str):
                    metadata[key] = val.split()
                if isinstance(_parameters_to_parse, int):
                    metadata[key] = int(re.match(r"^-?[0-9]+", val).group(0))
                if isinstance(_parameters_to_parse[key], float):
                    metadata[key] = float(re.search(r"^-?\d*\.?\d*[eE]?-?\d*", val).group(0))
                if isinstance(_parameters_to_parse[key], bool):
                    metadata[key] = re.match(r"^\.?([TFtf])[A-Za-z]*\.?", val).group(1).lower() in ['t']
                if isinstance(_parameters_to_parse[key], list):
                    metadata[key] = [float(y) for y in re.split(r"\s+", val.strip()) if not y.isalpha()]
        return metadata
