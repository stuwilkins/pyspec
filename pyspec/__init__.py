# Enforce proper version requirements
import sys
from spec import SpecDataFile
from spec import SpecScan

if sys.version[0:3] < '2.5':
    raise ImportError('Python Version 2.5 or above is required for pyspec.')

