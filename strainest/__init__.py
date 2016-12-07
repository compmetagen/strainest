import os

__version__ = "1.2"

def mummer_path():
    return os.path.join(os.path.dirname(__file__), "MUMmer323")

import mummer
import api
