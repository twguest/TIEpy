# -*- coding: utf-8 -*-

import numpy as np

from PIL import Image

def load_tif(directory):
    """
    load a tif file as an np array
    """
    return (np.asarray(Image.open(directory)))