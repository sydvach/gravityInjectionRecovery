#
#  recovery_maps.py
#
#
#  Created by Sydney Vach on 19.08.25.
#

import numpy as np
import pandas as pd
from astropy.io import fits
import glob

if __name__ == '__main__':
    
    path = '/Users/svach/GravityGaia/p115/2025-08-10/SCI_HD 14082B_search1/reduced/synthetic_reduced/*/*fast*.fits'
    
    fileNames = glob.glob(path)
    fileNames = fileNames[0:1]
    print(fileNames)
    for fileName in fileNames:
        with fits.open(fileName) as hdul:
            print(hdul[0].header)
