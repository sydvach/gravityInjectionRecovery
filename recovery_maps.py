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


def calculateRecoveryMaps(path, injectedParams):
    
    simPath = path + 'reduced/synthetic_reduced/*/*fast*.fits'
    
    fileNames = sorted(glob.glob(simPath))
#    print(fileNames)
#    breakpoint()
#    fileNames = fileNames[0:1]
    
#    results = {}
    poly = {
            '1': np.zeros(len(injectedParams)),
            '2': np.zeros(len(injectedParams)),
            '3': np.zeros(len(injectedParams)),
            'index': np.arange(len(injectedParams))
        }
    for fileName in fileNames:
        
        injectionIndex = int(fileName.split('/')[-2])
        with fits.open(fileName) as hdul:
            fiberRA = hdul[0].header['X_FIBER']
            fiberDec = hdul[0].header['Y_FIBER']
            injectedRA = fiberRA + injectedParams.deltaRA_mas.values[injectionIndex]
            injectedDec = fiberDec + injectedParams.deltaDec_mas.values[injectionIndex]
            for n in range(1,4):
                offset = np.sqrt((injectedRA - hdul[0].header[f'X_{n}'])**2 + (injectedDec - hdul[0].header[f'Y_{n}'])**2)
                if offset <= 3:
                    deltaContrast = np.abs(injectedParams.Contrast.values[injectionIndex] - hdul[0].header[f'CONTRA_{n}']) / injectedParams.Contrast.values[injectionIndex]
                    if deltaContrast <= 4: # this needs to be corrected
                        poly[str(n)][injectionIndex] = 1.
                        
#        results[fileName.split('/')[-2] + '/'+ fileName.split('/')[-1]] = poly
        
#    print(results)
    df = pd.DataFrame(poly)
    df.to_csv(fileName.split('/')[6]+f'_RA-{fiberRA}_Dec-{fiberDec}_recovery_results.csv', index=False)
#    return poly
    
                    
                        

            

if __name__ == '__main__':
    
    injectedParams = pd.read_csv('~/gravityInjectionRecovery/injection_parameters.csv')
    path = '/Users/svach/GravityGaia/p115/2025-08-09/SCI_G 80-21_search_cc/'
    calculateRecoveryMaps(path, injectedParams)
#    breakpoint()
    
    
