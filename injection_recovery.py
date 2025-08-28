from astropy.io import fits
import glob, os, sys
from datetime import datetime
from pathlib import Path
import shutil

import subprocess


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count

from getGravityObs import *

#def createSyntheticCompanion(wav, visFiber, uCoord, vCoord, visOnStar, contrast, deltaRA, deltaDec):
def createSyntheticCompanion(dataset, visOnStar, contrast, deltaRA, deltaDec):
    ## synthetic Companioin signal
    ## as per Pourré, N., et al.: A&A, 686, A258 (2024):
    ## Vsyn,comp = C*Vonstar * exp[(-i2pi/lambda)(OPD)]exp[i phi] : need to come back to check on if this phase term is properly accounted for!!!!!

    phaseShift = np.exp(-1j * ((deltaRA * dataset['ucoord']) + (deltaDec * dataset['vcoord'])))
    visFiber = getVisPhased(dataset)
    visSyntheticComp = visFiber + (contrast * visOnStar * phaseShift) * np.exp(np.angle(visOnStar))

    return getVisUnphased(visSyntheticComp, dataset)
        
def saveSyntheticReducedData(planetFile, visSyntheticComp, i, odir):
    
    with fits.open(planetFile, mode='update', memmap=False) as hdul:
        for hdu in hdul:
            if hdu.header.get('EXTNAME') == 'OI_VIS' and hdu.header.get('extver') == 10:
                ndit, nbase, nwave = np.shape(visSyntheticComp)
                reshapedVISDATA = visSyntheticComp.reshape(-ndit * nbase, nwave)
                
                hdu.data['VISDATA'][:] = reshapedVISDATA
                break
    
        hdul.writeto(planetFile, overwrite=True)
    
    return odir
    
def injectionRecovery(args):
    
    folderName, planetFiles, visOnStar, refPhase, injectionParams, outPath, starFiles = args
    i = int(folderName)

    deltaRA_mas = injectionParams.deltaRA_mas.values[i]
    deltaDec_mas = injectionParams.deltaDec_mas.values[i]
    contrast = injectionParams.Contrast.values[i]
    
    
    for planetFile in planetFiles:
        dataset = oiVis(planetFile, 10)
        
#        print(f'Fiber Located at RA {dataset["SXY"][0]}, Dec {dataset["SXY"][1]} (mas)')
        fiberRA = dataset['SXY'][0]#/1000.0/3600.0/180.0*np.pi
        fiberDec = dataset['SXY'][1]#/1000.0/3600.0/180.0*np.pi
        
        injectionRA = fiberRA + deltaRA_mas
        injectionDec = fiberDec + deltaDec_mas

        syntheticVISDATA = createSyntheticCompanion(dataset, visOnStar,
            contrast, injectionRA, injectionDec
        )

        lastOutPath = saveSyntheticReducedData(planetFile, syntheticVISDATA, i, outPath)

    subprocess.run(
        ['python', '/Users/svach/gravi_tools3/run_gravi_astrored_check.py'],
        cwd=outPath,
        check=True
    )
    subprocess.run(
        ['python', '/Users/svach/gravi_tools3/run_gravi_astrored_astrometry.py', '--reDo=TRUE'],
        cwd=outPath,
        check=True
    )
#    os.system('python /Users/svach/gravi_tools3/run_gravi_astrored_check.py')
#    os.system('python /Users/svach/gravi_tools3/run_gravi_astrored_astrometry.py --reDo=TRUE')
            
def run(filePaths, resetPoint=8, injectionPath = '/Users/svach/gravityInjectionRecovery/injection_parameters.csv'):
#    sys.stdout = open(os.devnull, 'w')

    injectionParams = pd.read_csv(injectionPath)[resetPoint:]

    tasks = []
    for filePath in filePaths:
        odir = filePath + 'synthetic_reduced/'
        
        
        os.makedirs(odir, exist_ok=True)
        reducedFiles = filePath + '*_astroreduced.fits'
        
        originalStarFiles, originalPlanetFiles = getFiles(reducedFiles)
        for folderName in injectionParams.folderName.values:
        
            outPath = f"{odir}{folderName}/"
            os.makedirs(outPath, exist_ok=True)
            for originalStarFile in originalStarFiles:
                shutil.copy2(originalStarFile, outPath)
            for originalPlanetFile in originalPlanetFiles:
                shutil.copy2(originalPlanetFile, outPath)
                
            reducedFiles = outPath + '*_astroreduced.fits'
            starFiles, planetFiles = getFiles(reducedFiles)
            
            visOnStar, refPhase = getVonStar(starFiles)
        
            args = [
                folderName, planetFiles, visOnStar, refPhase, injectionParams, outPath, starFiles
            ]
            tasks.append(args)

        # Run tasks in parallel
    nproc = min(cpu_count(), len(tasks))
    with Pool(processes=nproc) as pool:
        pool.map(injectionRecovery, tasks)


if __name__ == '__main__':
    
    filePaths = [
#                '/Users/svach/GravityGaia/p115/2025-08-10/SCI_HD 14082B_search1/reduced/',
#                        '/Users/svach/GravityGaia/p115/2025-08-09/SCI_HD 218396_search2/reduced/',
                            '/Users/svach/GravityGaia/p115/2025-08-09/SCI_G 80-21_search_cc/reduced/',
#                            '/Users/svach/GravityGaia/p115/2025-06-12/SCI_HD 111588_search3/reduced/'
 	]
    run(filePaths, resetPoint=0)


