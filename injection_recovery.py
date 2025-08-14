from astropy.io import fits
import glob, os
from datetime import datetime
from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count

from injection_params import make_injection_params

mas2rad = np.pi / (180.0 * 3600.0 * 1000.0)

def getFiles(filePath):
    #from exoGRAVITY pipeline: https://gitlab.obspm.fr/mnowak/exogravity/-/blob/master/tutorial/EX1_data_manipulation.ipynb?ref_type=heads
    fileNames = sorted(glob.glob(filePath))
    
    starFiles, planetFiles = [], []
    for fileName in fileNames:
        with fits.open(fileName) as hdul:
            hdr = hdul[0].header
            sObjX, sObjY = hdr["HIERARCH ESO INS SOBJ X"], hdr["HIERARCH ESO INS SOBJ Y"]
#        print("File {}:\n  Fiber offset: dRA = {} mas and dDEC = {} mas".format(filename, sObjX, sObjY))
            if (sObjX**2+sObjY**2) < 1:
                starFiles.append(fileName)
            else:
                planetFiles.append(fileName)
    
    # print("On star files:\n", "\n ".join(starFiles))
    # print("On planet files:\n", "\n ".join(planetFiles))
    return fileNames, starFiles, planetFiles

def assignStar_to_Planets(allFiles, starFiles, planetFiles):
    starSet, planetSet = set(starFiles), set(planetFiles)
    
    assigned = []
    lastStar = None
    for f in allFiles:
        if f in starSet:
            lastStar = f
        elif f in planetSet:
            assigned.append((f, lastStar if lastStar is not None else None))
    return assigned

def getWavelengthGrid(filename):
    # extract the wavelength grid
    oi_wav = fits.getdata(filename, extname = "OI_WAVELENGTH", extver = 10)
    wav = oi_wav.field("EFF_WAVE")

    return wav.copy(), len(wav)

def getVisPhased(oi_vis, nbase, nwav, wav):
    #from exoGRAVITY pipeline: 
    #https://gitlab.obspm.fr/mnowak/exogravity/-/blob/master/tutorial/EX1_data_manipulation.ipynb?ref_type=heads
    
    # also in Section 10.26.1, page 84 of the pipeline manual
    #https://ftp.eso.org/pub/dfs/pipelines/instruments/gravity/gravity-pipeline-manual-1.9.0.pdf
    visData = oi_vis.field("VISDATA").reshape([-1, nbase, nwav]).copy()
    opdDisp = oi_vis.field("OPD_DISP").copy()
    phaseTelFc = oi_vis.field("PHASE_MET_TELFC").copy()
    phaseRef = oi_vis.field("PHASE_REF").copy()
    phaseFt = (phaseRef - 2*np.pi/wav*opdDisp - phaseTelFc).reshape([-1, nbase, nwav])
    visPhased = visData * np.exp(1j*phaseFt)
    return visPhased

def openReducedPlanet(planetFile):
    # extract the visdata and reshape it to a 3D array
    nbase = 6 # the number of baselines in GRAVITY
    oi_vis = fits.getdata(planetFile, extname = "OI_VIS", extver = 10)
    visData = oi_vis.field("VISDATA")
    wav, nwav = getWavelengthGrid(planetFile)
    visPhased = getVisPhased(oi_vis, nbase, nwav, wav)
    
    uCoord = oi_vis.field("UCOORD").reshape([-1, nbase]).copy()
    vCoord = oi_vis.field("VCOORD").reshape([-1, nbase]).copy()

    with fits.open(planetFile) as hdul:
        hdr = hdul[0].header
        sObjX, sObjY = hdr["HIERARCH ESO INS SOBJ X"], hdr["HIERARCH ESO INS SOBJ Y"]

    return wav, visPhased, uCoord, vCoord, sObjX, sObjY

def openReducedStar(starFile, visPlanet):
    nbase = 6 # the number of baselines
    oi_vis = fits.getdata(starFile, extname = "OI_VIS", extver = 10)
    wav, nwav = getWavelengthGrid(starFile)
    visPhased = getVisPhased(oi_vis, nbase, nwav, wav)

    #correct onstar shape
    visOnStar = np.mean(visPhased, axis=0)
    visOnStar = np.repeat(visOnStar[None, :, :], visPlanet.shape[0], axis=0)
    refPhase = np.angle(visOnStar)
    return visOnStar, refPhase

def createSyntheticCompanion(wav, visPointing, uCoord, vCoord, visOnStar, contrast, deltaRA, deltaDec):
    ## synthetic Companioin signal
    ## as per Pourré, N., et al.: A&A, 686, A258 (2024):
    ## Vsyn,comp = C*Vonstar * exp[(-i2pi/lambda)(OPD)]exp[i phi] : need to come back to check on if this phase term is properly accounted for!!!!!

    opd = (uCoord * deltaRA) + (vCoord * deltaDec)
    phaseShift = np.exp(-2j * np.pi * opd[:, :, None]  / wav[None, None, :])
    visSyntheticComp = visPointing + contrast * visOnStar * phaseShift
    
    return visSyntheticComp
        
def saveSyntheticReducedData(planetFile, starFile, visSyntheticComp, i, odir):
    outPath = f"{odir}{int(i)}/"
    os.makedirs(outPath, exist_ok=True)
    
    with fits.open(planetFile, mode='update', memmap=False) as hdul:
        for hdu in hdul:
            if hdu.header.get('EXTNAME') == 'OI_VIS':
                ndit, nbase, nwave = np.shape(visSyntheticComp)
                reshapedVISDATA = visSyntheticComp.reshape(ndit * nbase, nwave)
                hdu.data['VISDATA'][:] = reshapedVISDATA
                break
    
        hdul.writeto(outPath + Path(planetFile).stem + '.fits', overwrite=True)
    
    shutil.copy2(starFile, outPath)
    return outPath
    
def injectionRecovery(args):
    folderName, files, injectionParams, odir = args
    i = int(folderName)
    
    deltaRA_rad = injectionParams.deltaRA_rad.values[i]
    deltaDec_rad = injectionParams.deltaDec_rad.values[i]
    contrast = injectionParams.Contrast.values[i]
    
    lastOutPath = []
    
    for planetFile, starFile in files:
        wav, visPhased, uCoord, vCoord, fiberRA, fiberDec = openReducedPlanet(planetFile)
        fiberRA *= mas2rad
        fiberDec *= mas2rad
        
        visOnStar, refPhase = openReducedStar(starFile, visPhased)
        visPhased *= np.exp(-1j * refPhase)
        
        syntheticVISDATA = createSyntheticCompanion(
            wav, visPhased, uCoord, vCoord, visOnStar,
            contrast, deltaRA_rad + fiberRA, deltaDec_rad + fiberDec
        )
        
        lastOutPath = saveSyntheticReducedData(planetFile, starFile, syntheticVISDATA, i, odir)
        
        if lastOutPath:
            os.chdir(lastOutPath)
            os.system('python /Users/svach/gravi_tools3/run_gravi_astrored_check.py')
            os.system('python /Users/svach/gravi_tools3/run_gravi_astrored_astrometry.py --reDo=TRUE')
            
def run(filePaths, resetPoint=0, injectionPath = '/Users/svach/gravityInjectionRecovery/injection_parameters.csv'):

    injectionParams = pd.read_csv(injectionPath)[resetPoint:]
    injectionParams['deltaRA_rad'] = injectionParams['deltaRA_mas'] * mas2rad
    injectionParams['deltaDec_rad'] = injectionParams['deltaDec_mas'] * mas2rad
    for filePath in filePaths:
        odir = filePath + 'synthetic_reduced/'
        os.makedirs(odir, exist_ok=True)
        reducedFiles = filePath + '*_astroreduced.fits'
        
        allFiles, starFiles, planetFiles = getFiles(reducedFiles)
        files = assignStar_to_Planets(allFiles, starFiles, planetFiles)
        
        # Prepare parameter sets for multiprocessing
        args = [
            (folderName, files, injectionParams, odir)
            for folderName in injectionParams.folderName.values[17:]
        ]

        with Pool(cpu_count()) as pool:
            pool.map(injectionRecovery, args)

#def performInjectionRecovery(files, injectionParams, odir):
#    deltaRAs_mas = injectionParams.deltaRA_mas.values
#    deltaDecs_mas = injectionParams.deltaDec_mas.values
#    deltaRAs_rad = deltaRAs_mas * mas2rad
#    deltaDecs_rad = deltaDecs_mas * mas2rad
#    contrasts = injectionParams.Contrast.values
#    folderNames = injectionParams.folderName.values
#    for folderName in folderNames:
##        print('***************** Injected Signal *****************')
##        print(r'$\Delta RA (rad)$ = ', deltaRAs_rad[i], r'$\Delta RA (mas)$ = ', deltaRAs_mas[i])
##        print(r'$\Delta Dec (rad)$ = ', deltaDecs_rad[i], r'$\Delta Dec (mas)$ = ', deltaDecs_mas[i])
##        print(r'Contrast = ', contrasts[i])
#
#        for planetFile, starFile in files:
#            wav, visPhased, uCoord, vCoord, fiberRA, fiberDec = openReducedPlanet(planetFile)
#            fiberRA = fiberRA * mas2rad
#            fiberDec = fiberDec * mas2rad
#            
#            visOnStar, refPhase = openReducedStar(starFile, visPhased)
#            visPhased = visPhased * np.exp(-1j * refPhase)
#
#            i = int(folderName)
#
#            syntheticVISDATA = createSyntheticCompanion(wav,
#                                     visPhased, 
#                                     uCoord,
#                                     vCoord,
#                                     visOnStar,
#                                     contrasts[i], 
#                                     deltaRAs_rad[i] + fiberRA,
#                                     deltaDecs_rad[i] + fiberDec)
#            outPath = saveSyntheticReducedData(planetFile, starFile, syntheticVISDATA, i, odir)
#
#        os.chdir(outPath)
#        os.system('python /Users/svach/gravi_tools3/run_gravi_astrored_check.py')
#        os.system('python /Users/svach/gravi_tools3/run_gravi_astrored_astrometry.py --reDo=TRUE')
    # plt.plot(0,0, 'r*', markersize=20)
    # plt.plot(fibreRA, fibreDec, 'k.')
    # plt.plot(fibreRA + deltaRAs_rad, fibreDec + deltaDecs_rad, 'b.')
    # circle1 = plt.Circle((fibreRA, fibreDec), 65, color='k')
    # plt.gca().set_aspect(1.0)
    # plt.xlabel(r'$Delta$RA')
    # plt.ylabel(r'$Delta$Dec')
    # plt.show()
        

            
        
        

    

if __name__ == '__main__':
    
    filePaths = [
                '/Users/svach/GravityGaia/p115/2025-08-10/SCI_HD 14082B_search1/reduced/',
#                        '/Users/svach/GravityGaia/p115/2025-08-09/SCI_HD 218396_search2/reduced/',
#                            '/Users/svach/GravityGaia/p115/2025-08-09/SCI_G 80-21_search_cc/reduced/']
 	]
    run(filePaths)
#    injectionParams = pd.read_csv('/Users/svach/gravityInjectionRecovery/injection_parameters.csv')
#    for filePath in filePaths:
#        #create new directory for synthetic data if doesn't exist
#        odir = filePath + 'synthetic_reduced/'
#        os.makedirs(odir, exist_ok=True)
#        reducedFiles = filePath + '*_astroreduced.fits'
#        allFiles, starFiles, planetFiles = getFiles(reducedFiles)
#        files = assignStar_to_Planets(allFiles, starFiles, planetFiles)
#        
#    #    separations_mas = np.linspace(0, 65/2, 5)   # fibre FWHM is 65 mas(?)
#    #    contrasts = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
#    #    injectionParams = make_injection_params(separations_mas, contrasts)
#
#        performInjectionRecovery(files, injectionParams[17:], odir)
        
        # open_reduced_planet(planetFiles[0])




