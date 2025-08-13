from astropy.io import fits
import glob
from datetime import datetime

import numpy as np
import pandas as pd


def getFiles(filePath):
    #from exoGRAVITY pipeline: https://gitlab.obspm.fr/mnowak/exogravity/-/blob/master/tutorial/EX1_data_manipulation.ipynb?ref_type=heads
    filenames = glob.glob(filePath)
    filenames.sort()
    
    starFiles, planetFiles = [], []
    for filename in filenames:
        hdul = fits.open(filename)
        hdr = hdul[0].header
        sObjX, sObjY = hdr["HIERARCH ESO INS SOBJ X"], hdr["HIERARCH ESO INS SOBJ Y"] 
        print("File {}:\n  Fiber offset: dRA = {} mas and dDEC = {} mas".format(filename, sObjX, sObjY))
        if (sObjX**2+sObjY**2) < 1:
            starFiles.append(filename)    
        else:
            planetFiles.append(filename)
        hdul.close()
    
    # print("On star files:\n", "\n ".join(starFiles))
    # print("On planet files:\n", "\n ".join(planetFiles))
    return filenames, starFiles, planetFiles

def assignStar_to_Planets(allFiles, starFiles, planetFiles):
    starSet = set(starFiles)
    planetSet = set(planetFiles)

    assigned = []
    last_star = None
    for f in allFiles:
        if f in starSet:
            last_star = f
        elif f in planetSet:
            if last_star is not None:
                assigned.append((f, last_star))
            else:
                assigned.append((f, None))
    return assigned


def getWavelengthGrid(filename):
    # extract the wavelength grid
    oi_wav = fits.getdata(filename, extname = "OI_WAVELENGTH", extver = 10)
    wav = oi_wav.field("EFF_WAVE")
    nwav = len(wav)
    return wav, nwav

def getVisPhased(oi_vis, nbase, nwav, wav):
    #from exoGRAVITY pipeline: 
    #https://gitlab.obspm.fr/mnowak/exogravity/-/blob/master/tutorial/EX1_data_manipulation.ipynb?ref_type=heads
    
    # also in Section 10.26.1, page 84 of the pipeline manual
    #https://ftp.eso.org/pub/dfs/pipelines/instruments/gravity/gravity-pipeline-manual-1.9.0.pdf
    visData = oi_vis.field("VISDATA").reshape([-1, nbase, nwav])
    opdDisp = oi_vis.field("OPD_DISP")
    phaseTelFc = oi_vis.field("PHASE_MET_TELFC")
    phaseRef = oi_vis.field("PHASE_REF")
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
    
    uCoord = oi_vis.field("UCOORD")
    vCoord = oi_vis.field("VCOORD")
    uvCoord = np.column_stack((uCoord, vCoord))
    
    return wav, visPhased, uvCoord

def openReducedStar(starFile, visPlanet):
    # extract the visdata and reshape it to a 3D array
    nbase = 6 # the number of baselines in GRAVITY
    oi_vis = fits.getdata(starFile, extname = "OI_VIS", extver = 10)
    wav, nwav = getWavelengthGrid(starFile)
    visPhased = getVisPhased(oi_vis, nbase, nwav, wav)

    #correct onstar shape
    visOnStar = np.mean(visPhased, axis=0)
    visOnStar = np.repeat(visOnStar[None, :, :], visPlanet.shape[0], axis=0)
    return visOnStar

def calculateOPD(u):
    
def createSyntheticSignal(contrast, deltaRA, deltaDec, signal='planet')
    if signal == 'planet':
        


def performInjectionRecovery(files, injectionParams):
    for planetFile, starFile in files:
        wav, visPhased, uvCoord = openReducedPlanet(planetFile)
        visOnStar = openReducedStar(starFile, visPhased)
        
        

    

if __name__ == '__main__':

    filePath = '/Users/svach/GravityGaia/p115/2025-08-10/hd14082B/reduced/*_astroreduced.fits'
    allFiles, starFiles, planetFiles = get_star_planet_files(filePath)
    files = assignStar_to_Planets(allFiles, starFiles, planetFiles)

    injectionParams = pd.read_csv('injection_parameters.csv')
    
    performInjectionRecovery(files, injectionParams)
    # open_reduced_planet(planetFiles[0])



