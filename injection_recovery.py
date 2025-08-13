from astropy.io import fits
import glob, os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mas_to_rad = np.pi / (180.0 * 3600.0 * 1000.0)

def getFiles(filePath):
    #from exoGRAVITY pipeline: https://gitlab.obspm.fr/mnowak/exogravity/-/blob/master/tutorial/EX1_data_manipulation.ipynb?ref_type=heads
    filenames = glob.glob(filePath)
    filenames.sort()
    
    starFiles, planetFiles = [], []
    for filename in filenames:
        hdul = fits.open(filename)
        hdr = hdul[0].header
        sObjX, sObjY = hdr["HIERARCH ESO INS SOBJ X"], hdr["HIERARCH ESO INS SOBJ Y"] 
        # print("File {}:\n  Fiber offset: dRA = {} mas and dDEC = {} mas".format(filename, sObjX, sObjY))
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
    
    uCoord = oi_vis.field("UCOORD").reshape([-1, nbase])
    vCoord = oi_vis.field("VCOORD").reshape([-1, nbase])

    hdul = fits.open(planetFile)
    hdr = hdul[0].header
    sObjX, sObjY = hdr["HIERARCH ESO INS SOBJ X"], hdr["HIERARCH ESO INS SOBJ Y"]
    hdul.close()
    return wav, visPhased, uCoord, vCoord, sObjX, sObjY

def openReducedStar(starFile, visPlanet):
    nbase = 6 # the number of baselines
    oi_vis = fits.getdata(starFile, extname = "OI_VIS", extver = 10)
    wav, nwav = getWavelengthGrid(starFile)
    visPhased = getVisPhased(oi_vis, nbase, nwav, wav)

    #correct onstar shape
    visOnStar = np.mean(visPhased, axis=0)
    visOnStar = np.repeat(visOnStar[None, :, :], visPlanet.shape[0], axis=0)
    return visOnStar

def createSyntheticCompanion(wav, visPhased, uCoord, vCoord, visOnStar, contrast, deltaRA, deltaDec):
    ## synthetic Companioin signal
    ## as per Pourré, N., et al.: A&A, 686, A258 (2024):
    ## Vsyn,comp = C*Vonstar * exp[(-i2pi/lambda)(OPD)]exp[i phi] : need to come back to check on if this phase term is properly accounted for!!!!!

    opd = uCoord * deltaRA + vCoord * deltaDec
    # print(np.shape(opd))
    opd_2d = opd[:, None] 
    # print(np.shape(opd_2d))
    opd_3d = opd[:, :, None] 
    # print(np.shape(opd_3d))

    phase_term = np.exp(-2j * np.pi * opd_3d / wav[None, None, :])

    visSyntheticComp = visPointing + contrast * visOnStar * phaseShift
    return visSyntheticComp
        
def saveSyntheticReducedData(planetFile, visSyntheticComp, i):
    hdul = fits.open(planetFile)
    oi_vis = fits.getdata(planetFile, extname = "OI_VIS", extver = 10)
    visData = oi_vis.field("VISDATA")
    print(visData)
    for hdu in hdul:
        if hdu.header.get('EXTNAME') == 'OI_VIS':
#            oi_vis_hdu = hdu
            ndit, nbase, nwave = np.shape(visSyntheticComp)
            visSyntheticReshaped = visSyntheticComp.reshape(ndit * nbase, nwave)
            hdu.data['VISDATA'] = visSyntheticReshaped
            break
    print(hdul[5].data['VISDATA'])
    

def performInjectionRecovery(files, injectionParams):
    for planetFile, starFile in files:
        wav, visPhased, uvCoord, fibreRA, fibreDec = openReducedPlanet(planetFile)
        visOnStar = openReducedStar(starFile, visPhased)

        contrasts = injectionParams.Contrast.values
        deltaRAs_mas = injectionParams.deltaRA_mas.values
        deltaDecs_mas = injectionParams.deltaDec_mas.values
        
        deltaRAs_rad = deltaRAs_mas * mas_to_rad
        deltaDecs_rad = deltaDecs_mas * mas_to_rad
        for i in range(len(contrasts)):
            syntheticVISDATA = createSyntheticCompanion(wav, 
                                     visPhased, 
                                     uCoord,
                                     vCoord,
                                     visOnStar,
                                     contrasts[i], 
                                     deltaRAs_rad[i], 
                                     deltaDecs_rad[i])
                                     
            saveSyntheticReducedData(planetFile, syntheticVISDATA, i)
            

    
        # plt.plot(0,0, 'r*', markersize=20)
        # plt.plot(fibreRA, fibreDec, 'k.')
        # plt.plot(fibreRA + deltaRAs_rad, fibreDec + deltaDecs_rad, 'b.')
        # circle1 = plt.Circle((fibreRA, fibreDec), 65, color='k')
        # plt.gca().set_aspect(1.0)
        # plt.xlabel(r'$Delta$RA')
        # plt.ylabel(r'$Delta$Dec')
        # plt.show()
        

            
        
        

    

if __name__ == '__main__':

    filePath = '/Users/svach/GravityGaia/p115/2025-08-10/hd14082B/reduced/'
    #create new directory for synthetic data if doesn't exist
    os.makedirs(filePath + 'synthetic_reduced', exist_ok=True) 
    
    reducedFiles = filePath + '_astroreduced.fits'
    allFiles, starFiles, planetFiles = getFiles(reducedFiles)
    files = assignStar_to_Planets(allFiles, starFiles, planetFiles)

    injectionParams = pd.read_csv('injection_parameters.csv')
    
    performInjectionRecovery(files, injectionParams)
    # open_reduced_planet(planetFiles[0])



