#
#  getGravityObserables.py
#
#
#  Created by Sydney Vach on 20.08.25.
#

from astropy.io import fits
import glob, os
from datetime import datetime
from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count

# Code Mathias
def calculatePhaseRefArclength(phaseref_coeffs, wav):
    """ 
    Calculate the arclength for the polynomial fit to phaseRef 
    args:
      phaseref_coeffs: array of size (ndit, nchannel, 4) containing the coeffs of the poly fit to phase ref
      wav: wavelength axis (in m)
    """
    # calculate normalized wavelength grid used to fit phaseRef
    w0 = wav.mean()
    wref = (w0/wav - 1)/(np.max(wav) - np.min(wav))*w0
    # prepare array for values
    ndit = np.shape(phaseref_coeffs)[0]
    nchannel = np.shape(phaseref_coeffs)[1]
    phaseRefArclength = np.zeros([ndit, nchannel])
    # calculate arclength
    for dit in range(ndit):
        for c in range(nchannel):
            coeffs = phaseref_coeffs[dit, c, :]
            # calculate the coeffs of 1+(dy/dx)**2
            arclengthpolycoeff = np.array([1+coeffs[1]**2, 4*coeffs[1]*coeffs[2], 4*coeffs[2]**2+6*coeffs[1]*coeffs[3], 12*coeffs[2]*coeffs[3], 9*coeffs[3]**2])
            # integrate sqrt(1+(dx/dy)**2) do get arclength. Minus sign because wref is decreasing
            phaseRefArclength[dit, c] = -np.trapz(np.sqrt(np.polyval(arclengthpolycoeff[::-1], wref)), wref)
    return phaseRefArclength

def getFiles(filePath):
    #from exoGRAVITY pipeline: https://gitlab.obspm.fr/mnowak/exogravity/-/blob/master/tutorial/EX1_data_manipulation.ipynb?ref_type=heads
    fileNames = sorted(glob.glob(filePath))
    
    starFiles, planetFiles = [], []
    for fileName in fileNames:
        with fits.open(fileName) as hdul:
            hdr = hdul[0].header
            sObjX, sObjY = hdr["HIERARCH ESO INS SOBJ X"], hdr["HIERARCH ESO INS SOBJ Y"]
            if (sObjX**2+sObjY**2) < 1:
                starFiles.append(fileName)
            else:
                planetFiles.append(fileName)

    return starFiles, planetFiles

def getWavelengthGrid(filename):
    # extract the wavelength grid
    oi_wav = fits.getdata(filename, extname = "OI_WAVELENGTH", extver = 10)
    wav = oi_wav.field("EFF_WAVE")

    return wav.copy(), len(wav)
    

    
    
def oiVis(fileName, extension, fast = True):
#    FT_Start=self.ftstart
    
    oi_vis_data={
        "TIME":[],
        'VISDATA_FT':[],
        'VISDATA':[],
        'VISERR':[],
        'PHASE_REF':[],
        'PHASE_REF_COEFF':[],
        'FLAG':[],
        'REJECTION_FLAG':[],
        'UCOORD':[],
        'VCOORD':[],
        'MJD':[],
        'OPD_DISP':[],
        'PHASE_MET_TELFC':[],
        'SELF_REF':[],
        'SELF_REF_COEFF':[],
        }
        
    hdul=fits.open(fileName,memmap=True,mode="readonly")
    data=hdul['OI_VIS',10].data

    ndit = hdul[0].header['HIERARCH ESO DET2 NDIT']
    
    for k in oi_vis_data.keys():
        data_column=data.field(k).reshape((ndit,6,-1))
    #    if k == "TIME":
    #        oi_vis_data[k]+=[data_column[:,0,0]*1e-6+FT_Start]
        if data_column.shape[2] == 1:
            oi_vis_data[k]=data_column[:,:,0]
        else:
            oi_vis_data[k]=data_column
    hdul.close()
    
    try:
        oi_vis_data["SEPTRK"] = hdul[0].header["HIERARCH ESO FT KAL SEPTRK"]
    except:
        oi_vis_data["SEPTRK"]= 0
    
    oi_vis_data['NDIT'] = ndit
    oi_vis_data['DIT'] = hdul[0].header['HIERARCH ESO DET2 SEQ1 DIT']
    oi_vis_data["SXY"]=np.array([hdul[0].header["HIERARCH ESO INS SOBJ X"],hdul[0].header["HIERARCH ESO INS SOBJ Y"]])
    try:
        oi_vis_data["SXY_CTU"] = np.array([(hdul[0].header["HIERARCH ESO FT KAL SEPCTUX%i"%i],hdul[0].header["HIERARCH ESO FT KAL SEPCTUY%i"%i]) for i in 1+np.arange(4)])
    except:
        oi_vis_data["SXY_CTU"] = np.array([oi_vis_data["SXY"][-1] for i in 1+arange(4)])
        
    
    if oi_vis_data["SEPTRK"] > 0.5:
        oi_vis_data["sxy"] = oi_vis_data["SXY"] + 0.
    else:
        oi_vis_data["sxy"] = oi_vis_data["SXY_CTU"].mean(axis=0)
    
    oi_vis_data['wave']=fits.getdata(fileName,'OI_WAVELENGTH',extension).field('EFF_WAVE')
    oi_vis_data['phaseMet']= 2*np.pi / oi_vis_data['wave'] * oi_vis_data['OPD_DISP'] + oi_vis_data['PHASE_MET_TELFC']
    oi_vis_data['ucoord']= oi_vis_data["UCOORD"][:,:,None]*2*np.pi/oi_vis_data['wave']/1000/(180/np.pi*3600)
    oi_vis_data['vcoord']= oi_vis_data["VCOORD"][:,:,None]*2*np.pi/oi_vis_data['wave']/1000/(180/np.pi*3600)
    oi_vis_data['uvmax'] = np.max(np.sqrt( oi_vis_data['ucoord']**2 + oi_vis_data['vcoord']**2 ))
    oi_vis_data['ucoord_diff'] = np.diff(oi_vis_data['ucoord'])
    oi_vis_data['vcoord_diff'] = np.diff(oi_vis_data['vcoord'])
    
    oi_vis_data['arcLength']=calculatePhaseRefArclength(oi_vis_data["PHASE_REF_COEFF"], oi_vis_data['wave'])
    
    oi_vis_data['var'] = np.abs(oi_vis_data["VISERR"]/oi_vis_data["DIT"])**2/2
    oi_vis_data['weight'] = 1./ oi_vis_data['var']
#    oi_vis_data['flag']=[]
    
#    if fast == True: # this might need to be readressed. work from astrored_loadData.py line 427-470
#        fact = ndit
#        ndit_compact= ndit / fact
#        phase=oi_vis_data['ucoord']*oi_vis_data['sxy'][0] + oi_vis_data['vcoord']*oi_vis_data['sxy'][1]
#        
#        def vel_col(vec, fact=fact, mean=False):
#            if mean:
#                return vec.reshape(fact, -1, *vec.shape[1:]).mean(axis=0)
        
    
    return oi_vis_data
    
def getVisPhased(data):
    vis = data['VISDATA'] / data['DIT']
    vis *= np.exp(1j*(data["PHASE_REF"]-data["phaseMet"]))
    return vis
    
def getVisUnphased(simVis, data):
    vis = simVis * data['DIT']
    vis *= np.exp(-1j * (data["PHASE_REF"]-data["phaseMet"]))
    return vis
    
def getVonStar(starFiles, extension = 10):

    visOnStar = []
    for starFile in starFiles:
        data = oiVis(starFile, extension, fast = True)
        visOnStar.append(getVisPhased(data))

    visRef = np.concatenate(visOnStar).mean(axis = 0)
#    phaseFt = np.concatenate(phaseFt).mean(axis = 0)
    phase_zp = np.angle(visRef)
    return visRef, phase_zp#, phaseFt
    
    
    
