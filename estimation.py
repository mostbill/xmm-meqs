'''
author: snowball@USTC
date: 2020.3.20
last update: 2020.3.20
fuction: derive SN ratios and determine if a change is significant
'''
import math
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np; np.set_printoptions(suppress=True)

from scipy import interpolate

from tqdm import tqdm

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

from luminosity import groupSource
from match import multiepoch

def readRedshifts():
    # read redshifts for every source
    global redshifts
    redshifts = np.loadtxt('redshifts', delimiter = ' ')

def readBrevityLums():
    # read lums for brevity
    global brevity_uv_lums; global brevity_xr_lums
    brevity_uv_lums = np.loadtxt('uv_brevity_lums', delimiter = ' ')
    brevity_xr_lums = np.loadtxt('xr_brevity_lums', delimiter = ' ')

def readBrevityErrors():
    # read errors for brevity
    global brevity_uv_errors; global brevity_xr_errors
    brevity_uv_errors = np.loadtxt('uv_brevity_errors', delimiter = ' ')
    brevity_xr_errors = np.loadtxt('xr_brevity_errors', delimiter = ' ')

def singleLumsToFluxes(source):
    # lums to fluxes
    index = source
    _brevity_uv_lum = brevity_uv_lums[index]
    _brevity_xr_lum = brevity_xr_lums[index]
    _redshift = redshifts[index]

    D_Lum = cosmo.luminosity_distance(_redshift).to(u.m)
    _brevity_uv_flux = ((1.0 + _redshift) * 10.0**(_brevity_uv_lum)) / (4.0 * math.pi * D_Lum.value**2.0)
    _brevity_uv_flux = _brevity_uv_flux * (1.0 / (3.34 * 10.0**-15.0 * 2500.0**2.0)) # flux unit conversion
    gamma = 1.7
    _brevity_xr_flux = ((1.0 + _redshift)**(2-gamma) * 10.0**(_brevity_xr_lum)) / (4.0 * math.pi * D_Lum.value**2.0)
    _brevity_xr_flux = _brevity_xr_flux * ((4.84 * 10.0**17.0) * (2.25**(2.0-gamma) - 1.0)) / (10.0**4.0 * (2.0-gamma)) # flux unit conversion

    return _brevity_uv_flux, _brevity_xr_flux

def lumsToFluxes(sources):
    # compute brevity fluxes for all sources
    global brevity_uv_fluxes; brevity_uv_fluxes = []
    global brevity_xr_fluxes; brevity_xr_fluxes = []
    for i in tqdm(range(len(sources))):
        # compute brevity flux for every source (single/multi epochs)
        if(len(sources[i]) == 1):
            _brevity_uv_flux, _brevity_xr_flux = singleLumsToFluxes(sources[i])
            brevity_uv_fluxes.append(_brevity_uv_flux)
            brevity_xr_fluxes.append(_brevity_xr_flux)
        elif(len(sources[i]) > 1):
            for j in range(len(sources[i])):
                _brevity_uv_flux, _brevity_xr_flux = singleLumsToFluxes(sources[i][j])
                brevity_uv_fluxes.append(_brevity_uv_flux)
                brevity_xr_fluxes.append(_brevity_xr_flux)
    
    np.savetxt("uv_brevity_fluxes", brevity_uv_fluxes, fmt = '%E', delimiter = ' ')
    np.savetxt("xr_brevity_fluxes", brevity_xr_fluxes, fmt = '%E', delimiter = ' ')    

def computeSNRatio(sources):
    # compute SN ratio (for a change) for all sources
    # SN ratio := abs(s1-s2)/sqrt(n1^2+n2^2)
    uv_snratios = []; xr_snratios = []
    sources_locats = [] # save the former index of every SNratios' entries (uv/xr)
    for i in tqdm(range(len(sources))):
        # compute SNratio for every (multi epochs) source
        if(len(sources[i]) > 1):
            for j in range(len(sources[i])-1):
                index = sources[i][j]; index_next = sources[i][j+1]
                uv_delta_f = brevity_uv_fluxes[index_next]-brevity_uv_fluxes[index]
                uv_error_square = np.sqrt(brevity_uv_errors[index_next]**2.0 + brevity_uv_errors[index]**2.0)
                xr_delta_f = brevity_xr_fluxes[index_next]-brevity_xr_fluxes[index]
                xr_error_square = np.sqrt(brevity_xr_errors[index_next]**2.0 + brevity_xr_errors[index]**2.0)

                uv_snratios.append(np.abs(uv_delta_f)/uv_error_square)
                xr_snratios.append(np.abs(xr_delta_f)/xr_error_square) 
                sources_locats.append(index)            
    
    np.savetxt("uv_snratios", uv_snratios, fmt = '%f', delimiter = ' ')
    np.savetxt("xr_snratios", xr_snratios, fmt = '%f', delimiter = ' ')
    np.savetxt("snratios_sources_locats", sources_locats, fmt = '%d', delimiter = ' ')

def computeSync(sources):
    # compute Sync-index for every source
    # Sync-index := (X-axis: Delta_X-ray_lum, Y-axis: Delta_UV_lum)
    delta_uv_lums = []; delta_xr_lums = [];
    sources_locats = [] # save the former index of every deltas' entries (uv/xr)
    for i in tqdm(range(len(sources))):
        # compute SNratio for every (multi epochs) source
        if(len(sources[i]) > 1):
            for j in range(len(sources[i])-1):
                index = sources[i][j]; index_next = sources[i][j+1]
                _delta_uv_lums = brevity_uv_lums[index_next]-brevity_uv_lums[index]
                _delta_xr_lums = brevity_xr_lums[index_next]-brevity_xr_lums[index]

                delta_uv_lums.append(_delta_uv_lums)
                delta_xr_lums.append(_delta_xr_lums) 
                sources_locats.append(index)        
    
    np.savetxt("delta_uv_lums", delta_uv_lums, fmt = '%f', delimiter = ' ')
    np.savetxt("delta_xr_lums", delta_xr_lums, fmt = '%f', delimiter = ' ')
    np.savetxt("delta_sources_locats", sources_locats, fmt = '%d', delimiter = ' ') 

def readSNratios():
    # read the locats for sources 
    global uv_snratios; global xr_snratios; global source_locats
    uv_snratios = np.loadtxt('uv_snratios', delimiter=' ')
    xr_snratios = np.loadtxt('xr_snratios', delimiter=' ')
    source_locats = np.loadtxt('snratios_sources_locats', delimiter=' ')

def computeR2(sources):
    # compute Pearson R2
    xr_sn_threshold = 6.327
    uv_sn_threshold = 20

    all_coef = []
    for i in tqdm(range(len(sources))):
        if(len(sources[i]) > 1):
            # only for multi-epoch sources
            brevity_uv_lums_vector = []
            brevity_xr_lums_vector = []
            for j in range(len(sources[i])):
                try:
                    index = np.where(source_locats == sources[i][j])
                    if(float(uv_snratios[index]) > uv_sn_threshold and float(xr_snratios[index]) > xr_sn_threshold):
                        for k in range(len(sources[i])):
                            brevity_uv_lums_vector.append(brevity_uv_lums[sources[i][k]])
                            brevity_xr_lums_vector.append(brevity_xr_lums[sources[i][k]])
                        break
                except:
                    pass
            if(len(brevity_uv_lums_vector)>0):
                all_coef.append(np.corrcoef([brevity_uv_lums_vector, brevity_xr_lums_vector])[0][1])

    np.savetxt("20_sn_coef", all_coef, fmt = '%f', delimiter = ' ')

def main():
    workpath = '/Users/snowball/astro/workspace/XMM-OM/'
    lists = fits.open(workpath + 'new.fits')
    flags = multiepoch(lists) # indicate sources which belong to a single source

    sources = groupSource(flags)

    readBrevityLums(); readBrevityErrors(); readRedshifts()
    # lumsToFluxes(sources)
    # computeSNRatio(sources)
    # computeSync(sources)
    readSNratios(); computeR2(sources)

if __name__ == "__main__":
    main()