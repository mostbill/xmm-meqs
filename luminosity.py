'''
author: snowball@USTC
date: 2020.3.2
last update: 2020.8.20
fuction: derive UV and X-ray luminosities (for brevity) within catalogues
'''
import math
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

from scipy import interpolate

from tqdm import tqdm

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

from match import multiepoch

# colors from JT
black = (0, 0, 0)
orange = (230/255., 159/255., 0)
blue = (86/255., 180/255., 233/255.)
green = (0, 158/255., 115/255.)
yellow = (240/255., 228/255., 66/255.)
dblue = (0, 114/255., 178/255.)
vermillion = (213/255., 94/255., 0)
purple = (204/255., 121/255., 167/255.)

uv_bands = {'UVW2': 1894, 'UVM2': 2205, 'UVW1': 2675, 'U': 3275, 'B': 4050, 'V': 5235}

def groupSource(flags): 
    # group entries to a single source
    sources = []
    source = []
    for i in tqdm(range(len(flags))):
        source.append(i) # indicator of entries in the full list
        if(i+1 < len(flags)):
            if(flags[i+1] != flags[i]):
                sources.append(source)
                source = []
        else:
            sources.append(source)
    return(sources)

def readUVFlux(lists):
    # read luminosities from fits
    global u_fluxes; global v_fluxes; global b_fluxes; global uvw1_fluxes; global uvm2_fluxes; global uvw2_fluxes
    u_fluxes = lists[1].data['U_AB_FLUX'][:]; v_fluxes = lists[1].data['V_AB_FLUX'][:]; b_fluxes = lists[1].data['B_AB_FLUX'][:]
    uvw1_fluxes = lists[1].data['UVW1_AB_FLUX'][:]; uvm2_fluxes = lists[1].data['UVM2_AB_FLUX'][:]; uvw2_fluxes = lists[1].data['UVW2_AB_FLUX'][:]
 
def readUVError(lists):
    # read luminosities' errors from fits
    global u_errors; global v_errors; global b_errors; global uvw1_errors; global uvm2_errors; global uvw2_errors
    u_errors = lists[1].data['U_AB_FLUX_ERR'][:]; v_errors = lists[1].data['V_AB_FLUX_ERR'][:]; b_errors = lists[1].data['B_AB_FLUX_ERR'][:]
    uvw1_errors = lists[1].data['UVW1_AB_FLUX_ERR'][:]; uvm2_errors = lists[1].data['UVM2_AB_FLUX_ERR'][:]; uvw2_errors = lists[1].data['UVW2_AB_FLUX_ERR'][:]
     
def readXrFlux(lists):
    # read luminosities from fits
    global ep_4_fluxes
    ep_4_fluxes = lists[1].data['EP_4_FLUX'][:]

def readXrError(lists):
    # read luminosities' errors from fits
    global ep_4_errors
    ep_4_errors = lists[1].data['EP_4_FLUX_ERR'][:]

def readRedshift(lists):
    # read redshifts from fits
    global redshifts
    redshifts = lists[1].data['Z'][:]
    redshifts = np.array(redshifts)

    np.savetxt('redshifts', redshifts, fmt = '%f', delimiter= ' ')

def readRaDec(lists):
    # read ra, dec from fits
    global ra, dec
    ra = lists[1].data['RA_1'][:]; dec = lists[1].data['DEC_1'][:]
    ra = np.array(ra); dec = np.array(dec)

    np.savetxt('ra', ra, fmt = '%f', delimiter= ' ')
    np.savetxt('dec', dec, fmt = '%f', delimiter= ' ')

def readID(lists):
    # read srcid, detid from fits
    global srcid, detid
    srcid = lists[1].data['SRCID'][:]; detid = lists[1].data['DETID'][:]
    srcid = np.array(srcid); detid = np.array(detid)

    np.savetxt('srcid', srcid, fmt = '%f', delimiter= ' ')
    np.savetxt('detid', detid, fmt = '%f', delimiter= ' ')   

def readBrevityLums():
    # read lums for brevity
    global brevity_uv_lums; global brevity_xr_lums
    brevity_uv_lums = np.loadtxt(workpath + 'uv_brevity_lums', delimiter=' ')
    brevity_xr_lums = np.loadtxt(workpath + 'xr_brevity_lums', delimiter = ' ')

def readMeanSED():
    # read mean SEDs of QSOs from Richards et al. 2006
    workpath = '/Users/snowball/astro/workspace/XMM-OM/'
    lists = fits.open(workpath + 'mean_SED.fits')
    global Logfreqs; global Logalllums
    Logfreqs = lists[1].data['LogFreq'][:]; Logalllums = lists[1].data['LogLall']
    Logalllums = Logalllums - Logfreqs # vL -> L

    plt.xlabel(r'log $\nu_{emit}$ (Hz)', fontsize = 12)
    plt.ylabel(r'log $L_{UV}(ergs\ s^{-1}\ Hz^{-1})$', fontsize = 12)
    plt.plot(Logfreqs, Logalllums, linewidth=2, color = 'black')
    # plt.ylim(28, 32)

def singleAngstromToHz(Angstroms, redshift):
    # for a single entry, invert angstrom to Hz
    c = 3*10**8 # velocity of light
    Angstroms = np.array(Angstroms)
    Hzs = (c/(Angstroms*10**-10)) * (1 + redshift)
    Hzs = np.log10(Hzs)
    return Hzs

def singleFluxToMag(fluxes, Angstroms, redshift):
    # use cosmological distance, compute luminosity from flux
    D_Lum = cosmo.luminosity_distance(redshift).to(u.m)
    Angstroms = np.array(Angstroms)
    fluxes = np.array(fluxes).reshape(6,)
    fluxes = fluxes * (3.34 * 10**-15) * Angstroms**2 # flux unit conversion
    lums = fluxes * (4.0 * math.pi * D_Lum.value**2)/(1+redshift)
    lums = np.log10(lums)
    return lums

def readSourceLums(source):
    # read source luminsoities and frequencies
    index = source
    uv_bands = [1894, 2205, 2675, 3275, 4050, 5235]
    frequencies = singleAngstromToHz(uv_bands, redshifts[index])
    fluxes = []
    fluxes.append(uvw2_fluxes[index]); fluxes.append(uvm2_fluxes[index]); fluxes.append(uvw1_fluxes[index])
    fluxes.append(u_fluxes[index]); fluxes.append(b_fluxes[index]); fluxes.append(v_fluxes[index])
    lums = singleFluxToMag(fluxes, uv_bands, redshifts[index])

    return frequencies, lums

def readSrouceErrors(source):
    # (in flux's unit) read source flux errors
    index = source
    uv_flux_errors = []; xr_flux_errors = []
    uv_flux_errors.append(uvw2_errors[index]); uv_flux_errors.append(uvm2_errors[index]); uv_flux_errors.append(uvw1_errors[index])
    uv_flux_errors.append(u_errors[index]); uv_flux_errors.append(b_errors[index]); uv_flux_errors.append(v_errors[index])
    xr_flux_errors.append(ep_4_errors[index]) 

    return uv_flux_errors, xr_flux_errors

def plotSet():
    plt.clf()
    plt.figure(figsize=(6, 6))          # Make a new figure window
    plt.xlabel(r'$log\ \nu_{emit}\ (Hz)$', fontsize = 12)
    plt.ylabel(r'$log\ L_{UV}\ (erg\ s^{-1}\ Hz^{-1})$', fontsize = 12)
    band_brevity = [2500]
    freq_brevity = singleAngstromToHz(band_brevity, 0)
    plt.axvline(freq_brevity, linestyle = '--')

def plotSetAlpha():
    plt.clf()
    plt.figure(figsize=(6, 6))          # Make a new figure window
    # plt.ylabel(r'$\alpha_{ox}$', fontsize = 12)
    plt.ylabel(r'$\alpha_{OX}$', fontsize = 12)
    # plt.ylabel(r'$log\ L_{X}\ (erg\ s^{-1}\ Hz^{-1})$', fontsize = 12)
    plt.xlabel(r'$log\ L_{UV}\ (erg\ s^{-1}\ Hz^{-1})$', fontsize = 12)

def redshiftDist(sources):
    # plot redshifts span of the dataset
    plt.figure(figsize=(6, 6))
    plt.xlabel(r'$log\ L_{UV}\ (erg\ s^{-1}\ Hz^{-1})$', fontsize = 12)
    plt.ylabel(r'$log\ z$', fontsize = 12)

    global z_lums; z_lums = []
    global z_log10_zs; z_log10_zs = []

    for i in tqdm(range(len(sources))):
        if(len(sources[i]) == 1):
            singlePlotZ(sources[i])
        elif(len(sources) > 1):
            for j in range(len(sources[i])):
                singlePlotZ(sources[i][j])

    plt.scatter(z_lums, z_log10_zs, s = 0.5)
    plt.savefig('z_distribution.pdf', dpi=2000)

def singlePlotZ(source):
    # plot redshift for source with only 1 epoch
    index = source
    z_lums.append(brevity_uv_lums[index])
    z_log10_zs.append(np.log10(redshifts[index]))
    
def plotSED(sources):
    # plot SED for the whole dataset
    for i in tqdm(range(len(sources))):
        if(len(sources[i]) == 1):
            sourcePlotSED(sources[i])
        elif(len(sources[i]) > 1):
            for j in range(len(sources[i])):
                sourcePlotSED(sources[i][j])
        else:
            pass
            # multiPlotSED
    plt.xlim(14.6, 15.8)
    xmin, xmax = plt.xlim()
    plt.xticks(np.arange(xmin, xmax, step=0.2))
    plt.savefig('all_SED.pdf', dpi=2000)

def calcAlpha(sources):
    # calculate alpha for all detections
    alphas = []

    for i in tqdm(range(len(sources))):
        if(len(sources[i]) == 1):
            index = sources[i]
            alphas.append(0.3838 * (brevity_xr_lums[index] - brevity_uv_lums[index]))
        elif(len(sources[i]) > 1):
            for j in range(len(sources[i])):
                index = sources[i][j]
                alphas.append(0.3838 * (brevity_xr_lums[index] - brevity_uv_lums[index]))

    np.savetxt("all_alphas", alphas, fmt = '%f', delimiter = ' ')

def plotSingleAlpha(source):
    # plot alpha for single-epoch source
    index = source
    try:
        alpha = 0.3838 * (brevity_xr_lums[index] - brevity_uv_lums[index])
        # plt.scatter(brevity_uv_lums[index], alpha, marker = '.', color = 'black', zorder = 1, alpha = 0.2)
        plt.scatter(brevity_uv_lums[index], brevity_xr_lums[index], marker = '.', color = 'black', zorder = 1, alpha = 0.2)
        # pass
    except:
        pass

def plotMultiAlpha(source):
    # plot alpha for multi-epochs source
    alphas = []; _uv_lums = []; _xr_lums = []
    # start_time = []
    try:
        for i in range(len(source)):
            alphas.append(0.3838 * (brevity_xr_lums[source[i]] - brevity_uv_lums[source[i]]))
            _uv_lums.append(brevity_uv_lums[source[i]])
            _xr_lums.append(brevity_xr_lums[source[i]])
            # start_time.append(mjd_start[source[i]])
        alphas = np.array(alphas); _uv_lums = np.array(_uv_lums); _xr_lums = np.array(_xr_lums)
        _uv_lums = _uv_lums[np.argsort(start_time)]; _xr_lums = _xr_lums[np.argsort(start_time)]; alphas = alphas[np.argsort(start_time)]
        plt.plot(_uv_lums, alphas, color = 'black', linewidth = 0.5, zorder = -1)
        plt.scatter(_uv_lums[len(_uv_lums)-1], alphas[len(alphas)-1], marker = '.', color = 'black', s=10)
        # plt.plot(_uv_lums, _xr_lums, color = 'black', linewidth = 0.5, zorder = -1)
    except:
        pass

def sourcePlotSED(source):
    # plot SED for one single source
    index = source
    uv_bands = [1894, 2205, 2675, 3275, 4050, 5235]
    frequencies = singleAngstromToHz(uv_bands, redshifts[index])
    fluxes = []
    fluxes.append(uvw2_fluxes[index]); fluxes.append(uvm2_fluxes[index]); fluxes.append(uvw1_fluxes[index])
    fluxes.append(u_fluxes[index]); fluxes.append(b_fluxes[index]); fluxes.append(v_fluxes[index])
    lums = singleFluxToMag(fluxes, uv_bands, redshifts[index])

    if(np.count_nonzero(~np.isnan(lums)) == 1):
        # with 1 frequency point
        plt.scatter(frequencies, lums, marker = '.', color = 'black', zorder = 1, alpha = 0.2)
    elif(np.count_nonzero(~np.isnan(lums)) > 1):
        # with >1 frequency point
        plt.plot(frequencies, lums, zorder = -1)

def interpolateSED(source, flag):
    # interpolate for brevity point in SED
    frequencies, lums = readSourceLums(source)
    uv_errors, _ = readSrouceErrors(source) # units in flux
    brevity_frequency = singleAngstromToHz(2500, 0)

    if flag == 1: # brevity point fall on the regular interval
        f = interpolate.interp1d(frequencies[1:3], lums[1:3])
        error_brevity = np.sqrt(uv_errors[1]**2 + uv_errors[2]**2)
    elif flag == 2: # brevity point fall on the bigger interval
        f = interpolate.interp1d(frequencies[2:4], lums[2:4])
        error_brevity = np.sqrt(uv_errors[2]**2 + uv_errors[3]**2)
    elif flag == 3:
        f = interpolate.interp1d(frequencies[3:5], lums[3:5])
        error_brevity = np.sqrt(uv_errors[3]**2 + uv_errors[4]**2)
    elif flag == 4:
        f = interpolate.interp1d(frequencies[4:], lums[4:])
        error_brevity = np.sqrt(uv_errors[4]**2 + uv_errors[5]**2)

    lum_brevity = f(brevity_frequency)
    return lum_brevity, error_brevity

def extrapolateSED(source):
    # extrapolate for brevity point in SED
    frequencies, lums = readSourceLums(source)
    uv_errors, _ = readSrouceErrors(source) # units in flux
    brevity_frequency = singleAngstromToHz(2500, 0)

    indexes = np.argsort(np.abs(frequencies-brevity_frequency)) # find closest point in source
    for i in range(len(indexes)):
        # if index = np.nan, find next band
        if(~np.isnan(lums[indexes[i]])):
            index = indexes[i]       
    error_brevity = uv_errors[index] # units in flux
    _index = np.nanargmin(np.abs(Logfreqs-frequencies[index])) # find closest (to point above) point in the mean SED
    _Logalllums = Logalllums; _Logalllums = np.array(_Logalllums)
    _Logalllums = _Logalllums - (_Logalllums[_index] - lums[index])
    
    f = interpolate.interp1d(Logfreqs, _Logalllums) # mean SED is not continuous, so use interpolate
    lum_brevity = f(brevity_frequency)
    return lum_brevity, error_brevity

def trueInterpolate(source):
    # check if interpolate-able or not
    index = source
    if(2500*(1+redshifts[index]) < 2675):
        if(~(np.isnan(uvm2_fluxes[index])) and ~(np.isnan(uvw1_fluxes[index]))):
            return 1
    elif(2500*(1+redshifts[index]) >= 2675 and 2500*(1+redshifts[index]) < 3275):
        if(~(np.isnan(uvw1_fluxes[index])) and ~(np.isnan(u_fluxes[index]))):
            return 2  
    elif(2500*(1+redshifts[index]) >= 3275 and 2500*(1+redshifts[index]) < 4050):
        if(~(np.isnan(u_fluxes[index])) and ~(np.isnan(b_fluxes[index]))):
            return 3 
    elif(2500*(1+redshifts[index]) >= 4050 and 2500*(1+redshifts[index]) < 5235):
        if(~(np.isnan(b_fluxes[index])) and ~(np.isnan(v_fluxes[index]))):
            return 4 

    return False

def directXrSED(source):
    # directly compute X-ray luminosity for brevity
    _, error_brevity = readSrouceErrors(source) # units in flux

    index  = source
    gamma = 1.7
    brevity_frequency = 4.84 * 10**17
    brevity_flux = ((ep_4_fluxes[index] * 10**4) / brevity_frequency) * ((2 - gamma) / (2.25**(2-gamma)-1))
    D_Lum = cosmo.luminosity_distance(redshifts[index]).to(u.m)
    lum_brevity = brevity_flux * (4 * math.pi * D_Lum.value**2)/(1+redshifts[index])**(2-gamma)
    lum_brevity = np.log10(lum_brevity)

    return lum_brevity, error_brevity

def computeUVBrevityFlux(sources):
    # compute UV flux for brevity for every source (with error) and save as file
    brevity_lums = []
    brevity_errors = []
    for i in tqdm(range(len(sources))):
        if(len(sources[i]) == 1):
            # only one-epoch point, surely go for extrapolate
            lum_brevity, error_brevity = extrapolateSED(sources[i])
            brevity_lums.append(lum_brevity); brevity_errors.append(error_brevity)
        elif(len(sources[i]) > 1):
            for j in range(len(sources[i])):
                if(trueInterpolate(sources[i][j])): # brevity point is interpolate-able
                    lum_brevity, error_brevity = interpolateSED(sources[i][j], trueInterpolate(sources[i][j]))
                    brevity_lums.append(lum_brevity); brevity_errors.append(error_brevity)
                else:
                    lum_brevity, error_brevity = extrapolateSED(sources[i][j])
                    brevity_lums.append(lum_brevity); brevity_errors.append(error_brevity)
        else:
            pass
            # multiPlotSED

    np.savetxt(workpath + "uv_brevity_lums", brevity_lums, fmt = '%f', delimiter = ' ')
    np.savetxt(workpath + "uv_brevity_errors", brevity_errors, fmt = '%E', delimiter = ' ')

def computeXrBrevityFlux(sources):
    # compute Xr flux for brevity for every source
    brevity_lums = []
    brevity_errors = []
    for i in tqdm(range(len(sources))):
        if(len(sources[i]) == 1):
            lum_brevity, error_brevity = directXrSED(sources[i])
            brevity_lums.append(lum_brevity)
            brevity_errors.append(error_brevity)
        elif(len(sources[i]) > 1):
            for j in range(len(sources[i])):
                lum_brevity, error_brevity = directXrSED(sources[i][j])
                brevity_lums.append(lum_brevity)
                brevity_errors.append(error_brevity)

    np.savetxt("xr_brevity_lums", brevity_lums, fmt = '%f', delimiter = ' ')
    np.savetxt("xr_brevity_errors", brevity_errors, fmt = '%E', delimiter = ' ')

def plotAlphaOX(sources):
    # compute and plot alpha_ox for every source
    for i in tqdm(range(len(sources))):
        if(len(sources[i]) == 1):
            plotSingleAlpha(sources[i])
        elif(len(sources[i]) > 1):
            plotMultiAlpha(sources[i])

def readSNratios():
    # read the locats for sources 
    global uv_snratios; global xr_snratios; global source_locats
    uv_snratios = np.loadtxt(workpath + 'uv_snratios', delimiter=' ')
    xr_snratios = np.loadtxt(workpath + 'xr_snratios', delimiter=' ')
    source_locats = np.loadtxt(workpath + 'snratios_sources_locats', delimiter=' ')
    
def plotSelectedAlphaOX(sources):
    # plot alpha_ox for selected (SN ratio value) source
    xr_sn_threshold = 0
    uv_sn_threshold = 0

    for i in tqdm(range(len(sources))):
        if(len(sources[i]) == 1):
            # plotSingleAlpha(sources[i]) # useless
            pass
        elif(len(sources[i]) > 1):
            for j in range(len(sources[i])):
                try:
                    index = np.where(source_locats == sources[i][j])
                    if(float(uv_snratios[index]) > uv_sn_threshold and float(xr_snratios[index]) > xr_sn_threshold):
                        plotMultiAlpha(sources[i])
                        index = -1
                except:
                    pass

def linearRegression(x, y): 
    # linear least squares fit 
    # y = a0 + a1*x
    x = np.array(x); y = np.array(y)
    N = len(x) 
    sumx = sum(x) 
    sumy = sum(y) 
    sumx2 = sum(x**2) 
    sumxy = sum(x*y) 

    # print(sumx); print(sumy)
    # print(N, sumx[0], sumx2[0], sumy[0], sumxy[0])

    try:
        A = np.mat([[N, sumx], [sumx, sumx2]]) 
        b = np.array([sumy, sumxy]) 
    except:
        A = np.mat([[N, float(sumx[0])], [float(sumx[0]), float(sumx2[0])]]) 
        b = np.array([float(sumy[0]), float(sumxy[0])])         

    # A = np.mat([[518, 15530.786], [15530.786, 465840.395]])
    # b = np.array([-782.632, -23494.554])
    a0, a1 = np.linalg.solve(A, b)

    return a0, a1

def doLinearFitting():
    alphas = np.loadtxt(workpath + 'all_alphas', delimiter=' ')
    _alphas = []; _uv_brevity_lums = [] # with no value of -inf

    for i in range(len(alphas)):
        if(alphas[i]>-999): # not -inf
            _alphas.append(alphas[i])
            _uv_brevity_lums.append(brevity_uv_lums[i])
    
    a0, a1 = linearRegression(_uv_brevity_lums, _alphas)
    print(a0, a1)
    bootstrapLF(_uv_brevity_lums, _alphas, a0, a1, 4000) # 4000 times of bootstrap, more trials than this will not bring any more knowledge

def doSelectedLinearFitting(sources):
    # do linear fitting to selected sources
    xr_sn_threshold = 6.3273
    uv_sn_threshold = 20

    alphas = np.loadtxt(workpath + 'all_alphas', delimiter=' ')
    _alphas = []; _uv_brevity_lums = []

    for i in tqdm(range(len(sources))):
        if(len(sources[i]) == 1):
            pass
        elif(len(sources[i]) > 1):
            for j in range(len(sources[i])):
                try:
                    index = np.where(source_locats == sources[i][j])
                    if(float(uv_snratios[index]) > uv_sn_threshold and float(xr_snratios[index]) > xr_sn_threshold):
                        if(alphas[index]>-999):
                            _alphas.append(alphas[index]); _uv_brevity_lums.append(brevity_uv_lums[index])
                            index = -1
                except:
                    pass   
    
    linearRegression(_uv_brevity_lums, _alphas)

def bootstrapLF(x, y, m0, m1, times):
    # carry out error estimation through bootstrap method

    temp_m0 = []; temp_m1 = []

    for i in range(times):
        index = np.random.choice(len(x), len(x))
        _x = np.array(x)[index]
        _y = np.array(y)[index]

        _m0, _m1 = linearRegression(_x, _y)
        temp_m0.append((_m0-m0)**2.0); temp_m1.append((_m1-m1)**2.0)

    print('sigma_m0**2 = ' + str(np.sum(temp_m0)/times) + '\n' + 'sigma_m1**2 = ' + \
        str(np.sum(temp_m1)/times))

    print('+sqrt(sigma_m0**2) = ' + str(np.sqrt(np.sum(temp_m0)/times)) + '\n' + '+sqrt(sigma_m1**2) = ' + \
        str(np.sqrt(np.sum(temp_m1)/times)))
        
def readMJDStart():
    # read MJD start time of every observation
    global mjd_start; mjd_start = []
    mjd_start = np.loadtxt(workpath + 'mjd_start', delimiter = ' ')
    
def makeCSV(sources):
    # make a catalog of all sources with all data available
    _ra = []; _dec = []; _srcid = []; _detid = []
    _uv_brevity_lums = []; _xr_brevity_lums = []; _alphas = []
    _uv_brevity_fluxes = []; _xr_brevity_fluxes = []; _uv_brevity_errors = []; _xr_brevity_errors = []
    _source_numbers = []; _epochs = []
    # flags = [] # flags to show that the sources are synchronized or not
    for i in tqdm(range(len(sources))):
        if(len(sources[i]) == 1):
            _source_numbers.append(int(i+1))
            _epochs.append(int(1))
        elif(len(sources[i]) > 1):
            for j in range(len(sources[i])):
                j=i
                _source_numbers.append(int(j+1))
                _epochs.append(int(len(sources[i])))

    _ra = ra; _dec = dec; _srcid = srcid; _detid = detid
    _uv_brevity_lums = np.loadtxt(workpath + 'uv_brevity_lums', delimiter=' '); _xr_brevity_lums = np.loadtxt(workpath + 'xr_brevity_lums', delimiter=' ')
    _uv_brevity_fluxes = np.loadtxt(workpath + 'uv_brevity_fluxes', delimiter=' '); _xr_brevity_fluxes = np.loadtxt(workpath + 'xr_brevity_fluxes', delimiter=' ')
    _uv_brevity_errors = np.loadtxt(workpath + 'uv_brevity_errors', delimiter=' '); _xr_brevity_errors = np.loadtxt(workpath + 'xr_brevity_errors', delimiter=' ')
    _alphas = np.loadtxt(workpath + 'all_alphas', delimiter=' ')

    data = Table([_ra, _dec, _source_numbers, _epochs, _alphas, _uv_brevity_fluxes, _uv_brevity_errors, _uv_brevity_lums, \
        _xr_brevity_fluxes, _xr_brevity_errors, _xr_brevity_lums, _srcid, _detid], \
            names = ['ra', 'dec', 'source', 'epoch number', 'alpha', 'uv_brevity_flux', 'uv_brevity_error', 'uv_brevity_luminosity', \
                'xr_brevity_flux', 'xr_brevity_error', 'xr_brevity_luminosity', 'SRCID', 'DETID'])
    ascii.write(data, 'meqs.csv', format='csv', fast_writer=False, overwrite = True)

def main():
    # main
    global workpath; workpath = '/Users/snowball/astro/workspace/XMM-OM/'
    lists = fits.open(workpath + 'new.fits')
    flags = multiepoch(lists) # indicate sources which belong to a single source

    sources = groupSource(flags)

    # # readMJDStart() # only use to plot the sequence of the observations

    # readRedshift(lists) # read in all redshifts for every entry
    # # readBrevityLums(); redshiftDist(sources) # only use to plot z-distribution

    # readUVFlux(lists); readUVError(lists); readXrFlux(lists); readXrError(lists)
    # # plotSet() # only use to plot
    # readMeanSED() # read mean SEDs of QSOs from Richards et al. 2006
    # # plotSED(sources) # plot the SED for all sources
    
    # computeUVBrevityFlux(sources); computeXrBrevityFlux(sources) # compute flux for brevity for every source

    # readBrevityLums(); calcAlpha(sources) # only use to calculate all the Alphas

    # readBrevityLums(); doLinearFitting() # only use to do linear fitting to all the Alpha - L data
    # readSNratios(); readBrevityLums(); doSelectedLinearFitting(sources) # only use to do linear fitting to selected (based on S/N ratios) the Alpha - L data

    # plotSetAlpha()
    # readBrevityLums()
    # # plotAlphaOX(sources)
    # # # plt.savefig('multi_uv-xr.pdf', dpi = 2000)
    # # plt.savefig('all_uv-xr.pdf', dpi = 2000)

    # readSNratios()
    # plotSelectedAlphaOX(sources)
    # plt.savefig('multi_selected_Alpha_0.0.pdf', dpi = 2000)

    readRaDec(lists); readID(lists)
    makeCSV(sources)
    

if __name__ == "__main__":
    main()


