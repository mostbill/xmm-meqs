'''
author: snowball@USTC
date: 2019.12.5
last update: 2020.3.2
fuction: match for multi-epoch sources in XMM-OM SUSS4.1
'''
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np; np.set_printoptions(suppress=True)
from tqdm import tqdm

def match_for_single_epoch(_ra, _dec, src_number):
    #select a single-epoch and match
    _ra = np.array([_ra]); _dec = np.array([_dec])
    source = SkyCoord(ra = _ra*u.degree, dec = _dec*u.degree) # single source
    catalog = SkyCoord(ra = ra*u.degree, dec = dec*u.degree) # whole catalogue
    idx, d2d, d3d = catalog.match_to_catalog_sky(source)
    max_sep = 2 * u.arcsec
    sep_constraint = d2d < max_sep

    #set the flags value, if matched
    for i in range(len(ra)):
        if(flags[i] == 0 and sep_constraint[i]):
            #not matched & not matched to itself
            flags[i] = src_number

    return src_number
    # print(flags)

def multiepoch(lists):
    #readin catalog's ra and dec and numpy-ize
    global ra; global dec; global srcid
    ra = lists[1].data['RA_1'][:]; dec = lists[1].data['DEC_1'][:]; srcid = lists[1].data['SRCID'][:]; mjd_start = lists[1].data['MJD_START']
    np.savetxt("mjd_start", mjd_start, fmt = '%f', delimiter = ' ')
    ra = np.array(ra); dec = np.array(dec)
    global flags # indicate single/multiple entries within the same source
    flags = np.zeros(len(ra), dtype = int) # indicate if this source is already matched; 0 for not matched, SRCID if matched

    src_number = 0
    for i in tqdm(range(len(ra))):
        #update flagss
        if(flags[i] == 0):
            src_number += 1
        src_number = match_for_single_epoch(ra[i], dec[i], src_number)
        #srcdist.append(match_for_single_epoch(ra[i], dec[i], srcid[i]))

    len_unique_flags, counts_unique_flags = np.unique(flags, return_counts = True)

    np.savetxt("flags", flags, fmt = '%d', delimiter = ' ')
    n, _, _ = plt.hist(counts_unique_flags, bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    
    plt.savefig('counts_flags.png')
    plt.clf()

    return(flags)

def main():
    #readin lists information
    workpath = '/Users/snowball/astro/workspace/XMM-OM/'
    lists = fits.open(workpath + 'new.fits')

    multiepoch(lists)

if __name__ == "__main__":
    main()