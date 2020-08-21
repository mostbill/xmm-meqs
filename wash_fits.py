'''
author: snowball@USTC
date: 2019.12.19
last update: 2019.12.19
fuction: select multi- and single- epoch sample from fits
'''
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from progressbar import *


def writeToFits():
    c1 = fits.Column(name='RA', array=ra, format='E')
    c2 = fits.Column(name='DEC', array=dec, format='E')
    c3 = fits.Column(name='SRCID', array=srcid, format='K')
    c4 = fits.Column(name='OBSID', array=obsid, format='10A')
    # c4 = fits.Column(name='SRCDIST', array=srcdist, format='E')
    # c5 = fits.Column(name='FLAG', array=flag, format='K')
    t = fits.BinTableHDU.from_columns([c1, c2, c3, c4])
    t.writeto('3XMM_unique.fits')

workpath = '/Users/snowball/astro/workspace/XMM-OM/'
lists = fits.open(workpath + 'final')

srcid = lists[1].data['SRCID']
obsid = lists[1].data['OBSID_1']
_obsid = lists[1].data['OBS_ID']
# ra = lists[1].data['RA']
# dec = lists[1].data['DEC']

# unique_id = np.unique(srcid, return_index=True)[1]
# srcid = srcid[unique_id]; obsid = obsid[unique_id]; ra = ra[unique_id]; dec = dec[unique_id];

# writeToFits()

print(len(np.unique(srcid)))
print(len(np.unique(_obsid)))
# print(len(np.unique(srcid, return_index=True)[1])) # unique SRCID index
# print(np.unique(srcid))
# print(srcid[np.unique(srcid, return_index=True)[1]])
print(len(np.unique(obsid)))


n, _, _ =plt.hist(np.unique(srcid, return_counts = True)[1], bins = [1,2,3,4,5])
# plt.savefig(workpath + 'result.png')
print(n)