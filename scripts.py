import lightkurve as lk
import numpy as np
from urllib.request import HTTPError
from astropy.io import fits


def fetch_hlsps(lc):
    '''Fetches HLSPs, given a lightkurve.KeplerLightCurve

    Parameters
    ----------
    lc : lightkurve.KeplerLightCurve
        Light curve to find high level science products for

    Returns
    -------
    hlsps : lightkurve.LightCurveCollection
        Collection of HLSP light curves
     '''

    if not isinstance(lc, lk.KeplerLightCurve):
        raise ValueError('Please pass a lightkurve.KeplerLightCurve.')

    v_url = ('https://archive.stsci.edu/hlsps/k2sff/'
             'c{0:02d}/{1}/{2:05d}/hlsp_k2sff_k2_lightcurve_{3}-c{0:02d}_kepler_v1_llc.fits'
             ''.format(lc.campaign, (lc.targetid//100000)*100000, lc.targetid - (lc.targetid//100000)*100000, lc.targetid))
    e_url = ('https://archive.stsci.edu/hlsps/everest/v2/'
             'c{0:02d}/{1}/{2:05d}/hlsp_everest_k2_llc_{3}-c{0:02d}_kepler_v2.0_lc.fits'
             ''.format(lc.campaign, (lc.targetid//100000)*100000, lc.targetid - (lc.targetid//100000)*100000, lc.targetid))
    s_url = ('https://archive.stsci.edu/hlsps/k2sc/v2/'
             'c{0:02d}/{1}/hlsp_k2sc_k2_llc_{2}-c{0:02d}_kepler_v2_lc.fits'
             ''.format(lc.campaign, (lc.targetid//100000)*100000, lc.targetid))

    hlsps = []

    try:
        hdu = fits.open(v_url)[1].data
        lc1 = lk.KeplerLightCurve(hdu['T'], hdu['FCOR'], cadenceno=hdu['CADENCENO'], meta={'arclength': hdu['ARCLENGTH']},
                                  channel=lc.channel, campaign=lc.campaign, quarter=lc.quarter,
                                  mission=lc.mission, targetid=lc.targetid, ra=lc.ra, dec=lc.dec,
                                  label='K2SFF (Vanderburg & Johnson)')
        lc1 = lc1[np.isfinite(lc1.time)]
        hlsps.append(lc1)
    except HTTPError:
        pass

    try:
        hdu = fits.open(e_url)[1].data
        lc1 = lk.KeplerLightCurve(hdu['TIME'], hdu['FCOR'], cadenceno=hdu['CADN'], quality=hdu['QUALITY'],
                                  channel=lc.channel, campaign=lc.campaign, quarter=lc.quarter,
                                  mission=lc.mission, targetid=lc.targetid, ra=lc.ra, dec=lc.dec,
                                  label='EVEREST (Luger)')
        lc1 = lc1[((hdu['QUALITY'] & 24) == 0)]

        lc1 = lc1[np.isfinite(lc1.time)]
        hlsps.append(lc1)

    except HTTPError:
        pass

    try:
        hdu = fits.open(s_url)[1].data
        lc1 = lk.KeplerLightCurve(hdu['time'], hdu['flux'], hdu['error'], cadenceno=hdu['cadence'], quality=hdu['quality'],
                                  channel=lc.channel, campaign=lc.campaign, quarter=lc.quarter,
                                  mission=lc.mission, targetid=lc.targetid, ra=lc.ra, dec=lc.dec,
                                  label='K2SC (Aigrain)')
        lc1 = lc1[np.isfinite(lc1.time)]
        hlsps.append(lc1)
    except HTTPError:
        pass

    return lk.LightCurveCollection(hlsps)
