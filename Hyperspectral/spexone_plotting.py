#!/usr/bin/env python3
"""Colouring SPEXone reflectance.

:Author: Paul Tol
:Date: 2024-04-19

Corrected reflectance based on Martin Raspaud's 2015 Python version of Jacques
Descloitres (MODIS Rapid Response Project, NASA/GSFC/SSAI) 2010 C version
(1.7.1) of the MODIS algorithm.

https://lance.modaps.eosdis.nasa.gov/
https://gist.github.com/mraspaud/9a0964cca1ca375d97b4

The MODIS Corrected Reflectance Science Processing Algorithm
(CREFL_SPA) performs a simplified atmospheric correction with MODIS
visible, near-infrared, and short-wave infrared reflectance bands
(bands 1 through 16). It corrects for molecular (Rayleigh) scattering
and gaseous absorption (water vapor and ozone) using climatological
values for gas contents, based on the 6S Radiative Transfer Mode. It
requires no real-time input of ancillary data. The algorithm performs
no aerosol correction.
"""
import warnings
import subprocess
import numpy as np
import xarray as xr
from pathlib import Path
from matplotlib import pyplot as plt


def init_lut():
    """CREFL_SPA function."""
    # without aerosols, taur <= 0.4 in all bands everywhere
    tau = np.arange(4000) * 0.0001
    tau[0] = 0.0001
    a = [-0.57721566, 0.99999193, -0.24991055, 0.05519968, -0.00976004, 0.00107857]
    a.reverse()
    fintexp1 = np.polyval(a, tau) - np.log(tau)
    fintexp3 = (np.exp(-tau) * (1 - tau) + tau * tau * fintexp1) / 2.0
    sphalb = (3 * tau - fintexp3 * (4 + 2 * tau) + 2 * np.exp(-tau)) / (4 + 3 * tau)
    sphalb[0] = 0
    return sphalb


def chand(phi, muv, mus, taur):
    """CREFL_SPA function.

    :param float phi: Azimuthal difference between sun and observation [rad]
        (phi=0 in backscattering direction).
    :param float mus: Cosine of the sun zenith angle.
    :param float muv: Cosine of the observation zenith angle.
    :param float taur: Molecular optical depth.
    :return: Molecular path reflectance rhoray and two other quantities.
    :rtype: float, float, float
    """
    # With depolarization factor xdep = 0.0279:
    # xfd = (1-xdep/(2-xdep))/(1+2*xdep/(2-xdep)) = 2*(1-xdep)/(2+xdep) = 0.958725775
    xfd = 0.958725775
    xbeta2 = 0.5
    musn = mus[..., np.newaxis]
    muvn = muv[..., np.newaxis]
    phios = phi + np.pi
    xcos1 = 1.0
    xcos2 = np.cos(phios)[..., np.newaxis]
    xcos3 = np.cos(2 * phios)[..., np.newaxis]
    xph1 = 1 + (3 * musn * musn - 1) * (3 * muvn * muvn - 1) * xfd / 8
    xph2 = (
        -xfd
        * xbeta2
        * 1.5
        * musn
        * muvn
        * np.sqrt(1 - musn * musn)
        * np.sqrt(1 - muvn * muvn)
    )
    xph3 = xfd * xbeta2 * 0.375 * (1 - musn * musn) * (1 - muvn * muvn)
    p = np.array(
        [
            np.ones_like(mus),
            mus + muv,
            mus * muv,
            mus * mus + muv * muv,
            mus * mus * muv * muv,
        ]
    )
    as01 = np.array([0.33243832, 0.16285370, -0.30924818, -0.10324388, 0.11493334])
    as02 = np.array(
        [-6.777104e-02, 1.577425e-03, -1.240906e-02, 3.241678e-02, -3.503695e-02]
    )
    as01 = as01[:, np.newaxis, np.newaxis]
    as02 = as02[:, np.newaxis, np.newaxis]
    fs01 = np.sum(p * as01, axis=0)[..., np.newaxis]
    fs02 = np.sum(p * as02, axis=0)[..., np.newaxis]
    xlntaur = np.log(taur)
    fs0 = fs01 + fs02 * xlntaur
    fs1 = 0.19666292 - 5.439061e-02 * xlntaur
    fs2 = 0.14545937 - 2.910845e-02 * xlntaur
    trdown = np.exp(-taur / musn)
    trup = np.exp(-taur / muvn)
    xitm1 = (1 - trdown * trup) / (4 * (musn + muvn))
    xitm2 = (1 - trdown) * (1 - trup)
    xitot1 = xph1 * (xitm1 + xitm2 * fs0)
    xitot2 = xph2 * (xitm1 + xitm2 * fs1)
    xitot3 = xph3 * (xitm1 + xitm2 * fs2)
    rhoray = xitot1 * xcos1 + xitot2 * xcos2 * 2 + xitot3 * xcos3 * 2
    return rhoray, trdown, trup


def atm_variables(mus, muv, phi, height, modis_bands):
    """CREFL_SPA function."""
    # modis_bands are zero-based
    MAXAIRMASS = 18
    SCALEHEIGHT = 8000
    aH2O = np.array(
        [
            -5.60723,
            -5.25251,
            0,
            0,
            -6.29824,
            -7.70944,
            -3.91877,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
    )[modis_bands]
    bH2O = np.array(
        [
            0.820175,
            0.725159,
            0,
            0,
            0.865732,
            0.966947,
            0.745342,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
    )[modis_bands]
    # Values of aO3 and taur0 in (one-based) bands 9-16 below provided by B Murch and C Hu Univ South Florida IMaRS, obtained from SEADAS.
    # Jacques Descloitres values for bands 1-8 retained, but the SEADAS values are shown in commented out lines.
    # aO3 = np.array([0.0711, 0.00313, 0.0104, 0.0930, 0, 0, 0, 0.00244, 0.00383, 0.0225, 0.0663, 0.0836, 0.0485, 0.0395, 0.0119, 0.00263])[modis_bands]
    aO3 = np.array(
        [
            0.0715289,
            0,
            0.00743232,
            0.089691,
            0,
            0,
            0,
            0.001,
            0.00383,
            0.0225,
            0.0663,
            0.0836,
            0.0485,
            0.0395,
            0.0119,
            0.00263,
        ]
    )[modis_bands]
    # taur0 = np.array([0.0507, 0.0164, 0.1915, 0.0948, 0.0036, 0.0012, 0.0004, 0.3109, 0.2375, 0.1596, 0.1131, 0.0994, 0.0446, 0.0416, 0.0286, 0.0155])[modis_bands]
    taur0 = np.array(
        [
            0.05100,
            0.01631,
            0.19325,
            0.09536,
            0.00366,
            0.00123,
            0.00043,
            0.3139,
            0.2375,
            0.1596,
            0.1131,
            0.0994,
            0.0446,
            0.0416,
            0.0286,
            0.0155,
        ]
    )[modis_bands]
    UO3 = 0.319
    UH2O = 2.93
    sphalb0 = init_lut()
    musn = mus[..., np.newaxis]
    muvn = muv[..., np.newaxis]
    m = np.ma.masked_greater(1 / musn + 1 / muvn, MAXAIRMASS)
    psurfratio = np.exp(-height / SCALEHEIGHT)
    taur = taur0[np.newaxis, np.newaxis, :] * psurfratio[..., np.newaxis]
    rhoray, trdown, trup = chand(phi, muv, mus, taur)
    sphalb = sphalb0[(taur / 0.0001 + 0.5).astype(int)]
    sphalb = np.ma.masked_greater_equal(taur, 4000 * 0.0001)
    Ttotrayu = ((2 / 3 + muvn) + (2 / 3 - muvn) * trup) / (4 / 3 + taur)
    Ttotrayd = ((2 / 3 + musn) + (2 / 3 - musn) * trdown) / (4 / 3 + taur)
    tO3 = np.exp(-m * UO3 * aO3)
    tO3[:, :, aO3 == 0] = 1
    tH2O = np.exp(-np.exp(aH2O + bH2O * np.log(m * UH2O)))
    tH2O[:, :, bH2O == 0] = 1
    TtotraytH2O = np.ma.masked_invalid(Ttotrayu * Ttotrayd * tH2O)
    # tO2 = np.exp(-m*aO2) # commented out in original code
    tO2 = 1
    tOG = np.ma.masked_invalid(tO3 * tO2)
    return sphalb, rhoray, TtotraytH2O, tOG


def correct_reflectance(refl, sphalb, rhoray, TtotraytH2O, tOG):
    """CREFL_SPA function."""
    corr_refl = (refl / tOG - rhoray) / TtotraytH2O
    corr_refl /= 1 + corr_refl * sphalb
    corr_refl = np.clip(corr_refl, -0.01, 1.6)
    return corr_refl


def modis_colours(
    rad,
    irrad,
    waves,
    sza,
    saa,
    vza,
    vaa,
    height,
    *,
    reflect=False,
    cloud=True,
    scaling=None,
    enhancement=None,
):
    """Convert spectral radiance to RGB colours as in MODIS true colour images.

    :param ndarray rad: Array of spectral radiance [W m-2 sr-1 um-1] with last
        dimension: spectral sample.
    :param ndarray irrad: Array of spectral irradiance [W m-2 um-1], corrected
        for the Sun-Earth distance, with last dimension: spectral sample.
    :param ndarray waves: List of wavelength [nm] with one dimension: spectral
        sample.
    :param ndarray sza: Solar zenith angle [deg].
    :param ndarray saa: Solar azimuth angle [deg].
    :param ndarray vza: Viewing zenith angle [deg].
    :param ndarray vaa: Viewing azimuth angle [deg].
    :param ndarray height: Ground height [m].
    :param bool cloud: Brightness is enhanced less when there are clouds (True,
        default) and a bit more for cloud-free images.
    :param bool scaling: Optional scaling parameter, None (default) for MODIS
        value depending on cloud parameter. Value 1 does nothing.
    :param bool enhancement: Optional enhancement parameter, None (default) for
        MODIS value depending on cloud parameter. Value 0 does nothing.
    :param bool oci: flag to omit conversion to reflection when using oci
        reflectance as input
    :return: Array of colour components with dimensions like rad, but the last
        dimension consists instead of the components R, G, B.
    :rtype ndarray

    :note: Afterwards, the MODIS algorithm enhances lightness 5% and saturation
        25%, and sharpens the image by applying the ImageMagick command
        "-modulate 105,125 -sharpen 3x1" to the exported image.
    """
    mus = np.cos(np.deg2rad(np.clip(sza, 0, 86.5)))
    muv = np.cos(np.deg2rad(vza))
    if reflect:
        refl = rad
    else:
        refl = rad * np.pi / (irrad * mus[..., np.newaxis])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # Mean of empty slice is okay
        red = np.nanmean(refl[..., (waves >= 620) & (waves <= 670)], axis=-1)
        green = np.nanmean(refl[..., (waves >= 545) & (waves <= 565)], axis=-1)
        blue = np.nanmean(refl[..., (waves >= 459) & (waves <= 479)], axis=-1)
        refl_bands = np.dstack((red, green, blue))
        atm_var = atm_variables(
            mus, muv, np.deg2rad(saa - vaa), height, modis_bands=[0, 3, 2]
        )
        refl_bands = correct_reflectance(refl_bands, *atm_var)
        if scaling is None:
            scaling = 1 / 1.1 if cloud else 1 / 0.8
        if enhancement is None:
            enhancement = 3.7 if cloud else 4.5
        clrs = refl_bands * scaling
        # if a pixel colour is above 1, scale all colours of that pixel until largest is 1
        mx = np.clip(np.nanmax(clrs, axis=2), 1, None)
        clrs /= mx[..., np.newaxis]
        clrs = np.clip(clrs, 0, None)
    clrs = (1 + enhancement) * clrs / (1 + enhancement * clrs)  # brightness enhancement
    clrs = clrs.filled(fill_value=np.nan)  # masked values white
    return clrs

def spexone_rad_to_modis(l1c_dataset):

    # rgb dataset with as coordinates longitude, latitude and sensor_view_angle
    clrs = [] 
    
    # loop over the viewports and update the values in the rgb dataset
    for vp in l1c_dataset["sensor_view_angle"]:
        clrs += [modis_colours(
            l1c_dataset["i_polsample"].where(l1c_dataset.sensor_view_angle == vp, drop=True).squeeze().values,
            (l1c_dataset["polarization_f0"] / l1c_dataset.attrs["sun_earth_distance"]**2).where(l1c_dataset.sensor_view_angle == vp, drop=True).squeeze().values,
            l1c_dataset["polarization_wavelength"].where(l1c_dataset.sensor_view_angle == vp, drop=True).squeeze().values,
            l1c_dataset[f"solar_zenith_angle"].where(l1c_dataset.sensor_view_angle == vp, drop=True).squeeze().values,
            l1c_dataset[f"solar_azimuth_angle"].where(l1c_dataset.sensor_view_angle == vp, drop=True).squeeze().values,
            l1c_dataset[f"sensor_zenith_angle"].where(l1c_dataset.sensor_view_angle == vp, drop=True).squeeze().values,
            l1c_dataset[f"sensor_azimuth_angle"].where(l1c_dataset.sensor_view_angle == vp, drop=True).squeeze().values,
            l1c_dataset["height"].values,
            scaling=0.7,
            enhancement=6,
        )]  # adjusted parameters for SPEXone

        
    rad_rgb = xr.Dataset(
        data_vars=dict(
            rgb=(["number_of_views", "bins_along_track", "bins_across_track", "color_dim"], np.array(clrs))
        ),
        coords=dict(
            longitude=(["bins_along_track", "bins_across_track"], l1c_dataset["longitude"].values),
            latitude=(["bins_along_track", "bins_across_track"], l1c_dataset["latitude"].values),
            sensor_view_angle=("number_of_views", l1c_dataset["sensor_view_angle"].values),
            color_dim=("color_dim", ["red", "green", "blue"]),
        )
    )

    return rad_rgb