"""
Copyright 2021 Marco Dal Molin et al.

This file is part of the HydroAnalysis modelling framework. For details about
it, visit the page https://hydroanalysis.readthedocs.io/

HydroAnalysis is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HydroAnalysis is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with HydroAnalysis. If not, see <https://www.gnu.org/licenses/>.

CODED BY: Marco Dal Molin
DESIGNED BY: Marco Dal Molin

********************************************************************************

This file contains the python code to calculate the meteorological indices
presented in table 2 of Addor et al. (2017).

This code represent the translation of the R code "clim_indices.R"
hosted on Github in the "camels" repository of "naddor".

References
Addor, N., Newman, A. J., Mizukami, N., and Clark, M. P.: The CAMELS data set:
catchment attributes and meteorology for large-sample studies, Hydrol. Earth
Syst. Sci., 21, 5293-5313, https://doi.org/10.5194/hess-21-5293-2017, 2017.
https://github.com/naddor/camels/blob/master/clim/clim_indices.R

"""

import numpy as np
import pandas as pd
from .utils import check_data, calculate_season
from scipy import stats
from scipy.optimize import least_squares


def calculate_p_mean(precipitation, quality):
    """
    This function calculates the signature "mean daily precipitation".

    Parameters
    ----------
    precipitation : np.array
        Array of precipitation measurements. It is assumed that it represent
        daily data.
    quality : np.array
        Array containing the quality code for the precipitation measurements. It
        is assumed that it is concomitant to the precipitation time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the signature.
    """

    good_quality_data = check_data(precipitation=precipitation,
                                   quality=quality)

    if not good_quality_data:
        return None

    sig = precipitation[quality == 0].mean()

    return float(sig)


def calculate_pet_mean(pet, quality):
    """
    This function calculates the signature "mean daily pet".

    Parameters
    ----------
    pet : np.array
        Array of pet measurements. It is assumed that it represent daily data.
    quality : np.array
        Array containing the quality code for the pet measurements. It is
        assumed that it is concomitant to the pet time series. Data with good
        quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the signature.
    """

    good_quality_data = check_data(pet=pet,
                                   quality=quality)

    if not good_quality_data:
        return None

    sig = pet[quality == 0].mean()

    return float(sig)


def calculate_aridity(precipitation, pet, quality):
    """
    This function calculates the signature "aridity".

    Parameters
    ----------
    precipitation : np.array
        Array of precipitation measurements. It is assumed that it represent
        daily data.
    pet : np.array
        Array containing the pet time series. It is assumed that it is
        concomitant to the precipitation time series.
    quality : np.array
        Array containing the quality code for the precipitation measurements. It
        is assumed that it is concomitant to the precipitation time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the signature.
    """

    good_quality_data = check_data(precipitation=precipitation,
                                   pet=pet,
                                   quality=quality)

    if not good_quality_data:
        return None

    sig = calculate_pet_mean(pet, quality) /\
        calculate_p_mean(precipitation, quality)

    return float(sig)


def calculate_p_seasonality(precipitation, quality, date, temperature):
    """
    This function calculates the signature "p_seasonality".

    Parameters
    ----------
    precipitation : np.array
        Array of precipitation measurements. It is assumed that it represent
        daily data.
    quality : np.array
        Array containing the quality code for the precipitation measurements. It
        is assumed that it is concomitant to the precipitation time series. Data
        with good quality is "0", data with bad quality is "1"
    date : pandas.core.indexes.datetimes.DatetimeIndex
        Date series. It is assumed that it is concomitant to the precipitation
        time series
    temperature : np.array
        Array containing the temperature time series. It is assumed that it is
        concomitant to the precipitation time series.

    Returns
    -------
    float
        Value of the signature.
    """

    good_quality_data = check_data(precipitation=precipitation,
                                   quality=quality,
                                   temperature=temperature)

    if not good_quality_data:
        return None

    if len(precipitation) != len(date):
        raise ValueError('precipitation and date have different length: {} vs {}'.format(
            len(precipitation), len(date)))

    # Calculate the day of the year
    t_julian = date.dayofyear-1  # 1 Jan is zero

    # Create a DataFrame
    prec = pd.DataFrame(data=np.array([precipitation, quality]).transpose(),
                        index=date,
                        columns=['P', 'QC'])

    mean_month_prec = prec[prec.QC == 0].groupby(
        prec.index[prec.QC == 0].month).mean()

    # Get a first guess of the phase -> month with the most precipitation
    sp_first_guess = 90 - (mean_month_prec.idxmax()['P'] - 1)*30
    sp_first_guess = sp_first_guess + 360 if sp_first_guess < 0 else sp_first_guess

    # Fit the two sine functions
    def fit_p(pars, x, y):
        prec = y.mean() * (1 + pars[0]*np.sin(2*np.pi*(x-pars[1])/365.25))
        return y - prec

    def fit_t(pars, x, y):
        temp = y.mean() + pars[0]*np.sin(2*np.pi*(x-pars[1])/365.25)
        return y - temp

    prec_pars = least_squares(fun=fit_p,
                              x0=[0.4, sp_first_guess],
                              bounds=([-1, 0], [1, 365.25]),
                              args=(t_julian[quality == 0], precipitation[quality == 0]))

    temp_pars = least_squares(fun=fit_t,
                              x0=[5, 270],
                              args=(t_julian[quality == 0], temperature[quality == 0]))

    # Explicit the parameters
    delta_p = prec_pars.x[0]
    sp = prec_pars.x[1]
    delta_t = temp_pars.x[0]
    st = temp_pars.x[1]

    sig = delta_p * np.sign(delta_t) * np.cos(2 * np.pi * (sp - st / 365.25))

    return(float(sig))


def calculate_frac_snow(precipitation, temperature, quality, threshold=0.0):
    """
    This function calculates the signature "frac_snow".

    Parameters
    ----------
    precipitation : np.array
        Array of precipitation measurements. It is assumed that it represent
        daily data.
    temperature : np.array
        Array containing the temperature time series. It is assumed that it is
        concomitant to the precipitation time series.
    quality : np.array
        Array containing the quality code for the precipitation measurements. It
        is assumed that it is concomitant to the precipitation time series. Data
        with good quality is "0", data with bad quality is "1"
    threshold : float
        Threshold temperature to separate snow from rainfall

    Returns
    -------
    float
        Value of the signature.
    """

    good_quality_data = check_data(precipitation=precipitation,
                                   temperature=temperature,
                                   quality=quality)

    if not good_quality_data:
        return None

    if not isinstance(threshold, float):
        raise TypeError('threshold is of type {}'.format(type(threshold)))

    p = precipitation[quality == 0]
    t = temperature[quality == 0]

    sig = np.sum(p[t < threshold])/np.sum(p)

    return float(sig)


def calculate_high_prec_freq_time(precipitation, quality, date):
    """
    This function calculates the signature "high_prec_freq".

    Parameters
    ----------
    precipitation : np.array
        Array of precipitation measurements. It is assumed that it represent
        daily data.
    quality : np.array
        Array containing the quality code for the precipitation measurements. It
        is assumed that it is concomitant to the precipitation time series. Data
        with good quality is "0", data with bad quality is "1"
    date : pandas.core.indexes.datetimes.DatetimeIndex
        Date series. It is assumed that it is concomitant to the precipitation
        time series

    Returns
    -------
    dict
        'hp_freq' : High precipitation frequency (number per year)
        'hp_dur' : Mean high precipitation duration
        'hp_time' : Season with more high precipitation
    """

    good_quality_data = check_data(precipitation=precipitation,
                                   quality=quality)

    if not good_quality_data:
        return None

    if len(precipitation) != len(date):
        raise ValueError('precipitation and date have different length: {} vs {}'.format(
            len(precipitation), len(date)))

    # Calculate the season
    seasons = calculate_season(date=date)

    # Flag the high flows
    hp = precipitation[quality == 0] > 5*np.mean(precipitation[quality == 0])
    high_season = seasons[quality == 0][hp]

    if any(hp):
        seq = [x[x != 0]
               for x in np.split(hp, np.where(hp == 0)[0]) if len(x[x != 0])]
        dur = np.mean([len(x) for x in seq])
        freq = (hp.sum()/len(precipitation[quality == 0]))*365.25
        freq = float(freq)
        dur = float(dur)
        time = stats.mode(high_season).mode[0]
    else:
        freq = None
        dur = None

    sig = {'hp_freq': freq,
           'hp_dur': dur,
           'hp_time': time}

    return sig


def calculate_low_prec_freq_time(precipitation, quality, date):
    """
    This function calculates the signature "low_prec_freq".

    Parameters
    ----------
    precipitation : np.array
        Array of precipitation measurements. It is assumed that it represent
        daily data.
    quality : np.array
        Array containing the quality code for the precipitation measurements. It
        is assumed that it is concomitant to the precipitation time series. Data
        with good quality is "0", data with bad quality is "1"
    date : pandas.core.indexes.datetimes.DatetimeIndex
        Date series. It is assumed that it is concomitant to the precipitation
        time series

    Returns
    -------
    dict
        'lp_freq' : Low precipitation frequency (number per year)
        'lp_dur' : Mean low precipitation duration
        'lp_time' : Season with more low precipitation
    """

    good_quality_data = check_data(precipitation=precipitation,
                                   quality=quality)

    if not good_quality_data:
        return None

    if len(precipitation) != len(date):
        raise ValueError('precipitation and date have different length: {} vs {}'.format(
            len(precipitation), len(date)))

    # Calculate the season
    seasons = calculate_season(date=date)

    # Flag the low flows
    lp = precipitation[quality == 0] < 1  # Lower than 1 mm/day
    low_season = seasons[quality == 0][lp]

    if any(lp):
        seq = [x[x != 0]
               for x in np.split(lp, np.where(lp == 0)[0]) if len(x[x != 0])]
        dur = np.mean([len(x) for x in seq])
        freq = (lp.sum()/len(precipitation[quality == 0]))*365.25
        freq = float(freq)
        dur = float(dur)
        time = stats.mode(low_season).mode[0]
    else:
        freq = None
        dur = None

    sig = {'lp_freq': freq,
           'lp_dur': dur,
           'lp_time': time}

    return sig
