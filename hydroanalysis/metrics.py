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

This file contains the python code to calculate some metrics to evaluate model
performance.

This code is taken from the following sources:
- R code "hydro_accuracy.R" hosted on Github in the "camels" repository of
  "naddor".
- Paper McInerney et al., 2017

References

https://github.com/naddor/camels/blob/master/hydro/hydro_accuracy.R

McInerney, D., M. Thyer, D. Kavetski, J. Lerat, and G. Kuczera (2017),
Improving probabilistic prediction of daily streamflow by identifying Pareto
optimal approaches for modeling heteroscedastic residual errors, Water
Resour. Res., 53, 2199–2239, doi:10.1002/2016WR019168.
"""

import numpy as np
from .utils import check_data, check_data_2D


def calculate_nse(observed, simulated, quality):
    """
    This function calculates the Nash-Sutcliffe efficiency.

    Parameters
    ----------
    observed : np.array
        Array of streamflow measurements.
    simulated : np.array
        Array of streamflow simulations. It is assumed that it is concomitant
        to the observed time series.
    quality : np.array
        Array containing the quality code for the observed measurements. It
        is assumed that it is concomitant to the observed time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the metric.
    """

    good_quality_data = check_data(observed=observed,
                                   simulated=simulated,
                                   quality=quality)

    if not good_quality_data:
        return None

    metric = 1 - np.sum(np.power(simulated[quality == 0] - observed[quality == 0], 2)) /\
        np.sum(np.power(observed[quality == 0] -
               np.mean(observed[quality == 0]), 2))

    return float(metric)


def calculate_rmse(observed, simulated, quality):
    """
    This function calculates the root mean square error.

    Parameters
    ----------
    observed : np.array
        Array of streamflow measurements.
    simulated : np.array
        Array of streamflow simulations. It is assumed that it is concomitant
        to the observed time series.
    quality : np.array
        Array containing the quality code for the observed measurements. It
        is assumed that it is concomitant to the observed time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the metric.
    """

    good_quality_data = check_data(observed=observed,
                                   simulated=simulated,
                                   quality=quality)

    if not good_quality_data:
        return None

    metric = np.sqrt(np.sum(np.power(simulated[quality == 0] - observed[quality == 0], 2)) /
                     len(observed[quality == 0]))

    return float(metric)


def calculate_nrmse(observed, simulated, quality):
    """
    This function calculates the normalized root mean square error (values
    normalized by the mean of the observed values).

    Parameters
    ----------
    observed : np.array
        Array of streamflow measurements.
    simulated : np.array
        Array of streamflow simulations. It is assumed that it is concomitant
        to the observed time series.
    quality : np.array
        Array containing the quality code for the observed measurements. It
        is assumed that it is concomitant to the observed time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the metric.
    """

    # No need to check since it is done by calculate_rmse
    rmse = calculate_rmse(observed=observed,
                          simulated=simulated,
                          quality=quality)

    metric = rmse/np.mean(observed[quality == 0])

    return float(metric)


def calculate_dv(observed, simulated, quality):
    """
    This function calculates the error in the water balance.

    Parameters
    ----------
    observed : np.array
        Array of streamflow measurements.
    simulated : np.array
        Array of streamflow simulations. It is assumed that it is concomitant
        to the observed time series.
    quality : np.array
        Array containing the quality code for the observed measurements. It
        is assumed that it is concomitant to the observed time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the metric.
    """

    good_quality_data = check_data(observed=observed,
                                   simulated=simulated,
                                   quality=quality)

    if not good_quality_data:
        return None

    metric = (np.sum(simulated[quality == 0]) - np.sum(observed[quality == 0])) /\
        np.sum(observed[quality == 0])

    return float(metric)


def calculate_kge(observed, simulated, quality):
    """
    This function calculates the Kling Gupta efficiency

    Parameters
    ----------
    observed : np.array
        Array of streamflow measurements.
    simulated : np.array
        Array of streamflow simulations. It is assumed that it is concomitant
        to the observed time series.
    quality : np.array
        Array containing the quality code for the observed measurements. It
        is assumed that it is concomitant to the observed time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    dict
        'kge' : Kling Gupta efficiency
        'r' : Correlation coefficient
        'alpha' : Ratio between sim and obs standard deviations (calculated
                  dividing by n-1)
        'beta' : Ratio between sim and obs means
    """

    good_quality_data = check_data(observed=observed,
                                   simulated=simulated,
                                   quality=quality)

    if not good_quality_data:
        return None

    r = np.corrcoef(observed[quality == 0], simulated[quality == 0])[1, 0]
    alpha = np.std(simulated[quality == 0], ddof=1) /\
        np.std(observed[quality == 0], ddof=1)
    beta = np.mean(simulated[quality == 0]) /\
        np.mean(observed[quality == 0])

    kge = 1 - np.sqrt(np.power(r-1, 2) +
                      np.power(alpha-1, 2) +
                      np.power(beta-1, 2))

    metrics = {'kge': float(kge),
               'r': float(r),
               'alpha': float(alpha),
               'beta': float(beta)}

    return metrics


def calculate_reliability(modelled, observed, quality):
    """
    This function calculates the reliability.

    Parameters:
    ----------
    modelled : np.array
        2D array of streamflow simulations. Each row is a realization. It is
        assumed that it is concomitant to the observed time series.
    observed : np.array
        Array of streamflow measurements.
    quality : np.array
        Array containing the quality code for the observed measurements. It
        is assumed that it is concomitant to the observed time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the metric.
    """

    good_quality_data = check_data(modelled=modelled,
                                   quality=quality,
                                   observed=observed)

    if not good_quality_data:
        return None

    # Calculate pqq
    pqq = calculate_pqq(modelled, observed, quality)

    # Get the vertical distance from the diagonal
    diagonal = np.linspace(
        start=0, stop=1, num=modelled[:, quality == 0].shape[1])
    vertical_distance = pqq - diagonal

    # Calculate the metric
    metric = (2/modelled[:, quality == 0].shape[1]) * \
        np.sum(np.abs(vertical_distance))

    return float(metric)


def calculate_pqq(modelled, observed, quality):
    """
    This function calculates the pqq plot. It returns an array of values.

    Parameters
    ----------
    modelled : np.array
        2D array of streamflow simulations. Each row is a realization. It is
        assumed that it is concomitant to the observed time series.
    observed : np.array
        Array of streamflow measurements.
    quality : np.array
        Array containing the quality code for the observed measurements. It
        is assumed that it is concomitant to the observed time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    numpy.ndarray
        Value of the metric.
    """

    good_quality_data = check_data(modelled=modelled,
                                   quality=quality,
                                   observed=observed)

    if not good_quality_data:
        return None

    # Calculate pqq
    num_samples = modelled.shape[0]
    pqq = np.sort(
        ((modelled[:, quality == 0] < observed[quality == 0]).sum(axis=0))/num_samples)

    return pqq


def calculate_precision(modelled, observed, quality):
    """
    This function calculates the precision.

    Parameters:
    ----------
    modelled : np.array
        2D array of streamflow simulations. Each row is a realization. It is
        assumed that it is concomitant to the observed time series.
    observed : np.array
        Array of streamflow measurements.
    quality : np.array
        Array containing the quality code for the observed measurements. It
        is assumed that it is concomitant to the observed time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the metric.
    """

    good_quality_data = check_data(modelled=modelled,
                                   quality=quality,
                                   observed=observed)

    if not good_quality_data:
        return None

    num_timesteps = modelled[:, quality == 0].shape[1]
    num = np.sum(np.std(modelled[:, quality == 0], axis=0))/num_timesteps
    den = np.sum(observed[quality == 0])/num_timesteps

    return float(num/den)


def calculate_volumetric_bias(modelled, observed, quality):
    """
    This function calculates the volumetric bias.

    Parameters:
    ----------
    modelled : np.array
        2D array of streamflow simulations. Each row is a realization. It is
        assumed that it is concomitant to the observed time series.
    observed : np.array
        Array of streamflow measurements.
    quality : np.array
        Array containing the quality code for the observed measurements. It
        is assumed that it is concomitant to the observed time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the metric.
    """

    good_quality_data = check_data(modelled=modelled,
                                   quality=quality,
                                   observed=observed)

    if not good_quality_data:
        return None

    mod_mean = np.mean(modelled[:, quality == 0], axis=0)

    metric = np.abs((np.sum(observed[quality == 0]) - np.sum(mod_mean)) /
                    np.sum(observed[quality == 0]))

    return float(metric)
