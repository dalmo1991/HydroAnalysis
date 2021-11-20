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
"""


import numpy as np
import warnings
import inspect


def check_data(**kwargs):
    """
    This function checks if all the input arguments:
    - are 1D np.ndarray
    - have the same shape
    - there are data points with good quality code (0)
    """

    for k in kwargs:
        # Check if array
        if not isinstance(kwargs[k], np.ndarray):
            raise TypeError('{} is of type {}'.format(k, type(kwargs[k])))
        # Check if shape
        if len(kwargs[k].shape) != 1:
            raise ValueError(
                '{} must be 1D. Shape :  {}'.format(k, kwargs[k].shape))

    for k1 in kwargs:
        for k2 in kwargs:
            if k1 == k2:
                continue
            if kwargs[k1].shape != kwargs[k2].shape:
                raise ValueError('{} and {} have different shape: {}, {}'.format(k1,
                                                                                 k2,
                                                                                 kwargs[k1].shape,
                                                                                 kwargs[k2].shape))

    # Check if at least some data have good quality
    good_quality_data = True

    if 'quality' in kwargs:
        if sum(kwargs['quality'] == 0) < 1:
            good_quality_data = False
            warnings.warn('Skipped because of no data')

    return good_quality_data


def calculate_hydro_year(date, first_month=10):
    """
    This function calculates the hydrological year from a date. The
    hydrological year starts on the month defined by the parameter first_month.

    Parameters
    ----------
    date : pandas.core.indexes.datetimes.DatetimeIndex
        Date series
    first_month : int
        Number of the first month of the hydrological year

    Returns
    -------
    numpy.ndarray
        Hydrological year time series
    """

    hydrological_year = date.year
    hydrological_year[date.month >=
                      first_month] = date.year[date.month >= first_month] + 1

    return hydrological_year.values


def calculate_season(date, mapping=None):
    """
    This function calculates the season from a date. The mapping between month
    and season can be defined with the parameter mapping.


    Parameters
    ----------
    date : pandas.core.indexes.datetimes.DatetimeIndex
        Date series
    mapping : dict
        Dictionary with the month (int) as key and the season (str) as value

    Returns
    -------
    numpy.ndarray
        Season
    """

    if mapping is None:
        mapping = {
            1: 'Winter',
            2: 'Winter',
            3: 'Spring',
            4: 'Spring',
            5: 'Spring',
            6: 'Summer',
            7: 'Summer',
            8: 'Summer',
            9: 'Fall',
            10: 'Fall',
            11: 'Fall',
            12: 'Winter'
        }

    season = date.month.map(mapping).values

    return np.array(season)


def check_data_2D(**kwargs):
    """
    This function checks if all the input arguments:
    - are np.ndarray
    - have the same length
    """

    for k in kwargs:
        if not isinstance(kwargs[k], np.ndarray):
            raise TypeError('{} is of type {}'.format(k, type(kwargs[k])))

    for k1 in kwargs:
        for k2 in kwargs:
            if k1 == k2:
                continue
            if kwargs[k1].shape[-1] != kwargs[k2].shape[-1]:
                raise ValueError('{} and {} have different shape: {}, {}'.format(k1,
                                                                                 k2,
                                                                                 kwargs[k1].shape,
                                                                                 kwargs[k2].shape))

    good_quality_data = True

    if 'quality' in kwargs:
        if sum(kwargs['quality'] == 0) < 1:
            good_quality_data = False
            warnings.warn('Skipped because of no data')

    return good_quality_data


def calculate_multiple_metrics(metrics, observed, simulated, quality):
    """
    This function calculates multiple metrics.

    Parameters
    ----------
    metrics : list
        List of functions to apply. The functions must accept only the
        following parameters, respecting the order
        - observed time series
        - simulated time series
        - quality time series
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
        Dictionary with the metrics
    """

    values = {}

    for fun in metrics:
        output = fun(observed, simulated, quality)

        if isinstance(output, float):
            values[fun.__name__] = output
        elif isinstance(output, dict):
            for key in output.keys():
                values['{}_{}'.format(fun.__name__, key)] = output[key]
        else:
            raise TypeError('Output of the function must be either a float or'
                            'a dictionary. Current type {}'.format(type(output)))

    return values


def calculate_multiple_signatures(signatures, streamflow, quality,
                                  precipitation, hydro_year):
    """
    This function calculates multiple metrics.

    Parameters
    ----------
    signatures : list
        List of functions to apply. The functions can accept only the
        following parameters (the order is not required but the name must be
        respected)
        - "streamflow" time series
        - "quality" time series
        - "precipitation" time series
        - "hydro_year" time series
    streamflow : np.array
        Array of streamflow measurements. It is assumed that it represent daily
        data.
    quality : np.array
        Array containing the quality code for the streamflow measurements. It
        is assumed that it is concomitant to the streamflow time series. Data
        with good quality is "0", data with bad quality is "1"
    precipitation : np.array
        Array of streamflow measurements. It is assumed that it is concomitant
        to the streamflow time series.
    hydro_year : np.array
        Array expressing the hydrological year of the measurements. It is
        assumed that it is concomitant to the streamflow time series.

    Returns
    -------
    dict
        Dictionary with the signatures
    """

    values = {}
    possible_arguments_names = ['streamflow',
                                'quality', 'precipitation', 'hydro_year']
    possible_arguments = [streamflow, quality, precipitation, hydro_year]

    for fun in signatures:

        fun_arguments = [p.name for p in inspect.signature(
            fun).parameters.values()]

        kwarg = {}
        for arg_name, arg in zip(possible_arguments_names, possible_arguments):
            if arg_name in fun_arguments:
                kwarg[arg_name] = arg

        output = fun(**kwarg)

        if isinstance(output, float):
            values[fun.__name__] = output
        elif isinstance(output, dict):
            for key in output.keys():
                values['{}_{}'.format(fun.__name__, key)] = output[key]
        else:
            raise TypeError('Output of the function must be either a float or'
                            'a dictionary. Current type {}'.format(type(output)))

    return values
