"""
Copyright 2021 Marco Dal Molin et al.

This file is part of SuperflexPy.

SuperflexPy is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SuperflexPy is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with SuperflexPy. If not, see <https://www.gnu.org/licenses/>.

This file is part of the SuperflexPy modelling framework. For details about it,
visit the page https://superflexpy.readthedocs.io

CODED BY: Marco Dal Molin
DESIGNED BY: Marco Dal Molin

This file contains the implementation of Element classes with different levels
of specialization.
"""

import numpy as np
import warnings

def check_data(**kwargs):
    """
    This function checks if all the input arguments:
    - are 1D np.ndarray
    - have the same shape
    - there are good quality data
    """

    current_length = None

    for k in kwargs:
        # Check if array
        if not isinstance(kwargs[k], np.ndarray):
            raise TypeError('{} is of type {}'.format(k, type(kwargs[k])))
        # Check if shape
        if len(kwargs[k].shape) != 1:
            raise ValueError('{} must be 1D. Shape :  {}'.format(k, kwargs[k].shape))

    # This is not efficient but in most of the cases we are dealing with a few arrays
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

def calculate_hydro_year(date, first_month = 10):
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

    hydrological_year = []

    # TODO: this can probably be done better with pandas
    for el in date:
        if el.month in range(first_month, 13):
            hydrological_year.append(el.year + 1)
        else:
            hydrological_year.append(el.year)

    return np.array(hydrological_year)

def calculate_season(date, mapping = None):
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
                1 : 'Winter',
                2 : 'Winter',
                3 : 'Spring',
                4 : 'Spring',
                5 : 'Spring',
                6 : 'Summer',
                7 : 'Summer',
                8 : 'Summer',
                9 : 'Fall',
                10 : 'Fall',
                11 : 'Fall',
                12 : 'Winter'
            }

    # TODO: this can probably be done better with pandas
    season = []
    for el in date:
        season.append(mapping[el.month])

    return np.array(season)