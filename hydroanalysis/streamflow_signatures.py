"""
This file contains the python code to calculate the hydrological signatures
presented in table 3 of Addor et al. (2017).

This code represent the translation of the R code "hydro_signatures.R"
hosted on Github in the "camels" repository of "naddor".

References
Addor, N., Newman, A. J., Mizukami, N., and Clark, M. P.: The CAMELS data set:
catchment attributes and meteorology for large-sample studies, Hydrol. Earth
Syst. Sci., 21, 5293-5313, https://doi.org/10.5194/hess-21-5293-2017, 2017.

https://github.com/naddor/camels/blob/master/hydro/hydro_signatures.R

Copyright 2021 Marco Dal Molin et al.

This file is part of HydroAnalysis.

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

This file is part of the HydroAnalysis modelling framework. For details about it,
visit the page https://hydroanalysis.readthedocs.io/

CODED BY: Marco Dal Molin
DESIGNED BY: Marco Dal Molin
"""

from .utils import check_data, calculate_hydro_year, calculate_season

def calculate_q_mean(streamflow, quality):
    """
    This function calculates the signature "mean daily discharge".

    Parameters
    ----------
    streamflow : np.array
        Array of streamflow measurements. It is assumed that it represent daily
        data.
    quality : np.array
        Array containing the quality code for the streamflow measurements. It
        is assumed that it is concomitant to the streamflow time series. Data
        with good quality is "0", data with bad quality is "1"

    Returns
    -------
    float
        Value of the signature.
    """

    good_quality_data = check_data(
        streamflow=streamflow,
        quality=quality
    )

    if not good_quality_data:
        return None

    sig = streamflow[quality == 0].mean()

    return float(sig)

def calculate_runoff_ratio(streamflow, quality, precipitation):
    """
    This function calculates the signature "runoff_ratio".

    Parameters
    ----------
    streamflow : np.array
        Array of streamflow measurements. It is assumed that it represent daily
        data.
    quality : np.array
        Array containing the quality code for the streamflow measurements. It
        is assumed that it is concomitant to the streamflow time series. Data
        with good quality is "0", data with bad quality is "1"
    precipitation : np.array
        Array of precipitation. It is assumed that it is concomitant to the
        streamflow time series.

    Returns
    -------
    float
        Value of the signature.
    """

    good_quality_data = check_data(
        streamflow=streamflow,
        quality=quality,
        precipitation=precipitation
    )

    if not good_quality_data:
        return None

    sig = streamflow[quality == 0].mean() / precipitation[quality == 0].mean()

    return float(sig)