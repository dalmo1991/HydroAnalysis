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

This file contains the python code to calculate the hydrological signatures
presented in table 3 of Addor et al. (2017).

This code represent the translation of the R code "hydro_signatures.R"
hosted on Github in the "camels" repository of "naddor".

References
Addor, N., Newman, A. J., Mizukami, N., and Clark, M. P.: The CAMELS data set:
catchment attributes and meteorology for large-sample studies, Hydrol. Earth
Syst. Sci., 21, 5293-5313, https://doi.org/10.5194/hess-21-5293-2017, 2017.

https://github.com/naddor/camels/blob/master/hydro/hydro_signatures.R
"""


from .utils import check_data
import pandas as pd
import numpy as np


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


def calculate_stream_elas(streamflow, quality, precipitation, hydro_year):
    """
    This function calculates the signature "stream_elas". Note that the analysis
    is done also if some values are missing (i.e., they are just removed from
    the time series). This goes against the code of  Addor et al. (2017) that
    only considers time series with only good quality.

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
        Array of streamflow measurements. It is assumed that it is concomitant
        to the streamflow time series.
    hydro_year : np.array
        Array expressing the hydrological year of the measurements. It is
        assumed that it is concomitant to the streamflow time series.

    Returns
    -------
    dict
        'Sawicz' : signature calculated according to Sawicz et al., 2011,
                   HESS
        'Sankarasubramanian' : signature calculated according to
                               Sankarasubramanian et al., 2001, WRR
    """

    good_quality_data = check_data(streamflow=streamflow,
                                   quality=quality,
                                   precipitation=precipitation,
                                   hydro_year=hydro_year)

    if not good_quality_data:
        return None

    # Construct a pandas DataFrame to filter
    data = pd.DataFrame(data=np.array([streamflow, quality, precipitation]).transpose(),
                        index=hydro_year,
                        columns=['Q', 'QC', 'P'])

    # Take care of the quality code
    data = data[data['QC'] == 0]

    # Calculate the global means
    mean_tot = data.mean()

    # Calculate the annual means
    mean_year = data.groupby(data.index).mean()

    # Anomaly computed with respect to previous year (Sawicz et al., 2011, HESS)
    diff_prev_year = mean_year.diff()
    e_sawicz = ((diff_prev_year['Q']/mean_tot['Q']) /
                (diff_prev_year['P']/mean_tot['P'])).median()

    # Anomaly computed with respect to long-term mean (Sankarasubramanian et al., 2001, WRR)
    diff_mean = mean_year - mean_tot
    e_sanka = ((diff_mean['Q']/mean_tot['Q']) /
               (diff_mean['P']/mean_tot['P'])).median()

    sig = {'Sawicz': float(e_sawicz),
           'Sankarasubramanian': float(e_sanka)}

    return sig


def calculate_slope_fdc(streamflow, quality):
    """
    This function calculates the signature "slope_fdc".

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
    dict
        'Sawicz' : Signature calculated according to Sawicz et al 2011, Eq. 3
        'Yadav' : Signature calculated according to Yadav et al 2007, Table 3
        'McMillan' : Signature calculated according to McMillan et al 2017
        'Addor' : Signature calculated according to Addor et al 2017
    """

    good_quality_data = check_data(streamflow=streamflow,
                                   quality=quality)

    if not good_quality_data:
        return None

    quantiles = np.arange(start=0, stop=1.001, step=0.001)*100
    fdc = -np.sort(-np.percentile(streamflow[quality == 0], quantiles))
    q33 = fdc[quantiles == 33.0]
    q66 = fdc[quantiles == 66.0]
    q_median = np.median(streamflow[quality == 0])
    q_mean = streamflow[quality == 0].mean()

    # Calculate signatures
    if (q66 != 0) and (not np.isnan(q66)):
        slope_sawicz = float((np.log(q33) - np.log(q66))/(0.66 - 0.33))
        slope_yadav = float(((q33/q_mean) - (q66/q_mean))/(0.66 - 0.33))
        slope_mcmillan = float(
            (np.log(q33/q_median) - np.log(q66/q_median))/(0.66 - 0.33))
        q33_perc, q66_perc = np.percentile(streamflow[quality == 0], [33, 66])
        slope_addor = float(
            (np.log(q66_perc) - np.log(q33_perc))/(0.66 - 0.33))
    else:
        slope_sawicz = None
        slope_yadav = None
        slope_mcmillan = None
        slope_addor = None

    sig = {'Sawicz': slope_sawicz,
           'Yadav': slope_yadav,
           'McMillan': slope_mcmillan,
           'Addor': slope_addor}

    return sig


def calculate_baseflow_index(streamflow, quality, alpha=0.925, num_filters=3, num_reflect=30, returnBF=False):
    """
    This function calculates the signature "baseflow_index". It follows the
    implementation from Ladson et al., 2013.

    Code: https://raw.githubusercontent.com/TonyLadson/BaseflowSeparation_LyneHollick/master/BFI.R

    Parameters
    ----------
    streamflow : np.array
        Array of streamflow measurements. It is assumed that it represent daily
        data.
    quality : np.array
        Array containing the quality code for the streamflow measurements. It
        is assumed that it is concomitant to the streamflow time series. Data
        with good quality is "0", data with bad quality is "1"
    alpha : float
        Parameter of the Lyne & Hollick filter
    num_filters : int
        Number of filter passages
    num_reflect : int
        Number of time steps to reflect
    returnBF : bool
        True if you want to include the baseflow in the results. Doesn't work
        in case of missing data.

    Returns
    -------
    float
        Value of the signature.
    numpy.ndarray
        Time series of baseflow.
    """

    # Definition of auxiliary functions
    def forward_pass(q, alpha):
        """
        Note: this function assumes that the input is already checked for good
        quality.
        """

        qf = [q[0]]

        for i in range(1, len(q)):
            qf.append(alpha * qf[i-1] + 0.5 * (1 + alpha) * (q[i] - q[i-1]))

        qb = [qt-fl if fl > 0 else qt for (qt, fl) in zip(q, qf)]

        return (qf, qb)

    def backward_pass(qb, alpha):

        # Invert the order of qb
        qb_flipped = np.flip(qb)

        # Call forward_pass with arg q_b inverted
        qf_new_flipped, qb_new_flipped = forward_pass(qb_flipped, alpha)

        # Invert the output of forward_pass
        qf_new = np.flip(qf_new_flipped)
        qb_new = np.flip(qb_new_flipped)

        # Return qb and qf
        return(qf_new, qb_new)

    def missing_data(quality, num_reflect):

        # Find the sequences of homogeneous values (https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi)
        n = len(quality)
        # pairwise unequal (string safe)
        y = np.array(quality[1:] != quality[:-1])
        # must include last element posi
        i = np.append(np.where(y), n - 1)
        lengths = np.diff(np.append(-1, i))           # run lengths
        # True means we have the data
        values = quality[i] == 0

        # Check if we have enough data
        if not np.logical_and(values, lengths > num_reflect).any():
            raise ValueError(
                'Must have at least {} consecutive valid values'.format(num_reflect))

        # Find the end of the sequence of valid data
        index = np.where(np.logical_and(values, lengths > num_reflect))[0]
        # It is plus +1 becouse when I do [:end] I don't take the last value
        ends = np.cumsum(lengths)[index]

        # Find the starts
        new_index = index[np.where(index != 0)[0]]-1
        starts = np.cumsum(lengths)[new_index]
        if 0 in index:
            starts = np.insert(starts, 0, 0)

        return (starts, ends)

    def get_bf(streamflow, alpha, num_filters, num_reflect):

        # Add reflected values
        q_reflect = np.zeros(2*num_reflect + len(streamflow))
        q_reflect[:num_reflect] = np.flip(streamflow[:num_reflect])
        q_reflect[num_reflect:len(q_reflect)-num_reflect] = streamflow
        q_reflect[len(
            q_reflect)-num_reflect:] = np.flip(streamflow[len(streamflow)-num_reflect:])

        # Run the filters
        for i in range(num_filters):
            if i % 2 == 0:  # Forward filter
                if i == 0:  # The input is q_reflect
                    qf, qb = forward_pass(q_reflect, alpha)
                else:  # The input is qb
                    qf, qb = forward_pass(qb, alpha)
            else:  # Backward filter
                qf, qb = backward_pass(qb, alpha)

        # Remove the reflected values
        qf = qf[num_reflect:len(q_reflect)-num_reflect]
        qb = qb[num_reflect:len(q_reflect)-num_reflect]

        return (qf, qb)

    good_quality_data = check_data(streamflow=streamflow,
                                   quality=quality)

    if not good_quality_data:
        return None

    if not isinstance(alpha, float):
        raise TypeError('alpha is of type {}'.format(type(alpha)))
    if not isinstance(num_filters, int):
        raise TypeError('num_filters is of type {}'.format(type(num_filters)))
    if not isinstance(num_reflect, int):
        raise TypeError('num_reflect is of type {}'.format(type(num_reflect)))
    if (alpha < 0) or (alpha > 1):
        raise ValueError(
            'alpha must be between 0 and 1. alpha = {}'.format(alpha))
    if (num_filters <= 0) or (num_filters % 2 == 0):
        raise ValueError(
            'num_filters must be positive and odd. num_filters = {}'.format(num_filters))
    if num_reflect < 0:
        raise ValueError(
            'num_reflect must be positive. num_reflect = {}'.format(num_reflect))
    if len(streamflow) < num_reflect:
        raise ValueError('the time series must be longer than num_reflect. {} vs {}'.format(
            len(streamflow), num_reflect))
    if not isinstance(returnBF, bool):
        raise TypeError('returnBF is of type {}'.format(type(returnBF)))

    if np.max(quality) == 0:  # Case 1: I have all the data
        qf, qb = get_bf(streamflow, alpha, num_filters, num_filters)
        sig = np.sum(qb)/np.sum(streamflow)
    else:  # Case 2: I have missing data
        starts, ends = missing_data(quality, num_reflect)
        w = []
        bfi = []
        for s, e in zip(starts, ends):
            chopped = streamflow[s:e]
            qf, qb = get_bf(chopped, alpha, num_filters, num_filters)
            bfi.append(np.sum(qb)/np.sum(chopped))
            w.append(len(chopped))

        # Calculate the weighted average
        sig = np.average(a=bfi, weights=w)

    if returnBF and np.max(quality) == 0:
        return (float(sig), qb)
    else:
        return float(sig)


def calculate_hfd_mean(streamflow, quality, hydro_year):
    """
    This function calculates the signature "hfd_mean". We consider only years
    with at least 360 days (this number is hard coded).

    Parameters
    ----------
    streamflow : np.array
        Array of streamflow measurements. It is assumed that it represent daily
        data.
    quality : np.array
        Array containing the quality code for the streamflow measurements. It
        is assumed that it is concomitant to the streamflow time series. Data
        with good quality is "0", data with bad quality is "1"
    hydro_year : np.array
        Array expressing the hydrological year of the measurements. It is
        assumed that it is concomitant to the streamflow time series.

    Returns
    -------
    dict
        'hfd_mean' : mean of the annual signature
        'hfd_std' : standard deviation of the annual signature
    """

    good_quality_data = check_data(streamflow=streamflow,
                                   quality=quality,
                                   hydro_year=hydro_year)

    if not good_quality_data:
        return None

    # Define the function to apply
    def calculate_hfd(x):
        x = x['Q'].values
        if len(x) < 360:
            return np.nan
        else:
            # The +1 is needed to get the same definition of Addor
            return sum(x.cumsum() < 0.5*x.sum()) + 1

    # Construct a pandas DataFrame to filter
    data = pd.DataFrame(data=np.array([streamflow, quality]).transpose(),
                        index=hydro_year,
                        columns=['Q', 'QC'])

    data = data[data['QC'] == 0]

    # Apply calculate_hfd
    hfd = data.groupby(data.index).apply(calculate_hfd)

    # Calculate the signature
    hfd_mean = float(hfd.mean())
    hfd_std = float(hfd.std())

    sig = {'hfd_mean': hfd_mean,
           'hfd_std': hfd_std}

    return sig


def calculate_percentile(streamflow, quality, percentile):
    """
    Generic function to calculate the percentile signature.

    Parameters
    ----------
    streamflow : np.array
        Array of streamflow measurements. It is assumed that it represent daily
        data.
    quality : np.array
        Array containing the quality code for the streamflow measurements. It
        is assumed that it is concomitant to the streamflow time series. Data
        with good quality is "0", data with bad quality is "1"
    percentile : float
        Percentile to calculate. It must be between 0 and 100.

    Returns
    -------
    float
        Value of the signature.
    """

    good_quality_data = check_data(streamflow=streamflow,
                                   quality=quality)

    if not good_quality_data:
        return None

    sig = np.percentile(streamflow[quality == 0], percentile)

    return float(sig)


def calculate_q_5(streamflow, quality):
    """
    This function calculates the signature "q_5".

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

    return calculate_percentile(streamflow=streamflow,
                                quality=quality,
                                percentile=5)


def calculate_q_95(streamflow, quality):
    """
    This function calculates the signature "q_95".

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

    return calculate_percentile(streamflow=streamflow,
                                quality=quality,
                                percentile=95)


def calculate_high_q_freq_dur(streamflow, quality):
    """
    This function calculates the signatures "high_q_freq" and "high_q_dur".

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
    dict
        'hq_freq' : High flow frequency (number per year)
        'hq_dur' : Mean high flow duration
    """

    good_quality_data = check_data(streamflow=streamflow,
                                   quality=quality)

    if not good_quality_data:
        return None

    # Flag the high flows
    hf = streamflow[quality == 0] > 9*np.median(streamflow[quality == 0])

    if any(hf):
        seq = [x[x != 0]
               for x in np.split(hf, np.where(hf == 0)[0]) if len(x[x != 0])]
        dur = np.mean([len(x) for x in seq])
        freq = (hf.sum()/len(streamflow[quality == 0]))*365.25
        freq = float(freq)
        dur = float(dur)
    else:
        freq = None
        dur = None

    sig = {'hq_freq': freq,
           'hq_dur': dur}

    return sig


def calculate_low_q_freq_dur(streamflow, quality):
    """
    This function calculates the signatures "low_q_freq" and "low_q_dur".

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
    dict
        'lq_freq' : Low flow frequency (number per year)
        'lq_dur' : Mean low flow duration
    """

    good_quality_data = check_data(streamflow=streamflow,
                                   quality=quality)

    if not good_quality_data:
        return None

    # Flag the high flows
    lf = streamflow[quality == 0] <= 0.2*streamflow[quality == 0].mean()

    if any(lf):
        seq = [x[x != 0]
               for x in np.split(lf, np.where(lf == 0)[0]) if len(x[x != 0])]
        dur = np.mean([len(x) for x in seq])
        freq = (lf.sum()/len(streamflow[quality == 0]))*365.25
        freq = float(freq)
        dur = float(dur)
    else:
        freq = None
        dur = None

    sig = {'lq_freq': freq,
           'lq_dur': dur}

    return sig


def calculate_zero_q_freq(streamflow, quality):
    """
    This function calculates the signature "zero_q_freq".

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
    # Function not implemented by Addor. Doing my version

    good_quality_data = check_data(streamflow=streamflow,
                                   quality=quality)

    if not good_quality_data:
        return None

    filtered_streamflow = streamflow[quality == 0]
    sig = len(
        filtered_streamflow[filtered_streamflow == 0])/len(streamflow[quality == 0])

    return float(sig)
