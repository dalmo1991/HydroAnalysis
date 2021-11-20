# HydroAnalysis

HydroAnalysis is a Python package to calculate indices and metrics useful in the everyday job of an hydrologist.

The package is the result of the re-organization of code that I have used during my research and that I have decided to publish because of its supposed usefulness.

## Generic interface

*Please refer to the docstrings of the single functions for documentation.*

Generally all the functions have the same interface

```python
def calculate_index_name(flux1, quality, flux2, fluxN, other_time_series):
    # code
```

where:
- `flux1`, `flux2`, `fluxN` are fluxes needed to calculate the index/metric
- `quality` is a vector that is used to signal time steps where the quality of the data is not good and, therefore, are not used to calculate the index/metric.
- `other_time_series`: are other time series needed to calculate the index/metric. Example could be the season.

## Files in the package

- meteo_indices.py: contains the functions to calculate the indices related to meteorological data
- metrics.py: contains the functions to calculate the metrics related to the hydrological data
- streamflow_signatures.py: contains the functions to calculate the signatures related to the hydrological data
- utils.py: contains the functions used by the other files

## Documentation

The documentation is work in progress. It is not complete yet but can be found in [ReadTheDocs](https://hydroanalysis.readthedocs.io/en/latest/).

## Examples

Examples of the usage of the package are not available yet. However, the code should be self-explanatory..just read carefully the docstrings.

## Testing

There is no systematic testing implemented here. However, the code has been manually tested against other available code (e.g., R code from Addor et al. (2017)) when possible.

## Your contribution

Please feel free to contribute to the package. If you have any suggestion or you want to contribute to the documentation, please contact the author and have a look at the documentation page explaining how to do it.