"""Generate the docstrings for the diagnostic routines
"""


docstrings = {}

docstrings['diag_effsample'] = \
    "Calculating the effective sample size of a particle filter.\n" \
    "\n    " \
    "Based on [1]_, it is defined as the" \
    "\n    " \
    "inverse of the sum of the squared particle filter weights:" \
    "\n    " \
    r":math:`N_{eff} = \frac{1}{\sum_{i=1}^{N} w_i^2}`" \
    "\n    " \
    r"where :math:`w_i` is the weight of particle with index i." \
    "\n    " \
    r"and :math:`N` is the number of particles." \
    "\n\n    " \
    r"If the :math:`N_{eff}=N`, all weights are identical," \
    "\n    " \
    "and the filter has no influence on the analysis." \
    "\n    "  \
    r"If :math:`N_{eff}=0`, the filter is collapsed.""\n\n    " \
    "This is typically called during the analysis step" \
    "\n    "  \
    "of a particle filter,\n    "\
    "e.g. in the analysis step of NETF and LNETF.\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Doucet, A., de Freitas, N., Gordon, N. (2001). \n    "\
    "       An Introduction to Sequential Monte Carlo Methods." \
    "\n    " \
    "       In: Doucet, A., de Freitas, N., Gordon, N. (eds)" \
    "\n    " \
    "       Sequential Monte Carlo Methods in Practice.\n    " \
    "       Statistics for Engineering and Information Science." \
    "\n    " \
    "       Springer, New York, NY." \
    "\n    " \
    "       https://doi.org/10.1007/978-1-4757-3437-9_1"

docstrings['diag_ensstats'] = \
    "Computing the skewness and kurtosis of" \
    "\n    " \
    "the ensemble of a given element of the state vector.\n" \
    "\n    " \
    "The definition used for kurtosis follows that used by [1]_.\n" \
    "\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Lawson, W. G., & Hansen, J. A. (2004).\n    " \
    "       Implications of stochastic and deterministic" \
    "\n    " \
    "       filters as ensemble-based\n    "\
    "       data assimilation methods in varying regimes" \
    "\n    " \
    "       of error growth.\n    "\
    "       Monthly weather review, 132(8), 1966-1981."

docstrings['diag_histogram'] = \
    "Computing the rank histogram of an ensemble.\n" \
    "\n    " \
    "A rank histogram is used to diagnose\n    " \
    "the reliability of the ensemble [1]_.\n    " \
    "A perfectly reliable ensemble should have\n    " \
    "a uniform rank histogram.\n\n    " \
    "The function can be called in the\n    " \
    "pre/poststep routine of PDAF\n    "\
    "both before and after the analysis step\n    " \
    "to collect the histogram information.\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Hamill, T. M. (2001).\n    " \
    "       Interpretation of rank histograms\n    " \
    "       for verifying ensemble forecasts.\n    " \
    "       Monthly Weather Review, 129(3), 550-560."

docstrings['diag_CRPS_nompi'] = \
    "A continuous rank probability score for" \
    "\n    " \
    "an ensemble without using MPI parallelisation.\n" \
    "\n    " \
    "The implementation is based on [1]_.\n" \
    "\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Hersbach, H. (2000), \n    "\
    "       Decomposition of the Continuous Ranked Probability" \
    "\n    " \
    "       Score for\n    " \
    "       Ensemble Prediction Systems,\n    " \
    "       Wea. Forecasting, 15, 559–570,\n    " \
    "       doi:10.1175/1520-0434(2000)015<0559:DOTCRP>2.0.CO;2"

docstrings['diag_CRPS'] = \
    "Obtain a continuous rank probability score for an ensemble.\n" \
    "\n    " \
    "The implementation is based on [1]_.\n\n    " \
    "References\n    " \
    "----------\n    " \
    ".. [1] Hersbach, H. (2000), \n    "\
    "       Decomposition of the Continuous Ranked Probability" \
    "\n    " \
    "       Score for\n    " \
    "       Ensemble Prediction Systems,\n    " \
    "       Wea. Forecasting, 15, 559–570,\n    " \
    "       doi:10.1175/1520-0434(2000)015<0559:DOTCRP>2.0.CO;2"
