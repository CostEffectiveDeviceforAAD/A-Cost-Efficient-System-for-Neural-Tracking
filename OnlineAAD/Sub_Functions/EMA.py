################################
## Exponential Moving Average ##
################################

import numpy as np


def EMA(corr_j, corr_t, EmaCorr_j, EmaCorr_t, i, tr):

    if i == 0:
        EmaCorr_j[tr - 14, i] = np.array(corr_j)
        EmaCorr_t[tr - 14, i] = np.array(corr_t)
    else:
        weight = 2 / (i + 2)
        EmaCorr_j[tr - 14, i] = np.add((weight * corr_j), ((1 - weight) * EmaCorr_j[tr - 14, i - 1]))
        EmaCorr_t[tr - 14, i] = np.add((weight * corr_t), ((1 - weight) * EmaCorr_t[tr - 14, i - 1]))

    return EmaCorr_j, EmaCorr_t
