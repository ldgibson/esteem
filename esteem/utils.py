import numpy as np
from scipy.special import gamma, gammainc


def boys_function(m, T):
    if np.ndim(T) > 0 and not isinstance(T, np.ndarray):
        T = np.array(T)
    else:
        pass
    if np.ndim(m) > 0 and not isinstance(m, np.ndarray):
        m = np.array(m)
    else:
        pass

    mp = m + (1 / 2)

    # Limit for T -> 0
    threshold = 1e-13
    if np.ndim(T) > 0:
        if np.any(np.abs(T) < threshold):
            y = np.zeros(len(T))
            idx = np.where(np.abs(T) < threshold)[0]
            y[idx] = 1 / (2 * m + 1)
            idx = np.where(np.abs(T) >= threshold)[0]
            y[idx] = gamma(mp) * gammainc(mp, T[idx]) / (2 * T[idx] ** (mp))
        else:
            y = gamma(mp) * gammainc(mp, T) / (2 * T ** (mp))
    else:
        if np.abs(T) < threshold:
            y = 1 / (2 * m + 1)
        else:
            y = gamma(mp) * gammainc(mp, T) / (2 * T ** (mp))
    return y
