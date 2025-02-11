# selection.py
# Utilities for component selection for lib_tedana and related files
# (c) 2025 Prantik Kundu, PhD

import gzip
import pickle

import numpy as np


def andb(arrs):
    result = np.zeros(arrs[0].shape)
    for aa in arrs:
        result += np.array(aa, dtype=int)
    return result


# Elbow detection
def getelbow(ks, val=False):
    # Elbow using linear projection method - moderate
    ks = np.sort(ks)[::-1]
    nc = ks.shape[0]
    coords = np.array([np.arange(nc), ks])
    p = coords - np.tile(np.reshape(coords[:, 0], (2, 1)), (1, nc))
    b = p[:, -1]
    b_hat = np.reshape(b / np.sqrt((b**2).sum()), (2, 1))
    proj_p_b = p - np.dot(b_hat.T, p) * np.tile(b_hat, (1, nc))
    d = np.sqrt((proj_p_b**2).sum(axis=0))
    k_min_ind = d.argmax()
    k_min = ks[k_min_ind]
    if val:
        return ks[k_min_ind]
    else:
        return k_min_ind


def getelbow2(ks, val=False):
    # Elbow using mean/variance method - conservative
    ks = np.sort(ks)[::-1]
    nk = len(ks)
    ds = np.array(
        [
            (
                ks[nk - 5 - ii - 1]
                > ks[nk - 5 - ii : nk].mean() + 2 * ks[nk - 5 - ii : nk].std()  # noqa
            )
            for ii in range(nk - 5)
        ][::-1],
        dtype=int,
    )
    dsum = []
    c_ = 0
    for d_ in ds:
        c_ = (c_ + d_) * d_
        dsum.append(c_)
    e2 = np.argmax(np.array(dsum))
    elind = np.max([getelbow(ks), e2])
    if val:
        return ks[elind]
    else:
        return elind


def getelbow3(ks, val=False):
    # Elbow using curvature - aggressive
    ks = np.sort(ks)[::-1]
    dKdt = ks[:-1] - ks[1:]
    dKdt2 = dKdt[:-1] - dKdt[1:]
    curv = np.abs((dKdt2 / (1 + dKdt[:-1] ** 2.0) ** (3.0 / 2.0)))
    curv[np.isnan(curv)] = -1 * 10**6
    maxcurv = np.argmax(curv) + 2
    if val:
        return ks[maxcurv]
    else:
        return maxcurv


def getfbounds(ne):
    F05s = [None, None, 18.5, 10.1, 7.7, 6.6, 6.0, 5.6, 5.3, 5.1, 5.0]
    F025s = [None, None, 38.5, 17.4, 12.2, 10, 8.8, 8.1, 7.6, 7.2, 6.9]
    F01s = [None, None, 98.5, 34.1, 21.2, 16.2, 13.8, 12.2, 11.3, 10.7, 10.0]
    return F05s[ne - 1], F025s[ne - 1], F01s[ne - 1]


def _interpolate(a, b, fraction):
    """Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a) * fraction


def scoreatpercentile(a, per, limit=(), interpolation_method="lower"):
    """
    This function is grabbed from scipy

    """
    values = np.sort(a, axis=0)
    if limit:
        values = values[(limit[0] <= values) & (values <= limit[1])]

    idx = per / 100.0 * (values.shape[0] - 1)
    if idx % 1 == 0:
        score = values[int(idx)]
    else:
        if interpolation_method == "fraction":
            score = _interpolate(values[int(idx)], values[int(idx) + 1], idx % 1)
        elif interpolation_method == "lower":
            score = values[int(np.floor(idx))]
        elif interpolation_method == "higher":
            score = values[int(np.ceil(idx))]
        else:
            raise ValueError(
                "interpolation_method can only be 'fraction', " "'lower' or 'higher'"
            )
    return score


def dice(A, B):
    denom = np.array(A != 0, dtype=int).sum(0) + (np.array(B != 0, dtype=int).sum(0))
    if denom != 0:
        AB_un = andb([A != 0, B != 0]) == 2
        numer = np.array(AB_un, dtype=int).sum(0)
        return 2.0 * numer / denom
    else:
        return 0.0


def rankvec(vals):
    asort = np.argsort(vals)
    ranks = np.zeros(vals.shape[0])
    ranks[asort] = np.arange(vals.shape[0]) + 1
    return ranks


def save_variables_to_file(variables, filename):
    with gzip.open(filename, "wb") as f:
        pickle.dump(variables, f)  # type: ignore


def load_variables_from_file(filename):
    with gzip.open(filename, "rb") as f:
        variables = pickle.load(f)  # type: ignore
    return variables
