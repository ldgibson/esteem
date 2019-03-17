import numpy as np
from scipy.misc import comb, factorial2

# from .basissets import BasisSet


def overlap_primitive(a, b, alpha, beta, A, B):
    """Calculates the overlap between two unnormalized primitives.

    Parameters
    ----------
    a, b : array-like of int
    alpha, beta : float
    A, B : array-like of float

    Returns
    -------
    Sp : float"""

    if isinstance(a, np.ndarray):
        pass
    else:
        a = np.array(a)
    if isinstance(b, np.ndarray):
        pass
    else:
        b = np.array(b)
    if isinstance(A, np.ndarray):
        pass
    else:
        A = np.array(A)
    if isinstance(B, np.ndarray):
        pass
    else:
        B = np.array(B)

    p = alpha + beta
    K_AB = np.exp(-alpha * beta / p * np.linalg.norm(A - B) ** 2)
    P = (alpha * A + beta * B) / p

    Sp = (np.pi / p) ** (3 / 2) * K_AB
    for i in range(3):
        Sp *= overlap_int1D(a[i], b[i], A[i], B[i], P[i], p)
    return Sp


def overlap_int1D(aw, bw, Aw, Bw, Pw, p):
    int1D = 0
    for i in range((aw + bw) // 2 + 1):
        int1D += f(aw, bw, Aw, Bw, Pw, i) * \
            factorial2(2 * i - 1) / (2 * p) ** i
    return int1D


def f(aw, bw, Aw, Bw, Pw, i):
    out = 0
    for j in range(max(0, 2 * i - aw), min(2 * i, bw) + 1):
        out += comb(aw, 2 * i - j) * comb(bw, j) * \
            (Pw - Aw) ** (aw - 2 * i + j) * (Pw - Bw) ** (bw - j)
    return out


def int_overlap(basis):
    """Calculates the overlap matrix for a given basis.

    Parameters
    ----------
    basis : BasisSet

    Returns
    -------
    overlap : np.ndarray
        Overlap matrix with shape M x M, where `M = len(basis)`.
    """
    M = len(basis)
    overlap = np.zeros((M, M))
    
    for i, basisfn_i in enumerate(basis):
        # Only calculate diagonal and upper off-diagonal values
        # to reduce number of iterations necessary.
        for j, basisfn_j in enumerate(basis[i:], i):
            ij_overlap = 0
            for prim_i in basisfn_i:
                for prim_j in basisfn_j:
                    ij_overlap += prim_i.d * prim_j.d * \
                        prim_i.N * prim_j.N * \
                        overlap_primitive(prim_i.a, prim_j.a,
                                          prim_i.alpha, prim_j.alpha,
                                          prim_i.A, prim_j.A)
            overlap[i, j] = ij_overlap
            if i != j:
                overlap[j, i] = ij_overlap
            else:
                pass
    return overlap
