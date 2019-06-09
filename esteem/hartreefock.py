import numpy as np
# from scipy.misc import comb
from scipy.signal import convolve2d
from scipy.special import comb, factorial2

from .utils import boys_function


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

    # if isinstance(a, np.ndarray) and not isinstance(a, np.matrix):
        # pass
    # else:
        # a = np.array(a)
    # if isinstance(b, np.ndarray) and not isinstance(a, np.matrix):
        # pass
    # else:
        # b = np.array(b)
    # if isinstance(A, np.ndarray) and not isinstance(a, np.matrix):
        # print(type(A))
        # pass
    # else:
        # A = np.array(A)
    # if isinstance(B, np.ndarray) and not isinstance(a, np.matrix):
        # print(type(B))
        # pass
    # else:
        # B = np.array(B)

    a = np.array(a).reshape(3)
    b = np.array(b).reshape(3)
    A = np.array(A).reshape(3)
    B = np.array(B).reshape(3)

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


def kinenergy_primitive(a, b, alpha, beta, A, B):
    Tp = 0
    for i in range(3):
        Tp += kinenergy_int1D(i, a, b, alpha, beta, A, B)
    return Tp


def kinenergy_int1D(w, a, b, alpha, beta, A, B):
    two = np.eye(len(b), dtype=np.int) * 2

    int1D = beta * (2 * b[w] + 1) *\
        overlap_primitive(a, b, alpha, beta, A, B) - \
        2 * beta ** 2 *\
        overlap_primitive(a, b + two[w, :], alpha, beta, A, B)

    # if np.all(b - two[w, :] >= 0):
    if b[w] >= 2:
        int1D -= (1 / 2) * b[w] * (b[w] - 1) * \
            overlap_primitive(a, b - two[w, :], alpha, beta, A, B)
    else:
        pass

    return int1D


def int_kinenergy(basis):
    """Calculates the kinenergy matrix for a given basis.

    Parameters
    ----------
    basis : BasisSet

    Returns
    -------
    kinenergy : np.ndarray
        Kinetic energy matrix with shape M x M, where `M = len(basis)`.
    """
    M = len(basis)
    kinenergy = np.zeros((M, M))

    for i, basisfn_i in enumerate(basis):
        # Only calculate diagonal and upper off-diagonal values
        # to reduce number of iterations necessary.
        for j, basisfn_j in enumerate(basis[i:], i):
            ij_kinenergy = 0
            for prim_i in basisfn_i:
                for prim_j in basisfn_j:
                    ij_kinenergy += prim_i.d * prim_j.d * \
                        prim_i.N * prim_j.N * \
                        kinenergy_primitive(prim_i.a, prim_j.a,
                                            prim_i.alpha, prim_j.alpha,
                                            prim_i.A, prim_j.A)
            kinenergy[i, j] = ij_kinenergy
            if i != j:
                kinenergy[j, i] = ij_kinenergy
            else:
                pass
    return kinenergy


def VRR(m, a, b, alpha, beta, A, B, C):
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
    P = (alpha * A + beta * B) / p
    K_AB = np.exp(-alpha * beta / p * np.linalg.norm(A - B) ** 2)
    T = p * np.linalg.norm(P - C) ** 2

    one = np.eye(len(a))
    two = 2 * one

    if np.all(a == 0) and np.all(b == 0):
        prim_int = (2 * np.pi / p) * K_AB * boys_function(m, T)
    elif np.any(a > 0):
        w = np.argmax(a > 0)

        prim_int = (P[w] - A[w]) * VRR(m, a - one[w, :], b,
                                       alpha, beta, A, B, C) \
            + (C[w] - P[w]) * VRR(m + 1, a - one[w, :], b,
                                  alpha, beta, A, B, C)

        if a[w] > 1:
            prim_int += (a[w] - 1) / (2 * p) *\
                (VRR(m, a - two[w, :], b, alpha, beta, A, B, C)
                 - VRR(m + 1, a - two[w, :], b, alpha, beta, A, B, C))
        else:
            pass

        if b[w] > 0:
            prim_int += (b[w] / (2 * p)) * \
                (VRR(m, a - one[w, :], b - one[w, :], alpha, beta, A, B, C)
                 - VRR(m + 1, a - one[w, :], b - one[w, :],
                       alpha, beta, A, B, C))
        else:
            pass

    elif np.any(b > 0):
        w = np.argmax(b > 0)

        prim_int = (P[w] - B[w]) * VRR(m, a, b - one[w, :],
                                       alpha, beta, A, B, C) \
            + (C[w] - P[w]) * VRR(m + 1, a, b - one[w, :],
                                  alpha, beta, A, B, C)

        if b[w] > 1:
            prim_int += (b[w] - 1) / (2 * p) *\
                (VRR(m, a, b - two[w, :], alpha, beta, A, B, C)
                 - VRR(m + 1, a, b - two[w, :], alpha, beta, A, B, C))
        else:
            pass
    else:
        pass
    return prim_int


def int_attraction(atoms, xyz, basis):
    M = len(basis)
    Vne = np.zeros((M, M))

    for i, basisfn_i in enumerate(basis):
        # Only calculate diagonal and upper off-diagonal values
        # to reduce number of iterations necessary.
        for j, basisfn_j in enumerate(basis[i:], i):
            ij_Vne = 0
            for atom, C in zip(atoms, xyz):
                for prim_i in basisfn_i:
                    for prim_j in basisfn_j:
                        ij_Vne -= atom * prim_i.d * prim_j.d * \
                            prim_i.N * prim_j.N * \
                            VRR(0, prim_i.a, prim_j.a,
                                prim_i.alpha, prim_j.alpha,
                                prim_i.A, prim_j.A, C)
            Vne[i, j] = ij_Vne
            if i != j:
                Vne[j, i] = ij_Vne
            else:
                pass
    return Vne


def eri_primitive(a, b, c, d, alpha, beta, gamma, delta, A, B, C, D):
    if isinstance(a, np.ndarray):
        pass
    else:
        a = np.array(a)
    if isinstance(b, np.ndarray):
        pass
    else:
        b = np.array(b)
    if isinstance(c, np.ndarray):
        pass
    else:
        c = np.array(c)
    if isinstance(d, np.ndarray):
        pass
    else:
        d = np.array(d)
    if isinstance(A, np.ndarray):
        pass
    else:
        A = np.array(A)
    if isinstance(B, np.ndarray):
        pass
    else:
        B = np.array(B)
    if isinstance(C, np.ndarray):
        pass
    else:
        C = np.array(C)
    if isinstance(D, np.ndarray):
        pass
    else:
        D = np.array(D)

    p = alpha + beta
    q = gamma + delta

    P = (alpha * A + beta * B) / p
    Q = (gamma * C + delta * D) / q
    W = (p * P + q * Q) / (p + q)

    AB = A - B
    CD = C - D

    # Precalculate [00|00]^(m) auxiliary integrals.
    T = p * q / (p + q) * np.linalg.norm(P - Q) ** 2
    m = np.arange(0, np.sum(a) + np.sum(b) + np.sum(c) + np.sum(d) + 1)
    K_AB = np.exp(-alpha * beta / p * np.linalg.norm(A - B) ** 2)
    K_CD = np.exp(-gamma * delta / q * np.linalg.norm(C - D) ** 2)
    ssss_m = 2 * np.pi ** (5 / 2) / (p * q * np.sqrt(p + q)) \
        * K_AB * K_CD * boys_function(m, T)

    PA = P - A
    WP = W - P
    QC = Q - C
    WQ = W - Q

    def eri_coeffs_1D(a, b, c, d, w):
        if b >= 1:
            t1 = eri_coeffs_1D(a + 1, b - 1, c, d, w)
            t2 = eri_coeffs_1D(a, b - 1, c, d, w)
            result = t1 + AB[w] * np.append(t2, 0).reshape(1, -1)
        elif d >= 1:
            t1 = eri_coeffs_1D(a, b, c + 1, d - 1, w)
            t2 = eri_coeffs_1D(a, b, c, d - 1, w)
            result = t1 + CD[w] * np.append(t2, 0).reshape(1, -1)
        elif a >= 1:
            t = eri_coeffs_1D(a - 1, 0, c, 0, w)
            result = PA[w] * np.append(t, 0).reshape(1, -1) + \
                WP[w] * np.append(0, t).reshape(1, -1)

            if a >= 2:
                t = eri_coeffs_1D(a - 2, 0, c, 0, w)
                v = (a - 1) / (2 * p) *\
                    (np.append(t, 0).reshape(1, -1) -
                     q / (p + q) * np.append(0, t).reshape(1, -1))
                result += np.append(v, 0).reshape(1, -1)
            else:
                pass
            if c >= 1:
                t = eri_coeffs_1D(a - 1, 0, c - 1, 0, w)
                v = c / (2 * (p + q)) * np.append(0, t).reshape(1, -1)
                result += np.append(v, 0).reshape(1, -1)
            else:
                pass
        elif c >= 1:
            t = eri_coeffs_1D(a, 0, c - 1, 0, w)
            result = QC[w] * np.append(t, 0).reshape(1, -1) + \
                WQ[w] * np.append(0, t).reshape(1, -1)

            if c >= 2:
                t = eri_coeffs_1D(a, 0, c - 2, 0, w)
                v = (c - 1) / (2 * q) *\
                    (np.append(t, 0).reshape(1, -1) -
                     q / (p + q) * np.append(0, t).reshape(1, -1))
                result += np.append(v, 0).reshape(1, -1)
            else:
                pass
            if a >= 1:
                t = eri_coeffs_1D(a - 1, 0, c - 1, 0, w)
                v = a / (2 * (p + q)) * np.append(0, t).reshape(1, -1)
                result += np.append(v, 0).reshape(1, -1)
            else:
                pass
        else:
            result = np.array([[1.0]])

        return result

    Cx = eri_coeffs_1D(a[0], b[0], c[0], d[0], 0)
    Cy = eri_coeffs_1D(a[1], b[1], c[1], d[1], 1)
    Cz = eri_coeffs_1D(a[2], b[2], c[2], d[2], 2)
    coeffs = convolve2d(convolve2d(Cx, Cy), Cz)
    integral = np.dot(coeffs.reshape(1, -1), ssss_m.reshape(-1, 1))

    return integral


def int_repulsion(basis):
    M = len(basis)
    ERI = np.zeros((M, M, M, M))

    for mu, basisfn_mu in enumerate(basis):
        for nu, basisfn_nu in enumerate(basis[:mu + 1]):
            for kap, basisfn_kap in enumerate(basis[:mu + 1]):
                for lam, basisfn_lam in enumerate(basis[:kap + 1]):
                    _ERI = 0
                    for k, prim_k in enumerate(basisfn_mu):
                        for l, prim_l in enumerate(basisfn_nu):
                            for n, prim_n in enumerate(basisfn_kap):
                                for o, prim_o in enumerate(basisfn_lam):
                                    _ERI += prim_k.d * prim_k.N * \
                                        prim_l.d * prim_l.N * \
                                        prim_n.d * prim_n.N * \
                                        prim_o.d * prim_o.N * \
                                        eri_primitive(prim_k.a, prim_l.a,
                                                      prim_n.a, prim_o.a,
                                                      prim_k.alpha,
                                                      prim_l.alpha,
                                                      prim_n.alpha,
                                                      prim_o.alpha,
                                                      prim_k.A, prim_l.A,
                                                      prim_n.A, prim_o.A)
                    ERI[mu, nu, kap, lam] = _ERI
                    ERI[nu, mu, kap, lam] = _ERI
                    ERI[mu, nu, lam, kap] = _ERI
                    ERI[nu, mu, lam, kap] = _ERI

                    ERI[kap, lam, mu, nu] = _ERI
                    ERI[kap, lam, nu, mu] = _ERI
                    ERI[lam, kap, mu, nu] = _ERI
                    ERI[lam, kap, nu, mu] = _ERI

    return ERI


def eerepulsion(ERI, P):
    M = P.shape[0]
    J = np.zeros((M, M))
    K = np.zeros((M, M))

    for mu in range(M):
        for nu in range(mu, M):
            _J = 0
            _K = 0
            for kap in range(M):
                for lam in range(M):
                    _J += P[kap, lam] * ERI[mu, nu, lam, kap]
                    # print("_J = {}".format(_J))
                    _K += 0.5 * P[kap, lam] * ERI[mu, kap, lam, nu]
                    # print("_K = {}".format(_K))
            J[mu, nu] = _J
            K[mu, nu] = _K
            if mu != nu:
                J[nu, mu] = _J
                K[nu, mu] = _K
    return (J, K)


def nucnucrepulsion(atoms, xyz):
    if isinstance(xyz, np.ndarray):
        pass
    else:
        xyz = np.array(xyz)
    Vnn = 0
    for i, at1 in enumerate(atoms[:-1]):
        for j, at2 in enumerate(atoms[i + 1:], i + 1):
            Vnn += at1 * at2 / np.linalg.norm(xyz[i, :] - xyz[j, :])
    return Vnn
