import numpy as np

from .testdata import MAT


def test_overlap_primitive():
    try:
        from ..hartreefock import overlap_primitive
    except(ImportError):
        raise ImportError("overlap_primitive function missing.")

    Sp_threshold = 1e-6

    overlapprimitivetest = MAT['overlapprimitivetest']

    A = overlapprimitivetest['A'][0]
    B = overlapprimitivetest['B'][0]
    alpha = overlapprimitivetest['alpha'][0, 0]
    beta = overlapprimitivetest['beta'][0, 0]
    spd = overlapprimitivetest['ab']

    spd_str = ['s', 'px', 'py', 'pz', 'dx2',
               'dxy', 'dxz', 'dy2', 'dyz', 'dzz']

    for i1 in range(10):
        for i2 in range(10):
            a = spd[i1, :]
            b = spd[i2, :]
            Sref = overlapprimitivetest['S'][i1, i2]
            S = overlap_primitive(a, b, alpha, beta, A, B)
            relerr = np.abs((S - Sref) / Sref)
            assert relerr < Sp_threshold,\
                "{}|{}  error {} --FAIL".format(spd_str[i1],
                                                spd_str[i2], relerr)
            
    return True


def test_int_overlap():
    try:
        from ..hartreefock import int_overlap
    except(ImportError):
        raise ImportError("int_overlap function missing.")
    from ..basissets import basisread, buildbasis

    S_threshold = 1e-7

    for test in MAT['testdata']:
        bdef = basisread(test.basisset)
        basis = buildbasis(test.atoms, test.xyz, bdef)

        S = int_overlap(basis)
        Sref = test.results.S
        S_err = np.max(np.abs(S - Sref))

        # if S_err > S_threshold:
        print("Error: {} for test ".format(S_err) +\
              "{} with molecule {}".format(test.testid, test.molecule))

        # assert S_err < S_threshold, "Error: {} for test ".format(S_err) +\
            # "{} with molecule {}".format(test.testid, test.molecule)
    raise Exception()
    return 
