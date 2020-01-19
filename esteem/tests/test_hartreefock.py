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
        # print("Error: {} for test ".format(S_err) +\
              # "{} with molecule {}".format(test.testid, test.molecule))

        assert S_err < S_threshold, "Error: {} for test ".format(S_err) +\
            "{} with molecule {}".format(test.testid, test.molecule)
    # raise Exception()
    return 


def test_int_kinenergy():
    try:
        from ..hartreefock import int_kinenergy
    except(ImportError):
        raise ImportError("int_kinenergy function missing.")
    from ..basissets import basisread, buildbasis

    T_threshold = 1e-7

    for test in MAT['testdata']:
        bdef = basisread(test.basisset)
        basis = buildbasis(test.atoms, test.xyz, bdef)

        T = int_kinenergy(basis)
        Tref = test.results.T
        T_err = np.max(np.abs(T - Tref))

        # if S_err > S_threshold:
        # print("Error: {} for test ".format(S_err) +\
              # "{} with molecule {}".format(test.testid, test.molecule))

        assert T_err < T_threshold, "Error: {} for test ".format(T_err) +\
            "{} with molecule {}".format(test.testid, test.molecule)
    # raise Exception()
    return 


def test_int_attraction():
    try:
        from ..hartreefock import int_attraction
    except(ImportError):
        raise ImportError("int_kinenergy function missing.")
    from ..basissets import basisread, buildbasis
    
    Vne_threshold = 1e-7

    for test in MAT['testdata']:
        bdef = basisread(test.basisset)
        basis = buildbasis(test.atoms, test.xyz, bdef)

        Vne = int_attraction(test.atoms, test.xyz, basis)
        Vne_ref = test.results.Vne
        Vne_err = np.max(np.abs(Vne - Vne_ref))

        # if S_err > S_threshold:
        # print("Error: {} for test ".format(S_err) +\
              # "{} with molecule {}".format(test.testid, test.molecule))

        assert Vne_err < Vne_threshold, "Error: {} for test".format(Vne_err) +\
            " {} with molecule {}".format(test.testid, test.molecule)
    # raise Exception()
    return 


def test_int_repulsion():
    try:
        from ..hartreefock import int_repulsion
    except(ImportError):
        raise ImportError("int_repulsion function missing.")
    from ..basissets import basisread, buildbasis
    
    ERI_threshold = 1e-7

    for test in MAT['testdata']:
        bdef = basisread(test.basisset)
        basis = buildbasis(test.atoms, test.xyz, bdef)

        ERI = int_repulsion(basis)
        ERI_ref = test.results.ERI
        ERI_err = np.max(np.abs(ERI - ERI_ref))

        # if S_err > S_threshold:
        print("Error: {} for test ".format(ERI_err) +\
              "{} with molecule {}".format(test.testid, test.molecule))

        # assert ERI_err < ERI_threshold, "Error: {} for test".format(ERI_err) +\
            # " {} with molecule {}".format(test.testid, test.molecule)
    raise Exception()
    return 


def test_eerepulsion():
    try:
        from ..hartreefock import eerepulsion
    except(ImportError):
        raise ImportError("eerepulsion function missing.")
    from ..basissets import basisread, buildbasis
    
    Vee_threshold = 1e-6

    for test in MAT['testdata']:
        bdef = basisread(test.basisset)
        basis = buildbasis(test.atoms, test.xyz, bdef)

        ERI = int_repulsion(basis)
        ERI_ref = test.results.ERI
        ERI_err = np.max(np.abs(ERI - ERI_ref))

        # if S_err > S_threshold:
        print("Error: {} for test ".format(ERI_err) +\
              "{} with molecule {}".format(test.testid, test.molecule))

        # assert ERI_err < ERI_threshold, "Error: {} for test".format(ERI_err) +\
            # " {} with molecule {}".format(test.testid, test.molecule)
    raise Exception()
    return 
