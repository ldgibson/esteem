import os.path

import numpy as np

from .testdata import MAT

test_path = os.path.dirname(os.path.realpath(__file__))
package_path = os.path.join(test_path, '..')


def test_basisfiles():
    basissets = ['STO-3G', '6-31G', '6-311G', 'cc-pVDZ']
    check_dir_path = os.path.join(package_path, 'basissets')
    if os.path.exists(check_dir_path):
        for basis in basissets:
            basisfile = basis + '.basis'
            assert os.path.exists(os.path.join(check_dir_path, basisfile)), \
                "Basis set files missing for {}.".format(basis)
    else:
        raise FileNotFoundError("'basissets/' directory does not exist.")

    return True


def test_basisread():
    try:
        from ..basissets import basisread
    except ModuleNotFoundError:
        raise ModuleNotFoundError("'basissets' module does not exist.")
    except ImportError:
        raise ImportError("'basisread' does not exist.")

    basisdef = basisread('6-31G')
    b61exp = np.array([7.8683, 1.8813, 0.5442])
    assert len(basisdef) == 10, "Basis set length incorrect."
    assert len(basisdef[6]) == 3
    assert np.all(np.abs(basisdef[6][1].exponents - b61exp) < 1e-4)
    return True


def test_buildbasis():
    from ..basissets import basisread
    try:
        from ..basissets import buildbasis
    except ModuleNotFoundError:
        raise ModuleNotFoundError("'basissets' module does not exist.")
    except ImportError:
        raise ImportError("'buildbasis' does not exist.")

    for test in MAT['testdata']:
        basisdef = basisread(test.basisset)
        basis = buildbasis(test.atoms, test.xyz, basisdef)
        basisref = test.results.basis

        M = len(basis)
        Mref = len(basisref)

        # Verify that the number of basis functions is correct.
        if M != Mref:
            raise Exception("Number of basis functions for " +
                            "{} is not correct: ".format(test.molecule) +
                            "{} given, {} expected.".format(M, Mref))
        else:
            pass

        # Verify that all the attributes exist.
        attributes = ['atom', 'A', 'a', 'alpha', 'd', 'N']
        for attr in attributes:
            for basisfn in basis:
                assert attr in vars(basisfn).keys(),\
                    "'{}' missing.".format(attr)
                assert np.all(vars(basisfn)[attr] != None),\
                    "'{}' value does not exist.".format(attr)

        # Verify that all the attributes have the correct type and size.
        for p, (basisfn, refbasisfn) in enumerate(zip(basis, basisref)):
            for attr in attributes:
                pass
                # if not isinstance(vars(basisfn)[attr],
                                  # type(vars(refbasisfn)[attr])):
                    # raise Exception("Wrong dtype for '{}' ".format(attr) +\
                                    # "field in basis function {}.".format(p))
                # else:
                    # pass
                # if len(vars(basisfn)[attr]) != len(vars(refbasisfn)[attr]):
                    # raise Exception("Wrong size for '{}' ".format(attr) +\
                                    # "field in basis function {}.".format(p))
                # else:
                    # pass

        # Verify that the ordering of the basis functions is correct.

        # Sort test and reference basis functions.

        # Verify that basis functions are correct.
        for p, (basisfn, refbasisfn) in enumerate(zip(basis, basisref)):
            assert basisfn.atom == refbasisfn.atom,\
                "Basis function {} is centered on ".format(p) +\
                "the wrong atom type."
            assert np.all(np.abs(basisfn.A - refbasisfn.A) < 1e-3),\
                "Basis function {} is centered on ".format(p) +\
                "the wrong atom position."
            assert np.all(basisfn.a == refbasisfn.a),\
                "Basis function {} has ".format(p) +\
                "incorrect cartesian exponents."
            assert np.all(np.abs(basisfn.alpha - refbasisfn.alpha)
                / refbasisfn.alpha < 1e-3),\
                "Basis function {} has ".format(p) +\
                "incorrect cartesian exponents."
            assert np.all(np.abs(basisfn.d - refbasisfn.d) < 1e-3),\
                "Basis function {} has ".format(p) +\
                "incorrect contraction coefficients."
            assert np.all(np.abs(basisfn.N - refbasisfn.N)
                / refbasisfn.N < 1e-5),\
                "Basis function {} has ".format(p) +\
                "incorrect incorrect primitive normalization constants."
    return True
