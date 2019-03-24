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

        # Verify that all the attributes have the correct type and size.

        # Verify that the ordering of the basis functions is correct.

        # Sort test and reference basis functions.

        # Verify that basis functions are correct.

    return True
