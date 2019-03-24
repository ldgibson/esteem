import os.path
import sys

import numpy as np
from scipy.io import loadmat

from ..basissets import basisread, buildbasis


test_path = os.path.realpath(__file__)
package_path = os.path.join(test_path, '..')

test_data = loadmat(os.path.join(test_path, 'motest.mat'))
print(test_data)


def test_basisfiles():
    basissets = ['STO-3G', '6-31G', '6-311G', 'cc-pVDZ']
    dir_path = os.path.realpath(__file__)
    check_dir_path = os.path.join(dir_path, '..', 'basissets')
    if os.path.exists(check_dir_path):
        for basis in basissets:
            basisfile = basis + '.basis'
            assert os.path.exists(os.path.join(check_dir_path, basisfile)), \
                "Basis set files missing for {}.".format(basis)
    else:
        raise Exception("'basissets/' directory does not exist.")

    return True


def test_basisread():
    basis = basisread('6-31G')
    b61exp = np.array([7.8683, 1.8813, 0.5422])
    assert len(basis) == 10, "Basis set length incorrect."
    assert len(basis[6]) == 3
    assert np.all(basis[6][1] - b61exp < 1e-4)
    return True


# def test_buildbasis():

