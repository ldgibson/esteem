import os.path as op

from .buildtestdata import buildtestdata


global MAT = buildtestdata()


def test_basisfiles():
    basissets = ['STO-3G', '6-31G', '6-311G', 'cc-pVDZ']
    dir_path = op.dirname(op.realpath(__file__))
    check_dir_path = op.join(dir_path, '..', 'basissets')
    if op.exists(check_dir_path):
        for basis in basissets:
            basisfile = basis + '.basis'
            assert op.exists(op.join(check_dir_path, basisfile)), \
                "Basis set files missing for {}.".format(basis)
    else:
        raise Exception("'basissets/' directory does not exist.")

    return True


def test_basisread():
    from ..basissets import basisread

    basis = basisread('6-31G')
    b61exp = np.array([7.8683, 1.8813, 0.5422])
    assert len(basis) == 10, "Basis set length incorrect."
    assert len(basis[6]) == 3
    assert np.all(basis[6][1] - b61exp < 1e-4)
    return True


def test_buildbasis():
    from ..basissets import buildbasis

    for test in MAT['testdata']:
        basisdef = basisread(test.basisset)
        atoms = test.atoms
        xyz = test.xyz
        basis = buildbasis(atoms, xyz, basisdef)
        basisref = test.results.basis

        M = len(basis)
        Mref = len(basisref)

        # Verify that the number of basis functions is correct.
        if M != Mref:
            raise Exception("Number of basis functions for"
                            " {} is not correct: ".format(test.molecule)
                            "{} given, {} expected.".format(M, Mref))
        else:
            pass

        # Verify that all the attributes exist.

        # Verify that all the attributes have the correct type and size.

        # Verify that the ordering of the basis functions is correct.

        # Sort test and reference basis functions.

        # Verify that basis functions are correct.

    return True
