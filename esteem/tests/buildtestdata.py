import os.path as op

import numpy as np
from scipy.io import loadmat

from ..basissets import BasisSet, BasisFunction


class TestCase:
    def __init__(self, testid=None, molecule=None, atoms=None,
                 xyz=None, charge=None, basisset=None, method=None,
                 tol_density=None, tol_energy=None,
                 exchange_functional=None, correlation_functional=None,
                 n_radial_points=None, n_angular_points=None,
                 results=None):
        self.testid = testid
        self.molecule = molecule
        self.atoms = atoms
        self.xyz = xyz
        self.charge = charge
        self.basisset = basisset
        self.method = method
        self.tol_density = tol_density
        self.tol_energy = tol_energy
        self.exchange_functional = exchange_functional
        self.correlation_functional = correlation_functional
        self.n_radial_points = n_radial_points
        self.n_angular_points = n_angular_points
        self.results = results
        return


class TestResults:
    def __init__(self, atoms, xyz, basis=None, S=None, T=None, Vne=None,
                 Vee=None, Vxc=None, Exc=None, ERI=None, C=None, P=None,
                 epsilon=None, E0=None, Etot=None):
        self.atoms = atoms
        self.xyz = xyz
        self.basis = basis
        self.S = S
        self.T = T
        self.Vne = Vne
        self.Vee = Vee
        self.Vxc = Vxc
        self.Exc = Exc
        self.ERI = ERI
        self.C = C
        self.P = P
        self.epsilon = epsilon
        self.E0 = E0
        self.Etot = Etot
        return

    def convert(self, results):
        self.basis = self.convert_to_basis_set(results[0].flatten())
        self.S = np.matrix(results[1])
        self.T = np.matrix(results[2])
        self.Vne = np.matrix(results[3])
        self.Vee = np.matrix(results[4])
        self.Vxc = np.matrix(results[5])
        self.Exc = results[6][0, 0]
        self.ERI = results[7]
        self.C = np.matrix(results[8])
        self.P = np.matrix(results[9])
        self.epsilon = results[10].flatten()
        self.E0 = results[11][0, 0]
        self.Etot = results[12][0, 0]
        return

    def convert_to_basis_set(self, basisarray):
        basis = BasisSet(self.atoms, self.xyz)
        for bf in basisarray:
            basis_func = BasisFunction()
            basis_func.atom = int(bf[0][0, 0])
            basis_func.A = bf[1].reshape(3)
            basis_func.a = bf[2].reshape(3)
            basis_func.alpha = bf[3].flatten()
            basis_func.d = bf[4].flatten()
            basis_func.N = bf[5].flatten()
            basis.append(basis_func)
        return basis


def buildtestdata():
    npmat = loadmat(op.join(op.dirname(op.realpath(__file__)), 'motest.mat'))
    optest = npmat['overlapprimitivetest'][0][0]
    mat = {}
    mat['overlapprimitivetest'] = {}
    fields = ['A', 'B', 'alpha', 'beta', 'ab', 'S']
    for i, field in enumerate(fields):
        mat['overlapprimitivetest'][field] = optest[i]

    mat['testdata'] = []
    for i, t in enumerate(npmat['testdata'].flatten()):
        test = TestCase()
        test.testid = i
        test.molecule = t[1][0]
        test.atoms = t[2].flatten()
        test.xyz = np.matrix(t[3])
        test.charge = t[4][0, 0]
        test.basisset = t[5][0]
        test.method = t[6][0]
        test.tol_density = t[7][0, 0]
        test.tol_energy = t[8][0, 0]

        if test.method == 'RKS':
            test.exchange_functional = t[9][0]
            test.correlation_functional = t[10][0]
            test.n_radial_points = t[11][0, 0]
            test.n_angular_points = t[12][0, 0]
        else:
            pass

        test.results = TestResults(test.atoms, test.xyz)
        test.results.convert(t[13][0, 0])
        mat['testdata'].append(test)
    return mat
