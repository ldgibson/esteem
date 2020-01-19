import numpy as np

from basissets import *
from esteem.hartreefock import *
from esteem.core import SelfConsistentField
from esteem.tests.testdata import MAT

testdata = MAT['testdata']
# atoms = [8, 1, 1]
# xyz = [[0, 0, 0.1272], [0, 0.7581, -0.5086], [0, -0.7581, -0.5086]]
atoms = [1, 2]
xyz = [[0, 0, 0], [0, 0, 1.4632]]
charge = 1
settings = {}
settings['basisset'] = 'STO-3G'
settings['tol_energy'] = 1e-8
settings['tol_density'] = 1e-8

basisset = '6-31G'


scf = SelfConsistentField(atoms, xyz, charge, settings)
scf.restricted_hartree_fock()
print(scf.epsilon)
print(scf.Etot)
# print(basisset, atoms, xyz)
# print(testdata[1].basisset, testdata[1].atoms, testdata[1].xyz)

# basisdef = basisread(testdata[1].basisset)
# basis = buildbasis(testdata[1].atoms, testdata[1].xyz.tolist(), basisdef)
# S = int_overlap(basis)
# Sref = testdata[1].results.S
# S_err = np.max(np.abs(S - Sref))
# print(S_err)
# print(testdata[1].xyz)
# bf = basis[0]
# for b in bf:
    # print(b)

# S = int_overlap(basis)
# T = int_kinenergy(basis)
# print(T)
# Vne = int_attraction(atoms, xyz, basis)
# eri_prim = eri_primitive(bf.A
# ERI = int_repulsion(basis)
# for i in range(2):
    # for j in range(2):
        # print(ERI[:, :, i, j])
# print(ERI)
# print("ERI[1, 0, 1, 0] = {}".format(ERI[1, 0, 1, 0]))
# print("ERI[4, 6, 2, 4] = {}".format(ERI[4, 6, 2, 4]))
# print("ERI[1, 3, 9, 14] = {}".format(ERI[1, 3, 9, 14]))
