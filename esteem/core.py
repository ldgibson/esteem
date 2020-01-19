import numpy as np
from scipy.linalg import eigh

from .basissets import basisread, buildbasis
from .hartreefock import *


class SelfConsistentField:
    def __init__(self, atoms, xyz, totalcharge, settings):
        self.atoms = atoms
        self.xyz = xyz
        self.totalcharge = totalcharge
        self.settings = settings

        self.basis = self.construct_basis(settings['basisset'])
        self.M = len(self.basis)
        self.S = np.zeros((self.M, self.M))
        self.T = np.zeros((self.M, self.M))
        self.Vne = np.zeros((self.M, self.M))
        self.Vnn = None
        self.J = np.zeros((self.M, self.M))
        self.K = np.zeros((self.M, self.M))
        self.Vee = np.zeros((self.M, self.M))
        self.ERI = np.zeros((self.M, self.M, self.M, self.M))
        self.Exc = None
        self.Vxc = np.zeros((self.M, self.M))
        self.rhoInt = None
        self.epsilon = np.zeros(self.M)
        self.C = np.zeros((self.M, self.M))
        self.P = np.zeros((self.M, self.M))
        self.E0 = None
        self.Etot = None
        return

    def construct_basis(self, basisset):
        basissetdef = basisread(basisset)
        basis = buildbasis(self.atoms, self.xyz, basissetdef)
        return basis

    def restricted_hartree_fock(self):
        N = np.sum(self.atoms) - self.totalcharge
        self.S = np.matrix(int_overlap(self.basis))
        self.T = np.matrix(int_kinenergy(self.basis))
        self.Vne = np.matrix(int_attraction(self.atoms, self.xyz, self.basis))
        self.ERI = int_repulsion(self.basis)
        self.Vnn = nucnucrepulsion(self.atoms, self.xyz)

        F = self.T + self.Vne
        eps, C = eigh(F, self.S)
        # eps = np.array([np.real(x) for x in eps if np.isreal(x)])
        self.epsilon = np.sort(eps.reshape(1, -1), axis=1)
        idx = eps.reshape(1, -1).argsort(axis=1)
        self.C = np.matrix(C[:, idx])

        for i in range(self.M):
            C_N = self.C[:, i].H * np.matrix(self.S) * self.C[:, i]
            C_N = C_N[0, 0]
            if C_N != 1.0:
                self.C[:, i] = self.C[:, i] / np.sqrt(C_N)
            else:
                pass

        self.P = 2 * self.C[:, :N // 2] * self.C[:, :N // 2].H

        self.E0 = np.trace(self.P * F)
        self.Etot = self.E0 + self.Vnn

        # SCF Procedure ------------------------------------------------------
        prev_energy = self.Etot
        prev_density = self.P
        conv_energy = 1
        conv_density = 1

        counter = 0
        converged = conv_energy < self.settings['tol_energy'] and \
            conv_density < self.settings['tol_density']

        while not converged:
            self.J, self.K = eerepulsion(self.ERI, self.P)
            self.Vee = np.matrix(self.J) - np.matrix(self.K)
            F = self.T + self.Vne + self.Vee  #np.matrix(self.J) - np.matrix(self.K)
            eps, self.C = eigh(F, self.S)
            eps = np.array([np.real(x) for x in eps if np.isreal(x)])
            self.epsilon = np.sort(eps.reshape(1, -1), axis=1)
            idx = eps.reshape(1, -1).argsort(axis=1)
            self.C = np.matrix(self.C[:, idx])

            for i in range(self.M):
                C_N = self.C[:, i].H * np.matrix(self.S) * self.C[:, i]
                C_N = C_N[0, 0]
                if np.abs(C_N - 1.0) > 1e-4:
                    self.C[:, i] /= np.sqrt(C_N)
                else:
                    pass

            self.P = 2 * self.C[:, :N // 2] * self.C[:, :N // 2].H

            A = self.T + self.Vne + 0.5 * self.Vee  # (self.J - self.K)
            self.E0 = np.trace(self.P * A)
            # self.E0 = np.trace(np.matmul(self.P, A))
            self.Etot = self.E0 + self.Vnn

            conv_energy = np.abs(prev_energy - self.Etot)
            conv_density = np.max(np.abs(prev_density - self.P))
            prev_energy = self.Etot
            prev_density = self.P
            converged = conv_energy < self.settings['tol_energy'] and \
                conv_density < self.settings['tol_density']

            print("i = {},\tE = {},\tE_conv = {},\tP_conv = {}"
                  .format(counter, self.Etot, conv_energy, conv_density))
            counter += 1
        return
