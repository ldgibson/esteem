import numpy as np
from scipy.linalg import eigh

from .basissets import basisread, buildbasis
from .hartreefock import *
from .pod import lowdin_orthogonalization


class SelfConsistentField:
    def __init__(self, atoms, xyz, totalcharge, settings, multiplicity=1):
        self.atoms = atoms
        self.xyz = xyz
        self.totalcharge = totalcharge
        self.settings = settings
        self.multiplicity = multiplicity

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
        self.F = None
        self.Fa = None
        self.Fb = None
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

            self.F = F

            self.E0 = np.trace(self.P * A)
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


    def unrestricted_hartree_fock(self):
        N = np.sum(self.atoms) - self.totalcharge

        self.S = np.matrix(int_overlap(self.basis))
        self.T = np.matrix(int_kinenergy(self.basis))
        self.Vne = np.matrix(int_attraction(self.atoms, self.xyz, self.basis))
        self.ERI = int_repulsion(self.basis)
        self.Vnn = nucnucrepulsion(self.atoms, self.xyz)

        N_unpaired = self.multiplicity - 1

        if self.multiplicity > 1:
            Nb = (N - N_unpaired)//2
            Na = Nb + N_unpaired
        else:
            Nb = N // 2
            Na = Nb

        F = self.T + self.Vne
        eps, C = eigh(F, self.S)
        self.epsilon_a = np.sort(eps.reshape(1, -1), axis=1)
        self.epsilon_b = self.epsilon_a
        idx = eps.reshape(1, -1).argsort(axis=1)
        self.Ca = np.matrix(C[:, idx])
        self.Cb = self.Ca.copy()

        for i in range(self.M):
            C_N_a = self.Ca[:, i].H * np.matrix(self.S) * self.Ca[:, i]
            C_N_a = C_N_a[0, 0]
            C_N_b = self.Cb[:, i].H * np.matrix(self.S) * self.Cb[:, i]
            C_N_b = C_N_b[0, 0]
            if C_N_a != 1.0:
                self.Ca[:, i] = self.Ca[:, i] / np.sqrt(C_N_a)
            else:
                pass
            if C_N_b != 1.0:
                self.Cb[:, i] = self.Cb[:, i] / np.sqrt(C_N_b)
            else:
                pass

        self.Pa = self.Ca[:, :Na] * self.Ca[:, :Na].H
        self.Pb = self.Cb[:, :Nb] * self.Cb[:, :Nb].H

        Aa = (1 / 2) * (self.T + self.Vne)
        Ab = (1 / 2) * (self.T + self.Vne)

        self.E0 = np.trace(self.Pa * Aa) + np.trace(self.Pb * Ab)
        self.Etot = self.E0 + self.Vnn

        # SCF Procedure ------------------------------------------------------
        prev_energy = self.Etot
        prev_density_alpha = self.Pa
        prev_density_beta = self.Pb
        conv_energy = 1
        conv_density_alpha = 1
        conv_density_beta = 1

        counter = 0
        converged = conv_energy < self.settings['tol_energy'] and \
            conv_density_alpha < self.settings['tol_density'] and \
            conv_density_beta < self.settings['tol_density']

        while not converged:
            self.Ja, self.Ka = eerepulsion(self.ERI, self.Pa)
            self.Jb, self.Kb = eerepulsion(self.ERI, self.Pb)

            self.Ka = 2 * self.Ka
            self.Kb = 2 * self.Kb

            self.J = self.Ja + self.Jb

            Fa = self.T + self.Vne + np.matrix(self.J) - np.matrix(self.Ka)
            Fb = self.T + self.Vne + np.matrix(self.J) - np.matrix(self.Kb)

            eps_a, self.Ca = eigh(Fa, self.S)
            eps_b, self.Cb = eigh(Fb, self.S)
            eps_a = np.array([np.real(x) for x in eps_a if np.isreal(x)])
            eps_b = np.array([np.real(x) for x in eps_b if np.isreal(x)])
            self.epsilon_a = np.sort(eps_a.reshape(1, -1), axis=1)
            self.epsilon_b = np.sort(eps_b.reshape(1, -1), axis=1)
            idx_a = eps_a.reshape(1, -1).argsort(axis=1)
            idx_b = eps_b.reshape(1, -1).argsort(axis=1)
            self.Ca = np.matrix(self.Ca[:, idx])
            self.Cb = np.matrix(self.Cb[:, idx])

            for i in range(self.M):
                C_N_a = self.Ca[:, i].H * np.matrix(self.S) * self.Ca[:, i]
                C_N_a = C_N_a[0, 0]
                C_N_b = self.Cb[:, i].H * np.matrix(self.S) * self.Cb[:, i]
                C_N_b = C_N_b[0, 0]
                if np.abs(C_N_a - 1.0) > 1e-4:
                    self.Ca[:, i] = self.Ca[:, i] / np.sqrt(C_N_a)
                else:
                    pass
                if np.abs(C_N_b - 1.0) > 1e-4:
                    self.Cb[:, i] = self.Cb[:, i] / np.sqrt(C_N_b)
                else:
                    pass

            self.Pa = self.Ca[:, :Na] * self.Ca[:, :Na].H
            self.Pb = self.Cb[:, :Nb] * self.Cb[:, :Nb].H

            Aa = (1 / 2) * (self.T + self.Vne + Fa)
            Ab = (1 / 2) * (self.T + self.Vne + Fb)

            self.Fa = Fa
            self.Fb = Fb

            self.E0 = np.trace(self.Pa * Aa) + np.trace(self.Pb * Ab)
            self.Etot = self.E0 + self.Vnn

            conv_energy = np.abs(prev_energy - self.Etot)
            conv_density_alpha = np.max(np.abs(prev_density_alpha - self.Pa))
            conv_density_beta = np.max(np.abs(prev_density_beta - self.Pb))
            prev_energy = self.Etot
            prev_density_alpha = self.Pa
            prev_density_beta = self.Pb
            converged = conv_energy < self.settings['tol_energy'] and \
                conv_density_alpha < self.settings['tol_density'] and \
                conv_density_beta < self.settings['tol_density']

            print("i = {},\tE = {},\tE_conv = {},\tP_conv_alpha = {},\tP_conv_beta = {}"
                  .format(counter, self.Etot, conv_energy, conv_density_alpha, conv_density_beta))
            counter += 1
        return
