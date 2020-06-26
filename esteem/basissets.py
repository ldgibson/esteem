import os

import numpy as np
from scipy.special import factorial2


class BasisFunctionPrimitive:
    def __init__(self, atom=None, atom_id=None,
                 A=None, a=None, alpha=None, d=None, N=None):
        self.atom = atom
        self.atom_id = atom_id
        self.A = A
        self.a = a
        self.alpha = alpha
        self.d = d
        self.N = N
        return

    def __str__(self):
        desc = "Atom: {}\n".format(self.atom) +\
            "Atom index: {}\n".format(self.atom_id) +\
            "A: {}\n".format(self.A) +\
            "a: {}\n".format(self.a) +\
            "alpha: {}\n".format(self.alpha) +\
            "d: {}\n".format(self.d) +\
            "N: {}\n".format(self.N)
        return desc


class BasisFunction(BasisFunctionPrimitive):
    def __init__(self, atom=None, atom_id=None,
                 A=None, a=None, alpha=None, d=None, N=None):
        super().__init__(atom, atom_id, A, a, alpha, d, N)
        self.n_primitives = None
        return

    def normalize(self):
        a = np.array(self.a)
        N = (2 / np.pi) ** (3 / 4) * 2 ** np.sum(a) * \
            self.alpha ** ((2 * np.sum(a) + 3) / 4) / \
            np.sqrt(np.prod(factorial2(2 * a - 1)))
        self.N = N
        self.n_primitives = len(self.alpha)
        return

    def __getitem__(self, index):
        primitive = BasisFunctionPrimitive()
        primitive.atom = self.atom
        primitive.atom_id = self.atom_id
        primitive.A = self.A
        primitive.a = self.a
        primitive.alpha = self.alpha[index]
        primitive.d = self.d[index]
        primitive.N = self.N[index]
        return primitive

    def __len__(self):
        return len(self.alpha)

    def __repr__(self):
        return "Basis function of {} primitives, ".format(self.n_primitives) +\
            "atom: {}, atom index: {}".format(self.atom, self.atom_id)


class BasisSet(list):
    """Container of all basis functions for each atom.

    Parameters
    ----------
    atoms : list of int
        List containing the atomic numbers of each atom.
    xyz : array-like of float
        XYZ coordinates of each atom in Bohr (a0). Shape is N x 3,
        where N is the number of atoms."""

    def __init__(self, atoms, xyz):
        super().__init__()
        self.atoms = atoms
        if isinstance(xyz, list):
            self.xyz = xyz.copy()
        elif isinstance(xyz, np.ndarray):
            self.xyz = xyz.tolist()
        else:
            self.xyz = list(xyz)
        return

    def extract_basis_functions(self, definition):
        shelltype_exponents = {'S': [[0, 0, 0]],
                               'P': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                               'D': [[2, 0, 0], [0, 2, 0], [0, 0, 2],
                                     [1, 1, 0], [1, 0, 1], [0, 1, 1]]}
        # Loop over all atoms
        for i, (atom, xyz) in enumerate(zip(self.atoms, self.xyz)):
            # Loop over all basis function definitions.
            for basis_def in definition[atom]:
                # Check if current definition is for 'SP'.
                if basis_def.shelltype == 'SP':
                    shells = ['S', 'P']
                else:
                    shells = [basis_def.shelltype]
                # Loop over each of the shells.
                # This loop will only occur twice iff the shelltype is 'SP'.
                for shell in shells:
                    # For any shell other than 'S', this will loop over
                    # each of the sets of cartesian exponents.
                    for cartesian_exponents in shelltype_exponents[shell]:
                        basisfunction = BasisFunction()
                        basisfunction.atom = atom
                        basisfunction.atom_id = i
                        basisfunction.A = xyz.copy()
                        basisfunction.a = cartesian_exponents
                        if basis_def.shelltype == 'SP':
                            if shell == 'S':
                                coeffs = basis_def.coefficients[0, :]
                            else:
                                coeffs = basis_def.coefficients[1, :]
                        else:
                            coeffs = basis_def.coefficients
                        basisfunction.alpha = basis_def.exponents
                        basisfunction.d = coeffs
                        basisfunction.normalize()
                        self.append(basisfunction)


class AtomicBasisSetDefinition:
    def __init__(self, shelltype=None, exponents=[], coefficients=[]):
        self.shelltype = shelltype
        self.exponents = exponents
        self.coefficients = coefficients
        return

    def __str__(self):
        print_str = 'Shell type ' +\
            '{} with {} primitives'.format(self.shelltype, len(self.exponents))
        return print_str


class BasisSetDefinition(dict):
    """Container for basis set definitions.

    Each key is an atomic number and the value is a list of atomic
    basis set definitions for each orbital type.

    Example
    -------
    >>> basisdef[1]  # List of orbital definitions for hydrogen.
    >>> basisdef[1][0]  # First primitive definition for hydrogen.
    """
    pass


def basisread(basis_set_string):
    elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                'Na', 'Mg', 'Al', 'S', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
                'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']
    basis_func_num_dict = {'STO-3G': [1, 2],
                           '6-31G': [2, 3],
                           '6-311G': [3, 4],
                           'cc-pVDZ': [5, 9]}

    if basis_set_string not in basis_func_num_dict.keys():
        raise Exception("Basis set not recognized.\n"
                        "Please choose from the following options:\n"
                        "\tSTO-3G, 6-31G, 6-311G, cc-pVDZ")
    else:
        pass

    num_basis_func_arr = basis_func_num_dict[basis_set_string]
    basissetdef = BasisSetDefinition()
    path_to_files = os.path.dirname(os.path.realpath(__file__))
    basis_set_file_name = os.path.join(path_to_files, 'basissets',
                                       basis_set_string + '.basis')
    with open(basis_set_file_name, 'r') as infile:
        atomic_number = 1
        for line in infile:
            # Skip comment lines
            if line.startswith('!'):
                continue
            # Skip separators
            elif line.startswith('****'):
                continue
            # Skip blank lines.
            elif line == '\n':
                continue
            else:
                pass

            line = line.split()
            # Record all basis function parameters for an element.
            if line[0] in elements:
                if line[0] == 'H' or line[0] == 'He':
                    num_basis_func = num_basis_func_arr[0]
                else:
                    num_basis_func = num_basis_func_arr[1]

                # List containing all of the primitive definitions.
                all_atom_defs = []
                # Repeat for number of basis functions.
                for i in range(num_basis_func):
                    atomic_basis = AtomicBasisSetDefinition()
                    info = next(infile)
                    info = info.split()
                    shell_type = info[0]
                    num_primitives = int(info[1])
                    # prefactor = info[2]

                    exponents = np.zeros(num_primitives)
                    if shell_type == 'SP':
                        contraction_coeffs = np.zeros((2, num_primitives))
                    else:
                        contraction_coeffs = np.zeros(num_primitives)

                    for j in range(num_primitives):
                        newline = next(infile)
                        exp_and_coeff = [float(n) for n in
                                         newline.strip().split()]
                        exponents[j] = exp_and_coeff[0]
                        if shell_type == 'SP':
                            contraction_coeffs[:, j] = exp_and_coeff[1:]
                        else:
                            contraction_coeffs[j] = exp_and_coeff[1]
                    atomic_basis.shelltype = shell_type
                    atomic_basis.exponents = exponents
                    atomic_basis.coefficients = contraction_coeffs
                    all_atom_defs.append(atomic_basis)

                basissetdef[atomic_number] = all_atom_defs
                # Increment atomic number
                atomic_number += 1
    return basissetdef


def buildbasis(atoms, xyz_a0, basissetdef):
    basis = BasisSet(atoms, xyz_a0)
    basis.extract_basis_functions(basissetdef)
    return basis
