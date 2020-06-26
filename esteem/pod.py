import numpy as np


def lowdin_orthogonalization(S):
    """Lowdin symmetrical orthogonalization for basis
    vectors in 3D Cartesian space."""
    
    # Get overlap matrix if not given
    # S = get_overlap(X, basis_dim)
    
    # Diagonalize overlap matrix
    eigs, U = np.linalg.eig(S)
    
    # Wwitch vector to diagonal matrix
    S_diag = np.matrix(np.diagflat(eigs))
    
    # Get sqrt of S_diag
    S_diag_sqrt = np.sqrt(S_diag)
    
    # Get S_inv_sqrt
    S_sqrt = U * S_diag_sqrt * U.H  # U.H is adjoint of U
    S_inv_sqrt = np.linalg.inv(S_sqrt)
    
    # Transform the basis
    # X_ = S_inv_sqrt * X
    
    return S_inv_sqrt


def get_overlap(X, basis_dim):
    S = np.matrix(np.zeros((basis_dim, basis_dim)))
    
    for i in range(basis_dim):
        for j in range(i+1, basis_dim):
            S[i, j] = np.dot(to_arr(X[i]),to_arr(X[j]))
            S[j, i] = S[i, j]
    np.fill_diagonal(S, 1.0)
    return S


def to_arr(x):
    return np.squeeze(np.asarray(x))
