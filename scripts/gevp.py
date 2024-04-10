## generalized eigenvalue problem

import numpy as np
import scipy

def gevp(C: np.ndarray, t0 = 0):
    """ Returns eigenvalues and eigenvectors of the GEVP (Generalized EigenValue Problem) """
    assert C.shape[0] == C.shape[1]
    N, T_ext = C.shape[0], C.shape[2]
    lam = np.matrix(np.zeros(shape=(N,T_ext)))
    C0 = np.matrix(C[:,:,t0])
    print(np.linalg.eigvals(C0))
    # x0, U = np.linalg.eigh(C0)
    # x0 = np.round(x0, decimals=16)
    # D = np.diag(x0)
    # sqrt_D = np.matrix(np.sqrt(D))
    # U = np.matrix(U)
    # print(sqrt_D.shape)
    # print(U.shape)
    # Q = U*sqrt_D*(U.H)
    # Q_inv = np.matrix(np.linalg.inv(Q))
    for t in range(T_ext):
        Ct = np.matrix(C[:,:,t])
        E, v = scipy.linalg.eigh(Ct, C0)
        #M = Q_inv * Ct * Q_inv.H
        #E, v = np.linalg.eigvalsh(M)
        #E_sorted = np.sort(E)
        lam[:, t] = E #E_sorted
    ####
    return lam
####

