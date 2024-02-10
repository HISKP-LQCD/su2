## resampling techniques from gauge configurations

import numpy as np
from . import uwerr

def correlated_confs_to_bts(Cg: np.ndarray, N_bts: int, seed=12345, output_file=None) -> np.ndarray:
    """Bootstrap samples from array of correlated configurations

    - The configurations are sampled every tau_int, (integrated autocorrelation time)
    - The bootstrap samples are drawn from the uncorrelated configurations

    Args:
        C (np.ndarray): "correlator" (observable computed for each configuration)
        N_bts (int): Number of bootstraps

    Returns:
        np.ndarray: Bootstrap samples
    """
    Ng = Cg.shape[0] ## total number of configurations
    tauint = int(uwerr.uwerr_primary(Cg, output_file=output_file)["tauint"]) ## integrated autocorrelation time
    if tauint == 0:
        tauint = 1 ## uncorrelated data
    ####
    Ng_uncorr = int(Ng/tauint) ## number of uncorrelated configurations
    Cg_uncorr = Cg[0:Ng:tauint,] #np.array([Cg[i] for i in range(0, Ng, tauint)]) ## uncorrelated values of the observabe
    np.random.seed(seed=seed) ## setting the seed
    bts_idx = np.random.randint(0, high=Ng_uncorr, size=N_bts, dtype=int) ## indices of bootstrap samples
    C_bts = Cg_uncorr[bts_idx,] #np.array([ for i in bts_idx])
    return C_bts
####

def uncorrelated_confs_to_bts(Cg: np.ndarray, bts_idx) -> np.ndarray:
    """Bootstrap samples from array of uncorrelated configurations

    Args:
        C (np.ndarray): "correlator" (observable computed for each configuration)
        bts_idx : list of indices of the bootstrap samples

    Returns:
        np.ndarray: Bootstrap samples
    """
    return np.array([Cg[i] for i in bts_idx])
####


if __name__ == "__main__":
    Cg = np.array(100*[list(np.random.normal(0.0, 0.4, 50))]).flatten()
    seed = 12345
    N_bts = 1000
    C_bts = correlated_confs_to_bts(Cg, N_bts=N_bts, seed=seed)

    print(np.average(Cg))
    print(np.average(C_bts))

    np.random.seed(seed=seed)
    Ng = Cg.shape[0]
    bts_idx = np.random.randint(0, high=Ng, size=N_bts, dtype=int)
    C_bts = uncorrelated_confs_to_bts(Cg, bts_idx=bts_idx)
    print(np.average(Cg))
    print(np.average(C_bts))
####
