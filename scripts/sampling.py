## resampling techniques from gauge configurations

import numpy as np
from scipy.stats import bootstrap

from . import uwerr

def uncorrelated_confs_to_bts(x, N_bts, seed=12345):
    """Bootstrap samples from array of data
    
    - generates N_b = N_bts/block_size averages
      y[i] = \sum_{j=1}^{block_size} x[i*block_size+j] / N_b
    - samples N_bts values z[i] (with repetition) from the y[i] and returns z
    
    Args:
        x (np.ndarray): time series. Bootstrapping is done on 1st index
        N_bts (int): Number of bootstraps
        block_size (int): size of the block
    
    Returns:
        np.ndarray: Bootstrap samples
    """
    np.random.seed(seed=seed) # rng = np.random.default_rng(seed=seed)
    N = x.shape[0]
    return np.array([np.average(x[np.random.randint(0,high=N,size=N_bts,dtype=int)], axis=0) for i in range(N_bts)])
####

def correlated_confs_to_bts(Cg: np.ndarray, N_bts: int, block_size=2, seed=12345, output_file=None) -> np.ndarray:
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
    Cg_uncorr = Cg[0:Ng:tauint] ## uncorrelated values
    return uncorrelated_confs_to_bts(x=Cg_uncorr, N_bts=N_bts, seed=seed)
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
    C_bts = uncorrelated_confs_to_bts(Cg, N_bts=N_bts)
    print(np.average(Cg))
    print(np.average(C_bts))
####
