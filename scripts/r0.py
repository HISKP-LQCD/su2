## find r0 from the Wilson loops

import os
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

import effective_curves
import sampling

def read_Wloops_1st_format(omeas_dir, i1, step, Ng, nmax_iter, L, T):
    """Reading the Wilson loops in the 1st type of format 

    Args:
        omeas_dir (str): directory storing the online measurements
        i1 (int): index of 1st configuration
        step (int): stap between indices of consedutive configurations
        Ng (int): number of configurations we want to extract
        nmax_iter (int): maximum number of attempts to read configurations (avoids getting stuck when files don't exist)
        L (int): Spatial extent of the lattice
        T (int): Temporal extent of the lattice

    Returns:
        numpy.ndarray: array with 3 indices: configuration, t (temporal extent of the loop), r (spatial extent of the loop)
    """
    data = np.zeros(shape=(Ng, T-1, L-1)) # t=0 or r=0 are not valid loops
    ig = 0
    for it in range(nmax_iter):
        path = "{d}/wilsonloop.{i:06}.dat".format(d=omeas_dir, i=i1 + step*it)
        if os.path.exists(path):
            df = pd.read_csv(path, sep=" ").drop(["t"], axis=1)
            data[ig,:] = df.to_numpy()
            ig += 1
        ####
        if ig == Ng:
            return data ## Ng configurations have been read
        ####
    ####
    err_msg = "Cannot read {Ng} configurations from {omeas_dir} starting from {i1} with step {step}".format(Ng=Ng, omeas_dir=omeas_dir, i1=i1, step=step)
    raise RuntimeError(err_msg)
####

def read_Wloops_format2(omeas_dir: str, L: int, T: int, beta: float, xi: float, ss: bool, sizeWloops: float, nape=0, alpha=1.0) -> np.ndarray:
    """Reading the Wilson loops in the 2nd type of format 

    Args:
        omeas_dir (str): directory storing the online measurements
        L (int): Spatial extent of the lattice
        T (int): Temporal extent of the lattice
        beta (float): beta in the Lagrangian
        xi (float): bare anisotropy
        ss (bool): spatial-spatial Wilson loops
        sizeWloops (float): ratio between size and maximum extent of measured loops
        nape (int, optional): number of APE smearing steps. Defaults to 0.
        alpha (float, optional): APE smearing parameter. Defaults to 1.0.

    Returns:
        numpy.ndarray: array with 3 indices: configuration, t (temporal extent of the loop), r (spatial extent of the loop)
    """
    ss_name = "coarse" if ss else "fine"
    path = "{d}/result2p1d.u1potential.rotated.Nt{T}.Ns{L}.b{beta:.6f}.xi{xi:.6f}.nape{nape}.alpha{alpha:.6f}{ss_name}distance".format(
        d=omeas_dir, T=T, L=L, beta=beta, xi=xi, nape=nape, alpha=alpha, ss_name=ss_name)
    L_ext = int(sizeWloops*L)
    T_ext = sizeWloops*L if ss else sizeWloops*T
    T_ext = int(T_ext)
    df = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None, usecols=[i for i in range(L_ext*T_ext)]).to_numpy()
    Ng = df.shape[0]
    data = np.zeros(shape=(Ng, T_ext, L_ext)) # t=0 or r=0 are not valid loops
    for t in range(T_ext):
        for r in range(L_ext):
            data[:, t, r] =  df[:, t*L_ext + r]
        ####
    ####
    return data
####


def get_potential_bts(W_bts: np.ndarray, t1:int, t2: int) -> np.ndarray:
    """Extracting the potential from the Wilson loop by fitting between t1 and t2 included

    Args:
        W_bts (np.ndarray): bootstraps for W(t) \sim e^{-V*t}. shape = (N_bts, T)
        t1 (int): start of the plateau of V_eff(t)
        t2 (int): end of the plateau of V_eff(t)

    Return:
        V_bts (np.ndarray): bootstraps for V
    """
    N_bts = W_bts.shape[0]
    V_eff = np.array([effective_curves.get_m_eff(W_bts[i,t1:(t2+1)], strategy="log") for i in range(N_bts)])
    dV_eff = np.std(a=V_eff, axis=0, ddof=1)
    V_bts = np.array([effective_curves.fit_eff_mass(m_eff=V_eff[i,:], dm_eff=dV_eff) for i in range(N_bts)])
    return V_bts
####

def fit_static_potential_QED2p1(V_bts_all: np.ndarray, r1: int, r2: int, p0=[1.0,1.0,1.0]) -> np.ndarray:
    """QED in 2+1 dimensions: fit the static potential V(r):

    V(r) = V_0 + b*log(r) + sigma*r,

    where V_0 is an unphysical offset constant

    Args:
        V_bts_all (np.ndarray): bootstraps for the potential. shape = (N_bts, L)
        r1 (int): start of the fit interval for V(r)
        r2 (int): end of the fit interval for F(r)
        p0 (list, optional): initial fit guess for "v0", "b" and "sigma" . Defaults to [1.0, 1.0, 1.0].

    Returns:
        np.ndarray: _description_
    """
    N_bts, L = V_bts_all.shape
    V_bts = V_bts_all[:,r1:(r2+1)]
    if r1 < 0 or r2 > L-1:
        raise ValueError("Invalid fit range: r1={r1}, r2={r2}, L={L}".format(r1=r1, r2=r2, L=L))
    ####
    x = [i for i in range(r1, r2+1)]
    dy = np.std(V_bts, axis=0, ddof=1)
    ansatz = lambda r, v0, b, sigma: v0 + b*np.log(r) + sigma*r
    par = np.array([curve_fit(ansatz, x, V_bts[i,:], sigma=dy, p0=p0)[0] for i in range(N_bts)])
    return par
####

def get_r0_QED2p1(b, sigma, c):
    """QED in 2+1 dimensions: finding the Sommer parameter r0 by the condition:

    r_0^2 * F(r_0) = c

    where F(r) = \frac{d}{dr}V(r). 
    In QED_{2+1} V(r) = V_0 + b*log(r) + sigma*r, hence r0 is found solving:

    sigma*r_0^2 + b*r_0 - c = 0 
    """
    r0 = (-b + np.sqrt(b*b + 4*sigma*c))/(2*sigma)
    return r0
####

## example code
if __name__ == "__main__":
    ## reading the configurations
    W_gauge = read_Wloops_1st_format("../build/omeas/", i1=10, step=10, Ng=5, nmax_iter=1000, L=16, T=16)
    Ng = W_gauge.shape[0]
    N_bts = 100
    bts_idx = np.random.randint(0, high=Ng, size=N_bts, dtype=int)
    W = sampling.uncorrelated_confs_to_bts(W_gauge, bts_idx=bts_idx)

    #W = read_Wloops_format2("../build/omeas/", L=16, T=16, beta=1.3, xi=1.0, ss=False, sizeWloops=0.5)
    V = []
    for ir in range(4):
       V.append(get_potential_bts(W[:,:,ir], t1=1, t2=3))
    ###
    V = np.array(V).transpose()
    print(V.shape)
    #
    ## fitting fake data
    L = 16
    N_bts = 1000
    r = np.array([i for i in range(0, L)])
    v0, b, sigma = 0.5, 0.2, 13
    noise = np.random.normal(1.0, 0.001, size=N_bts)
    ## ignore the warning for log(0), r=0 is not included in the fit
    V = np.array([noise[i]*(v0 + b*np.log(r) + sigma*r) for i in range(N_bts)])
    r1, r2 = 1, 3
    par = fit_static_potential_QED2p1(V_bts_all=V, r1=r1, r2=r2)
    r0 = get_r0_QED2p1(par[:,1], par[:,2], c=1.65)
    print("Input for    v0, b, sigma :", v0, b, sigma)
    print("Best fit for v0, b, sigma :", np.average(par, axis=0))
    print("r0 :", np.average(r0, axis=0))
####
