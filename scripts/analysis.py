""" Analysis routines """

import numpy as np
import matplotlib.pyplot as plt

from . import effective_curves as ec
from . import gevp

def bts_corr_plot(
        C_bts: np.ndarray, 
        t1: int, t2: int, 
        outfile: str, 
        yscale=None, label=None, linestyle="None"):
    """ Plot of the correlator from the bootstrap samples in the interval t1, t2 included """
    N_bts, T_ext = C_bts.shape
    print("## Potting the correlator")
    plt.errorbar(
        x=np.arange(t1, t2+1, 1), y=np.average(C_bts, axis=0)[t1:t2+1], 
        yerr=np.std(C_bts, axis=0)[t1:t2+1], 
        label=label,
        linestyle = linestyle,
        capsize=2,
        marker = "o",
        markerfacecolor='None'
        )
    plt.ylabel("C(t)")
    if label != None:
        plt.legend()
    ####
    if yscale == "log":
        plt.yscale("log")
    ####
    print("## output:", outfile)
    plt.savefig(outfile)
    plt.close()
####


def bts_meff(C_bts: np.ndarray, strategy: str, T: int, t1: int, t2: int, plot=False, outfile=None, **kwargs):
    """ 
    Bootstrap analysis of the effective mass. 
    Returns the bootstraps of the best fit value. 
    
    - The fit is done by fitting to a constant the effective mass curve in the interval [t1,t2] (both included)
    - The time extent T and the strategy for the effective mass curve determination are provided by the user
    - If outfile is passed, is used as path for the output plot
    
    """
    N_bts = C_bts.shape[0]
    print("## Finding the effective mass")
    M_eff = np.array([ec.get_m_eff(C_bts[ib,:], strategy=strategy, T=T) for ib in range(N_bts)])
    T_ext = M_eff.shape[1] ## time extent of the effective mass curve
    #times_eff = np.array([t for t in range(T-1)])
    M_eff_plat = M_eff[:,t1:(t2+1)]
    dM_eff_plat = np.std(M_eff_plat, axis=0, ddof=1)
    print("## Fitting the effective mass")
    print("T=",T, "t1=", t1, "t2=", t2)
    M_eff_fit = np.array([ec.fit_eff_mass(M_eff_plat[ib,:], dM_eff_plat) for ib in range(N_bts)])
    dM_eff_fit = np.std(M_eff_fit, axis=0, ddof=1)
    if plot:
        print("Plotting")
        M_eff_fit_avg = np.average(M_eff_fit, axis=0)
        ## plotting the fit band
        plt.fill_between(
            np.arange(t1, t2+1, 1), 
            M_eff_fit_avg+dM_eff_fit, 
            M_eff_fit_avg-dM_eff_fit,
            alpha=0.5
            )
        plt.plot(
            np.arange(t1, t2+1, 1), np.full(shape=(t2-t1+1), fill_value=M_eff_fit_avg), 
            linestyle="--")
        ## plotting the data
        tmin_plot, tmax_plot = 0, T_ext
        if "tmin_plot" in kwargs.keys():
            tmin_plot = kwargs["tmin_plot"]
        ####
        if "tmax_plot" in kwargs.keys():
            tmax_plot = kwargs["tmax_plot"]+1
        ####
        label = "$M={M} \\pm {dM}$".format(
            M=np.round(np.average(M_eff_fit, axis=0), 3),
            dM=np.round(np.std(M_eff_fit, axis=0, ddof=1), 3)
            )
        if "label" in kwargs.keys():
            label = kwargs["label"]
        ####
        plt.errorbar(
            x = np.arange(tmin_plot, tmax_plot),
            y = np.average(M_eff, axis=0)[tmin_plot:tmax_plot], 
            yerr=np.std(M_eff, axis=0)[tmin_plot:tmax_plot], 
            label=label,
            capsize=2, marker="o", markerfacecolor='None'
            )
        if outfile != None:
            print("output:", outfile) 
            plt.ylabel("$M_{\mathrm{eff}}(t)$")
            plt.legend()
            plt.savefig(outfile)
            plt.close()
    ####
    return M_eff_fit
####

def bts_gevp(C_bts: np.ndarray, t0: int, outfile = None):
    N_bts = C_bts.shape[0]
    assert C_bts.shape[1] == C_bts.shape[1]
    N, T_ext = C_bts.shape[2], C_bts.shape[3]
    res = np.array([gevp.gevp(C_bts[i,:], t0=t0) for i in range(N_bts)])
    return res
####
