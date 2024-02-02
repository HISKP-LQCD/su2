
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from scipy.optimize import root
from scipy.optimize import curve_fit

"""
Match the Lagrangian and Hamiltonian coupling beta using the curve y(beta):

Steps:

- Draw y(beta) for the Lagrangian and Hamiltonian, and fit according to the ansatz
- At each y0, find the pair (beta_L, beta_H) on the two curves
- Draw the obtained numerical curve beta_H(beta_L)

Parameters
----------
beta_L : np.array 
    Lagrangian couplings
y_L : np.array
    Lagrangian values of the observable y
dy_L: np.array
    Uncertainties on y_L
beta_H : np.array 
    Hamitonian couplings
y_H : np.array
    Hamiltonian values of the observable y
N_bts: int
    Number of bootstraps (propagation of the uncertainties on y_L
ansatz_beta_VS_y: function
    Fit ansatz for the functional dependence y(beta)
p0: list
    Initial guesses for the fit parameters
**kwargs: **dict
   Pointer to dictionary (see example)

Returns
-------
dict:
    A dictionary with the initial data and the bootstraps of the numerical matching curve coordinates

Example:
--------

kwargs = dict({
    "name_obs" : "mass gap",
    "label_obs" : "$M$",
    "output_beta_fit" : "beta_VS_mass_gap-fit.pdf",
    "output_matching_curve" : "matching_curve-from_mass_gap.pdf",
    "matching_curve_plot_title" : "$\\beta$ match using the mass gap"
})

match_LH(beta_L, y_L, dy_L,
         beta_H, y_H,
         N_bts, ansatz_beta_VS_y, p0,
         **kwargs)

"""
def match_LH(
        beta_L, y_L, dy_L,
        beta_H, y_H,
        N_bts, ansatz_beta_VS_y, p0,
        **kwargs):

    name_obs = kwargs["name_obs"]
    label_obs = kwargs["label_obs"]

    beta_fit_plot_title =  name_obs + " fit"
    output_beta_fit = kwargs["output_beta_fit"]
    output_matching_curve = kwargs["output_matching_curve"]
    matching_curve_plot_title = kwargs["matching_curve_plot_title"]

    plt.errorbar(
        beta_L, y_L, yerr = dy_L,
        linestyle="None", color="midnightblue", capsize=0.5,
        label = "Lagrangian data"
    )

    plt.scatter(
        beta_H, y_H,
        linestyle="None", color="red", s=0.2,
        label="Hamiltonian data"
    )

    ## Fitting y(beta)

    y_L_bts = [] ##, y_H_bts = [], []
    for _ in range(N_bts):
        y_L_bts.append(y_L + np.random.normal(0, dy_L))
        ##y_H_bts.append(y_H)    
    ####


    n_par = len(p0)
    pi_L_bts, pi_H_bts = [], [] # fit parameters
    ch2_red_L_bts, ch2_red_H_bts = [], [] # reduced chi square
    poly_L_bts, poly_H_bts = [], [] # best fit ansatz for each bootstrap
    xd_L_bts, yd_L_bts = [], []
    xd_H_bts, yd_H_bts = [], []

    for i in range(N_bts):
        pi_L, _ = curve_fit(ansatz_beta_VS_y, beta_L, y_L_bts[i], p0=p0, sigma = dy_L) # Lagrangian
        pi_H, _ = curve_fit(ansatz_beta_VS_y, beta_H, y_H, p0=p0) # Hamiltonian

        pi_L_bts.append(pi_L)
        pi_H_bts.append(pi_H)

        ch2_L_num = np.sum(((ansatz_beta_VS_y(beta_L, *pi_L) - y_L_bts[i])/dy_L) ** 2)
        ch2_L_den = (len(beta_L)-n_par)
        ch2_red_L = ch2_L_num/ch2_L_den
        ch2_red_H = np.sum(
            ((ansatz_beta_VS_y(beta_L, *pi_H) - y_L_bts[i]) ** 2)/(len(beta_L)-n_par)
        )

        ch2_red_L_bts.append(ch2_red_L)
        ch2_red_H_bts.append(ch2_red_H)

        beta_min, beta_max = min(beta_H), max(beta_H)
        Nd = 1000 # number of dense points

        xd_L = np.linspace(beta_min, beta_max, Nd)
        xd_L_bts.append(xd_L)
        yd_L = ansatz_beta_VS_y(xd_L, *pi_L)
        yd_L_bts.append(yd_L)


        xd_H = np.linspace(beta_min, beta_max, Nd)
        xd_H_bts.append(xd_H)
        yd_H = ansatz_beta_VS_y(xd_H, *pi_H)
        yd_H_bts.append(yd_H)
    ####

    xd_L_mean = np.mean(xd_L_bts, axis=0)
    xd_L_std = np.std(xd_L_bts, axis=0)
    yd_L_mean = np.mean(yd_L_bts, axis=0)
    yd_L_std = np.std(yd_L_bts, axis=0)

    xd_H_mean = np.mean(xd_H_bts, axis=0)
    xd_H_std = np.std(xd_H_bts, axis=0)
    yd_H_mean = np.mean(yd_H_bts, axis=0)
    yd_H_std = np.std(yd_H_bts, axis=0)

    pi_L_mean = np.mean(np.array(pi_L_bts), axis=0)
    pi_H_mean = np.mean(np.array(pi_H_bts), axis=0)

    ch2_red_L_mean = np.mean(ch2_red_L_bts, axis=0)
    ch2_red_L_std = np.std(ch2_red_L_bts, axis=0)
    ch2_red_H_mean = np.mean(ch2_red_H_bts, axis=0)

    plt.fill_between(
        xd_L_mean, yd_L_mean - yd_L_std, yd_L_mean + yd_L_std,
        linestyle="None",
        label = "Lagrangian fit",
        color="cornflowerblue"
    )

    plt.errorbar(
        xd_H_mean, yd_H_mean,
        label="Hamiltonian fit", #"Hamilt. fit: $\\chi^2={ch2_red_H}$".format(ch2_red_H=round(ch2_red_H,2)),
        color="salmon"
    )

    ax = plt.gca()
    ax.tick_params(axis="x", direction='in')
    ax.tick_params(axis="y", direction='in')
    plt.grid(linewidth=0.5)

    plt.title(beta_fit_plot_title)
    plt.xlabel("$\\beta$")
    plt.ylabel(label_obs)
    plt.legend()

    ax = plt.gca()
    #ax.set_yscale('log')

    plt.savefig(output_beta_fit)
    plt.cla()
    plt.clf()
    plt.close()

    ## matching the 2 curves

    # find x0 such that f(x0)=y0
    def find_x0(f, y_0, initial_guess=0.0):
        # Function g(x) = f(x) - y_0
        g = lambda x: f(x) - y_0

        # Find the root of g(x)
        result = root(g, initial_guess)

        if result.success:
            x_0 = result.x[0]
            return x_0
        else:
            raise ValueError("Root not found. Try a different initial_guess.")
        ####
    ####

    ## interval of matching on the y
    y_match = np.linspace(
        ansatz_beta_VS_y(min(beta_L), *pi_L_mean), ansatz_beta_VS_y(max(beta_H), *pi_H_mean), Nd)

    beta_L_match_bts, beta_H_match_bts = [], []
    for i in range(N_bts):
        beta_L_match, beta_H_match = [], []
        fi_L = lambda x: ansatz_beta_VS_y(x, *pi_L_bts[i])
        fi_H = lambda x: ansatz_beta_VS_y(x, *pi_H_bts[i])
        ####
        for p0 in y_match:
            b_L = find_x0(fi_L, p0)
            beta_L_match.append(b_L)
            #
            b_H = find_x0(fi_H, p0)
            beta_H_match.append(b_H)
        ####
        beta_L_match = np.array(beta_L_match)
        beta_L_match_bts.append(beta_L_match)
        #
        beta_H_match = np.array(beta_H_match)
        beta_H_match_bts.append(beta_H_match)
    ####

    beta_L_match_bts_mean = np.mean(beta_L_match_bts, axis=0)
    beta_L_match_bts_std = np.std(beta_L_match_bts, axis=0)
    beta_H_match_bts_mean = np.mean(beta_H_match_bts, axis=0)
    beta_H_match_bts_std = np.std(beta_H_match_bts, axis=0)

    plt.errorbar(
        beta_L_match_bts_mean, beta_H_match_bts_mean,
        xerr = beta_L_match_bts_std, yerr = beta_H_match_bts_std, 
        label="Matching curve from "+name_obs)


    ax = plt.gca()
    ax.tick_params(axis="x", direction='in')
    ax.tick_params(axis="y", direction='in')
    plt.grid(linewidth=0.5)

    plt.title(matching_curve_plot_title)
    plt.xlabel("$\\beta_L$")
    plt.ylabel("$\\beta_H$")
    plt.legend()
    plt.savefig(output_matching_curve)

    res = dict({
        "beta_L": beta_L, "y_L": y_L, "dy_L": dy_L,
        "beta_H": beta_H, "y_H": y_H,
        "beta_L_match_bts": beta_L_match_bts, "beta_H_match_bts": beta_H_match_bts
    })
    
    return res
####
