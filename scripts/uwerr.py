## python translation of the uwerrprimary function from "hadron": https://github.com/HISKP-LQCD/hadron

import pandas as pd
import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def gamma_error(Gamma, N, W, Lambda):
    gamma_err = np.zeros(W + 1)
    Gamma_loc = np.copy(Gamma)
    Gamma_loc[(W+1):(2 * W + W + 2)] = 0.0
    
    for t in range(W + 1):
        k = np.arange(max(1, t - W), t + W + 1)
        gamma_err[t] = np.sum((Gamma_loc[k + t] + Gamma_loc[np.abs(k - t)] - 2 * Gamma_loc[t] * Gamma_loc[k]) ** 2)
        gamma_err[t] = np.sqrt(gamma_err[t] / N)

    gamma_err[0] = 0.001
    return gamma_err


def plot_uwerr_primary(u, output_file, name_obs: str):
    y = u["data"]
    N = y.shape[0]
    t = np.array([i for i in range(N)])
    tauint, dtauint = u["tauintofW"], u["dtauintofW"]
    N_tauint = tauint.shape[0] # number of points for which we can plot Gamma and tau with errors
    Gamma, dGamma = u["Gamma"][0:N_tauint], u["dGamma"][0:N_tauint]
    t_tauint = np.array([i for i in range(N_tauint)])

    # Create a PDF file to save the plots
    ppdf = PdfPages(output_file, "ab")

    # Plot y
    plt.figure()
    plt.plot(t, y)
    plt.xlabel("t")
    plt.ylabel("y")
    plt.title(name_obs)
    ppdf.savefig()
    plt.close()

    # Plot histogram of y
    plt.figure()
    plt.hist(y, bins=int(np.sqrt(N)), density=True)
    plt.title('Statistical distribution')
    ppdf.savefig()
    plt.close()

    # # Plot Gamma with error bars
    plt.figure()
    plt.errorbar(t_tauint, Gamma, yerr=dGamma, fmt='o', capsize=3)
    plt.title("Autocorrelation function")
    plt.xlabel('t')
    plt.ylabel('$\\Gamma(t)$')
    ppdf.savefig()
    plt.close()

    # Plot tauint with error bars
    t_tauint = np.array([i for i in range(tauint.shape[0])])
    plt.figure()
    plt.errorbar(t_tauint, tauint, yerr=dtauint, fmt='o', capsize=3)
    plt.xlabel('t')
    plt.ylabel('$\\tau_{{int}}$')
    plt.title("Integrated autocorrelation time")
    ppdf.savefig()
    plt.close()

    ppdf.close()
####


def uwerr_primary(data, nrep=None, S=1.5, output_file=None):
    N = data.shape[0] #len(data)

    if nrep is None:
        nrep = [N]

    if any(np.array(nrep) < 1) or sum(nrep) != N:
        raise ValueError("Error, inconsistent N and nrep!")

    R = len(nrep)
    mx = np.mean(data)
    mxr = np.zeros(R)
    i0 = 0

    for r in range(R):
        i1 = i0 + nrep[r]
        mxr[r] = np.mean(data[i0:i1])
        i0 = i1

    Fb = np.sum(mxr * nrep) / N
    delpro = data - mx

    if S == 0:
        Wmax = 0
        Wopt = 0
        flag = False
    else:
        Wmax = int(np.floor(min(nrep) / 2))
        Gint = 0.0
        flag = True

    GammaFbb = np.zeros(Wmax + 1)
    GammaFbb[0] = np.mean(delpro ** 2)

    if GammaFbb[0] == 0:
        raise ValueError("Error, no fluctuations!")

    W = 1
    while W <= Wmax:
        GammaFbb[W] = 0.0
        i0 = 0

        for r in range(R):
            i1 = i0 + nrep[r]
            GammaFbb[W] += np.sum(delpro[i0:i1 - W] * delpro[i0 + W:i1])
            i0 = i1

        GammaFbb[W] /= (N - R * W)

        if flag:
            Gint += GammaFbb[W] / GammaFbb[0]

            if Gint < 0:
                tauW = 5e-16
            else:
                tauW = S / np.log((Gint + 1) / Gint)

            gW = np.exp(-W / tauW) - tauW / np.sqrt(W * N)

            if gW < 0:
                Wopt = W
                Wmax = min(Wmax, 2 * Wopt)
                flag = False

        W += 1

    if flag:
        print("Warning: Windowing condition failed!")
        Wopt = Wmax

    CFbbopt = GammaFbb[0] + 2 * np.sum(GammaFbb[1:(Wopt+1)])

    if CFbbopt <= 0:
        raise ValueError("Gamma pathological: error^2 < 0")

    GammaFbb += CFbbopt / N  # Correct for bias
    dGamma = gamma_error(Gamma=GammaFbb, N=N, W=Wmax, Lambda=100)
    CFbbopt = GammaFbb[0] + 2 * np.sum(GammaFbb[1:(Wopt+1)])  # Refined estimate
    sigmaF = np.sqrt(CFbbopt / N)  # Error of F
    tauintFbb = np.cumsum(GammaFbb) / GammaFbb[0] - 0.5  # Normalized autocorrelation time

    if R > 1:
        bF = (Fb - mx) / (R - 1)
        mx -= bF

        if np.abs(bF) > sigmaF / 4:
            print(f"A {bF/sigmaF:.1f} sigma bias of the mean has been cancelled")

        mxr -= bF * N / np.array(nrep)
        Fb -= bF * R

    value = mx
    dvalue = sigmaF
    ddvalue = dvalue * np.sqrt((Wopt + 0.5) / N)
    dtauintofW = tauintFbb[0:(Wmax+1)] * np.sqrt(np.arange(Wmax + 1) / N) * 2

    tauint = tauintFbb[Wopt]
    dtauint = tauint * 2 * np.sqrt((Wopt - tauint + 0.5) / N)

    # Q value for replica distribution if R >= 2
    if R > 1:
        chisqr = np.sum((mxr - Fb) ** 2 * np.array(nrep)) / CFbbopt
        Qval = 1 - chi2.cdf(chisqr / 2, (R - 1) / 2)
    else:
        Qval = None

    res = {
        "value": value,
        "dvalue": dvalue,
        "ddvalue": ddvalue,
        "tauint": tauint,
        "dtauint": dtauint,
        "Wopt": Wopt,
        "Wmax": Wmax,
        "tauintofW": tauintFbb[:Wmax + 1],
        "dtauintofW": dtauintofW[:Wmax + 1],
        "Qval": Qval,
        "S": S,
        "N": N,
        "R": R,
        "nrep": nrep,
        "data": data,
        "Gamma": GammaFbb,
        "dGamma": dGamma,
        "primary": 1,
    }

    if output_file != None:
        plot_uwerr_primary(u=res, output_file=output_file, name_obs="Observable")


    return res
####


# example code
# data = pd.read_csv("./markov_chain_data.txt", header=None).to_numpy().transpose()
# u = uwerr_primary(data[0,:], output_file="plots.pdf")


