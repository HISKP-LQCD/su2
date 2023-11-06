
import numpy as np
import scipy.optimize as opt


def check_data(x, ex, y, ey):
    b1 = (x.shape[1] == ex.shape[1])
    b2 = (len(y)==len(ey))
    b3 = (x.shape[1] == len(y))
    if not(b1 and b2 and b3):
        raise ValueError("Error: please check the data passed to the fit function")
    ####
####

#' fit a function f(x1,...,xn) with errors on both x and y
#' N_pts = number of points
#' x = np.array of shape (n, N_pts)
#' ex = np.array of size (n, N_pts)
#' y, dy = list of length N_pts
#' guess = list of N_par guesses (for the ansatz)
#' for the x_hat I use the mean values of the x (tested to be unstable if far from this)
#' maxiter = maximum number of iterations for the minimizer
def fit_xiexiyey(ansatz, x, ex, y, ey, guess, maxiter = 10000, method = "BFGS"):
    check_data(x=x, ex=ex, y=y, ey=ey) ## check that data is in the correct format
    n_var = x.shape[0] # number of variables
    N_pts = len(y) # number of data points

    N_par = len(guess) # number of parameters of the fit
    N_dof = N_pts - N_par # number of degrees of freedom

    # number of parameters of the ansatz
    np_ansatz = N_par

    # chi square residual function
    def ch2(p_all):
        ch2_res = 0
        p = p_all[0:np_ansatz]
        p_res = p_all[-(n_var * N_pts):]
        x_hat = np.array([[p_res[j + N_pts*i] for j in range(N_pts)] for i in range(n_var)])
        for i in range(N_pts):
            df1 = y[i] - ansatz(x_hat[:,i], p) # y_i - f_i
            df1 = df1 / ey[i]
            ch2_res = ch2_res + (df1 * df1)            
            for k in range(n_var):
                df_k = x[k, i] - x_hat[k, i]
                df_k = df_k / ex[k, i]
                ch2_res = ch2_res + (df_k * df_k)
            ####
        ####
        return ch2_res
    ####

    guess = list(guess)
    ## ansatz for x_hat=mean values of x
    for iv in range(n_var):
        guess += list(x[iv,:])
    ####
    mini = opt.minimize(fun = ch2, x0 = guess, method = method)

    ch2_value = ch2(mini.x)

    res = dict({})

    res["ansatz"] = ansatz
    res["N_par"] = N_par
    res["par"] = mini.x
    res["ch2"] = ch2_value
    res["dof"] = N_dof ## degrees of freedom
    res["ch2_dof"] = ch2_value / N_dof

    return(res)
####


