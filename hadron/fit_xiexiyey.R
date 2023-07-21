library(boot)


#' fit a function f(x1,...,xn) with errors on both x and y
#' N_pts = number of points
#' x = matrix of size n \times N_pts
#' ex = matrix of size n \times N_pts
#' y, dy = vectors of size N_pts
#' guess = vector of guesses:
#'   -  the first np_ansatz are the guesses for the ansatz itself
#'   -  while the remaining values are guesses for the matrix elements of x_hat, unrolled by row
#' maxiter = maximum number of iterations for the minimizer
fit_xiexiyey <- function(ansatz, x, ex, y, ey,
                         guess, maxiter = 10000, method = "BFGS") {
    n_var <- nrow(x) # number of variables
    N_pts <- length(y) # number of data points
    N_par <- length(guess) # number of parameters of the fit
    N_dof <- N_pts - N_par # number of degrees of freedom

    # number of parameters of the ansatz
    np_ansatz <- length(guess) - n_var * N_pts
    if (length(guess) <= n_var * N_pts) {
        print("Error: insufficient number of guess parameters")
        q()
    }

    # chi square residual function
    ch2 <- function(p_all) {
        ch2_res <- 0
        p <- p_all[1:np_ansatz]
        x_hat <- matrix(tail(p_all, n_var * N_pts), nrow = n_var, ncol = N_pts, byrow = TRUE)
        for (i in 1:N_pts) {
            df1 <- y[i] - ansatz(x_hat[, i], p) # y_i - f_i
            df1 <- df1 / ey[i]
            ch2_res <- ch2_res + (df1 * df1)

            for (k in 1:n_var) {
                df_k <- x[k, i] - x_hat[k, i]
                df_k <- df_k / ex[k, i]
                ch2_res <- ch2_res + (df_k * df_k)
            }
        }
        return(ch2_res)
    }

    mini <- optim(
        par = guess, # Initial guess for the minimum
        fn = ch2, # function to be minimized
        method = method # Optimization algorithm (e.g., "BFGS", "Nelder-Mead", etc.)
    )

    ch2_value <- mini[["value"]]

    res <- list()

    res[["ansatz"]] <- ansatz
    res[["N_par"]] <- N_par
    res[["par"]] <- mini[["par"]]
    res[["ch2"]] <- ch2_value
    res[["dof"]] <- N_dof # degrees of freedom
    res[["ch2_dof"]] <- ch2_value / N_dof

    return(res)
}

#' fit boostrap by boostrap (generated from y and dy)
#' x = matrix of size n \times N_pts
#' returns a list with bootstraps of parameters, their mean, stderr, chi^2, reduced chi^2
bootstrap.fit_xiyey <- function(
    ansatz, x, ex, y, ey, guess,
    maxiter = 10000, method = "BFGS", N_bts = 1000) {
    N_pts <- length(y) # number of points
    N_par <- length(guess)
    N_dof <- N_pts - N_par

    par_bts <- matrix(data = NA, N_bts, N_par)
    ch2_bts <- matrix(data = NA, N_bts, N_par)

    y_bts <- matrix(data = rnorm(N_bts * N_pts, y, ey), ncol = N_pts, byrow = TRUE)

    ansatz <- NA
    for (i in 1:N_bts) {
        mini <- fit_xiexiyey(
            ansatz = ansatz, x = x, ex = ex, y = y_bts[i, ], ey = ey,
            guess = guess,
            maxiter = maxiter, method = method
        )

        if (i == 1) {
            ansatz <- mini[["ansatz"]]
            N_dof <- mini[["N_dof"]]
        }

        par_bts[i, ] <- mini[["par"]]
        ch2_bts[i, ] <- mini[["ch2"]]
    }

    # mean and standard error on the bootstraps

    par_val <- apply(par_bts, 2, mean)
    par_sd <- apply(par_bts, 2, sd)

    ch2_val <- apply(ch2_bts, 2, mean)
    ch2_sd <- apply(ch2_bts, 2, sd)

    res <- list(
        ansatz = ansatz,
        N_bts = N_bts, N_pts = N_pts, N_par = N_par,
        x = x, ex = ex,
        par = list(bts = par_bts, val = par_val, dval = par_sd),
        ch2 = list(bts = ch2_bts, val = ch2_val, dval = ch2_sd),
        ch2_dof = list(
            bts = ch2_bts / N_dof, val = ch2_val / N_dof,
            dval = ch2_sd / N_dof
        )
    )

    return(res)
}


#' uses the values of the minimization to predict the yi
#' x = matrix of size n \times N_pts
bootstrap.fit_xiexiyey.predict <- function(mini, x) {
    n_var <- nrow(x) # number of variables
    N_pts <- ncol(x) # number of data points
    N_bts <- mini[["N_bts"]] # number of bootstraps
    np_ansatz <- length(mini[["par"]][["val"]]) - n_var * N_pts

    ansatz <- mini[["ansatz"]] # ansatz function
    y <- matrix(data = NA, nrow = N_bts, ncol = N_pts)

    par <- mini[["par"]][["bts"]]
    for (i in 1:N_bts) {
        p_all <- par[i, ]
        p <- p_all[1:np_ansatz]
        x_hat <- matrix(tail(p_all, n_var * N_pts), nrow = n_var, ncol = N_pts, byrow = TRUE)
        for (j in 1:N_pts) {
            y[i, j] <- ansatz(x_hat[, j], p)
        }
    }

    res <- list(
        mini = mini,
        x = x, bts = y,
        value = apply(y, 2, mean), dvalue = apply(y, 2, sd)
    )

    return(res)
}

#' fit_xiyey +  prediction for the values of x0
#' x0 = matrix of size n \times Nd
#' Nd is the number of points we want to interpolate/extrapolate to
bootstrap.fit_xiyey_and_predict <- function(ansatz, x, y, ey,
                                            guess, maxiter = 10000, method = "BFGS",
                                            N_bts = 1000,
                                            x0) {
    mini <- bootstrap.fit_xiyey(
        ansatz = ansatz, x = x, y = x, ey = ey,
        guess = guess, maxiter = maxiter, method = method,
        N_bts = N_bts
    )

    pred <- bootstrap.fit_xiexiyey.predict(mini = mini, x = x0)

    res <- list(prediction = pred)
    return(res)
}
