
library(boot)


#' fit a function f(x1,...,xn)
#' N_pts = number of points
#' x = matrix of size n \times N_pts
#' y, dy = vectors of size N_pts
#' maxiter = maximum number of iterations for the minimizer
fit_xiyey <-function(ansatz, x, y, ey,
                     guess, maxiter=10000, method = "BFGS"){
    
    N_pts <- length(y) # number of data points
    N_par <- length(guess) # number of parameters of the fit
    N_dof <- N_pts - N_par # number of degrees of freedom
                                        # chi square residual function
    ch2 <- function(par) {
        ch2_res <- 0
        for (i in 1:N_pts) {
            df <- y[i] - ansatz(x[,i], par) # y_i - f_i
            #print(x[,i])
            #print(y[i])
            #print(ansatz(x[,i], par))
            #q()
            df <- df / ey[i]
            df2 <- df*df
            ch2_res <- ch2_res + df2
        }
        return(ch2_res)
    }

    mini <- optim(
      par=guess,      # Initial guess for the minimum
      fn=ch2,         # function to be minimized
      method = method # Optimization algorithm (e.g., "BFGS", "Nelder-Mead", etc.)
    )

    ch2_value <- mini[["value"]]
    
    res <- list()

    res[["ansatz"]] <- ansatz
    res[["N_pts"]] <- N_pts
    res[["N_par"]] <- N_par
    res[["par"]] <- mini[["par"]]
    res[["ch2"]] <- ch2_value
    res[["dof"]] <- N_dof # degrees of freedom
    res[["ch2_dof"]] <- ch2_value/N_dof

    return(res)
}

#' fit boostrap by boostrap (generated from y and dy)
#' x = matrix of size n \times N_pts
#' returns a list with bootstraps of parameters, their mean, stderr, chi^2, reduced chi^2
bootstrap.fit_xiyey <-function(ansatz, x, y, ey, guess, maxiter=10000, method="BFGS", N_bts = 1000){
    N_var <- nrow(x)
    N_pts <- length(y) # number of points
    N_par <- length(guess)
    N_dof = N_pts - N_par

    par_bts <- matrix(data=NA, N_bts, N_par)
    ch2_bts <- matrix(data=NA, N_bts, N_par)
    
    y_bts <- matrix(data=rnorm(N_bts*N_pts, y, ey), ncol=N_pts, byrow=TRUE)

    ansatz <- NA
    for (i in 1:N_bts){
        mini <- fit_xiyey(
            ansatz = ansatz, x=x, y=y_bts[i,], ey=ey,
            guess = guess,
            maxiter=maxiter, method=method)
        
        if (i==1){
            ansatz <- mini[["ansatz"]]
            N_dof <- mini[["N_dof"]]
        }

        par_bts[i,] <- mini[["par"]]
        ch2_bts[i,] <- mini[["ch2"]]
    }

    # mean and standard error on the bootstraps

    par_val <- apply(par_bts, 2, mean)
    par_sd <- apply(par_bts, 2, sd)

    ch2_val <- apply(ch2_bts, 2, mean)
    ch2_sd <- apply(ch2_bts, 2, sd)

    res <- list(
        ansatz=ansatz,
        N_bts=N_bts, N_pts = N_pts, N_var = N_var, N_par=N_par,
        x=x, y=y, ey = ey, y_bts = y_bts,
        par=list(bts=par_bts, val=par_val, dval=par_sd),
        ch2=list(bts=ch2_bts, val=ch2_val, dval=ch2_sd),
        ch2_dof=list(
            bts=ch2_bts/N_dof, val=ch2_val/N_dof,
            dval=ch2_sd/N_dof)
        )

    return(res)
}


#' uses the values of the minimization to predict the yi
#' x = matrix of size n \times N_pts
bootstrap.fit_xiyey.predict <-function(mini, x){

    N_pts <- ncol(x) # number of data points
    N_bts <- mini[["N_bts"]] # number of bootstraps
 
    ansatz <- mini[["ansatz"]] # ansatz function
    y <- matrix(data=NA, nrow=N_bts, ncol=N_pts)
    
    par_bts <- mini[["par"]][["bts"]]
    for (i in 1:N_bts){
        for (j in 1:N_pts){
            y[i,j] <- ansatz(x[,j], par_bts[i,])
        }
    }

    res <- list(mini = mini,
                x = x, bts = y,
                value=apply(y, 2, mean), dvalue=apply(y, 2, sd)
                )
    return(res)
}

#' fit_xiyey +  prediction for the values of x0
#' x0 = matrix of size n \times Nd, where Nd is the number of points we want to interpolate/extrapolate to
bootstrap.fit_xiyey_and_predict <- function(ansatz, x, y, ey,
                                            guess, maxiter=10000, method="BFGS",
                                            N_bts = 1000,
                                            x0)
{
    mini <- bootstrap.fit_xiyey(ansatz=ansatz, x=x, y=x, ey=ey,
                                guess=guess, maxiter=maxiter, method=method,
                                N_bts = N_bts)

    pred <- bootstrap.predict_from_fit(mini=mini, x=x0)

    res <- list(prediction = pred)
    return(res)
}
