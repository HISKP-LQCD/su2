## reimplementation of hadron uwerr function in python

'''
  Time Series Analysis With Gamma Method
 
  Analyse time series data with the so called gamma method:
  https://arxiv.org/pdf/hep-lat/0306017.pdf
  

  Parameters:
    data_rep : list of nrep replica arrays of data to be analysed
           Each element must be of dimension (N x Nalpha) (i.e. N rows and Nalpha columns), where N is the total number of measurements and Nalpha is the number of (primary) observables
    f : callable function f(x, **kwargs)
        x is the data vector of length N_alpha.
        f may return a vector object ot numeric type
        If not given it is assumed that a primary quantity is analysed: f(x)=x
    S:  initial guess for the ratio tau/tauint, with tau the exponetial autocorrelation length.
    pl (bool): if True, the autocorrelation function, the integrated autocorrelation time as function of the integration cut-off and (for primary quantities) the time history of the observable are plotted
    **kwargs (dict): arguments passed to function f

  Returns:
    A dictionary with the following elements:
     - value
     - dvalue
     - ddvalue
     - tauint
     - dtauint
     - Wopt
     - Wmax
     - tauintofW
     - dtauintofW
     - Qval
     - S
     - N
     - R
     - nrep
     - data
     - Gamma
     - dGamma
     - primary

  Example:
'''
def uwerr(data_rep, f = lambda x: x, S, pl):
    nrep = len(data_rep)
    for
    
    N, N_alpha
