# use minimiser to improve the c_t points for functions near the AIC valley

import pickle
import numpy as np
from scipy.optimize import minimize

import sys
sys.path.insert(0,'../../..') # allows us to import undetected extinctions package

from undetected_extinctions.redetection_effort.redetection_effort import negLogLikelihoodSpline

# some parameter values
# ---

eps = 1e-6 # small value so that our redetection effort cannot be either 0 nor 1
range_c = 0.1 # how far above and below the c_t values to search for a better fit

results_dir = '../../../results/redetection_effort/tighten/' # where to put results


# which to tighten
# ---

data_dir = '../../../results/redetection_effort/backwards_stepping/' # where are data to use is
frst_no_params = 158 # range of no_params to perform tightening on
last_no_params = 166


# create a list of file name suffixes
# ---

suffixes = np.arange(frst_no_params, last_no_params+2, 2) # no. of parameters is used as suffix, increments by four
# done -- tighten_172.pkl  tighten_176.pkl  tighten_178.pkl  tighten_180.pkl  tighten_184.pkl
# dones = [172, 176, 178, 180, 184]
# suffixes = [ s for s in suffixes if s not in dones ]


# read in each results file and tighten
# ---

# loop over suffixes
for suffix in suffixes:

    print( 'tightening ' + str(suffix) )
    # read in info from that file

    fname = data_dir + 'backwards_stepping_' + str(suffix) + '.pkl'
    f = open(fname, 'rb')
    _ = pickle.load( f )   # explanatory string
    cs = pickle.load( f )   # c_t points defining the spline
    ts = pickle.load( f )   # t points defining the spline
    spps = pickle.load( f ) # redetections info for each species used in the fitting
    f.close()

    no_params = len(cs) + len(ts)

    # bounds to search for each cs

    bnds = list()
    for ci in cs:

        bnd = ( max(eps, ci-range_c), min(1-eps, ci+range_c) )
        bnds.append(bnd)

    # perform the minimisation and find the new c_t values along with their new AIC

    res = minimize( lambda cs: negLogLikelihoodSpline(cs, ts, spps), cs, method='L-BFGS-B', bounds=bnds )

    new_cs = res.x
    new_negLL = res.fun

    new_AIC = 2*no_params + 2*new_negLL

    # save it to a file

    fName = results_dir + 'tighten_' + str(no_params) + '.pkl'
    f = open(fName, 'wb')
    # a string explaining the pickle file
    ss  = 'Created by tighten.py.\n'
    ss += 'Contains the following:\n'
    ss += '0. ss, string: this string you are reading now.\n'
    ss += '1. cs, list of floats: c_t points defining the spline.\n'
    ss += '2. ts: list ints: t points defining the spline.\n'
    ss += '3. spps, dictionary: keys are names, and values are dictionary with frst, last, and redetns.\n'
    ss += '4. AIC, float: the AIC calculated for this spline.\n'
    pickle.dump( ss, f )
    pickle.dump( new_cs, f )
    pickle.dump( ts, f )
    pickle.dump( spps, f )
    pickle.dump( new_AIC, f )
    f.close()
