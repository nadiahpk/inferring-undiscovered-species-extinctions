# find the AIC of the initial fit

import pickle

import sys
sys.path.insert(0,'../../..') # allows us to import undetected extinctions package

from undetected_extinctions.redetection_effort.redetection_effort import negLogLikelihoodSpline


# -------------

# some parameters
# ---

# time range over which we'll do the redetection effort
t0 = 1822 # the Wallich collection
tf = 2015 # taken as our last date


# where databases are located
# ---

fname_initfit = '../../../results/redetection_effort/initial_fit/initial_fit.pkl'


# read in the information from our initial fit
# ---

f = open(fname_initfit, 'rb')
ss = pickle.load( f )   # explanatory string
cs = pickle.load( f )   # c_t points defining the spline
ts = pickle.load( f )   # t points defining the spline
spps = pickle.load( f ) # redetections info for each species used in the fitting
f.close()

# check the output
assert( ts[0] == t0 )  # t0 = 1822, the Wallich collection
assert( ts[-1] == tf ) # tf = 2015, taken as our last date


# calculate AIC
# ---

negLL = negLogLikelihoodSpline(cs, ts, spps)
AIC = 2*( len(cs) + len(ts) ) + 2*negLL
print('AIC = ', AIC)

