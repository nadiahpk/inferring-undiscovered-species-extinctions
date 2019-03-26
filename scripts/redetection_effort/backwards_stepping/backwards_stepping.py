# incrementally remove points from the redetection effort spline and calculate AIC for each

import pickle
import numpy as np

import sys
sys.path.insert(0,'../../..') # allows us to import undetected extinctions package

from undetected_extinctions.redetection_effort.redetection_effort import negLogLikelihoodSpline

# some parameters
# ---

# time range over which we'll do the redetection effort
t0 = 1822 # the Wallich collection
tf = 2015 # taken as our last date

results_dir = '../../../results/redetection_effort/backwards_stepping/' # where to put intermediate results


# where databases are located
# ---

# fname_startwith = '../../../results/redetection_effort/initial_fit/initial_fit.pkl'
fname_startwith = results_dir + 'backwards_stepping_150.pkl' # to restart backwards stepping


# obtain the initial fit
# ---

# read in data
f = open(fname_startwith, 'rb')
ss = pickle.load( f )   # explanatory string
cs = pickle.load( f )   # c_t points defining the spline
ts = pickle.load( f )   # t points defining the spline
spps = pickle.load( f ) # redetections info for each species used in the fitting
f.close()

# check the output
assert( ts[0] == t0 )  # t0 = 1822, the Wallich collection
assert( ts[-1] == tf ) # tf = 2015, taken as our last date

# calculate AIC
negLL = negLogLikelihoodSpline(cs, ts, spps)
no_params = len(cs) + len(ts)
AIC = 2*no_params + 2*negLL


# incrementally remove point that provides best AIC, store results
# ---

# storage of each of our steps
tsV = [ ts ] 
csV = [ cs ] 
AICV = [ AIC ] 
no_paramsV = [ no_params ]

while no_params > 120:

    # find the best point to remove

    best_AIC = np.inf # initialise

    for idx in range( 1, len(ts)-1 ): # we never remove endpoints of timeseries

        # new are the same as old ...
        new_ts = [ t for t in ts]
        new_cs = [ c for c in cs]
        # .. but with idx removed
        del new_ts[idx]
        del new_cs[idx]

        # calculate new AIC
        negLL = negLogLikelihoodSpline(new_cs, new_ts, spps)
        no_params = len(new_cs) + len(new_ts)
        new_AIC = 2*no_params + 2*negLL

        if new_AIC < best_AIC:
            best_ts = new_ts
            best_cs = new_cs
            best_AIC = new_AIC


    # store the best found for the next update

    ts = best_ts; cs = best_cs; AIC = best_AIC
    tsV.append( ts ); csV.append( cs ); AICV.append( AIC ); no_paramsV.append( no_params )

    print('no params: ' + str(no_params) + ' AIC: ' + str(AIC) )


    # save it to a file

    fName = results_dir + 'backwards_stepping_' + str(no_params) + '.pkl'
    f = open(fName, 'wb')
    # a string explaining the pickle file
    ss  = 'Created by backwards_stepping.py.\n'
    ss += 'Contains the following:\n'
    ss += '0. ss, string: this string you are reading now.\n'
    ss += '1. cs, list of floats: c_t points defining the spline.\n'
    ss += '2. ts: list ints: t points defining the spline.\n'
    ss += '3. spps, dictionary: keys are names, and values are dictionary with frst, last, and redetns.\n'
    ss += '4. AIC, float: the AIC calculated for this spline.\n'
    pickle.dump( ss, f )
    pickle.dump( cs, f )
    pickle.dump( ts, f )
    pickle.dump( spps, f )
    pickle.dump( AIC, f )
    f.close()
