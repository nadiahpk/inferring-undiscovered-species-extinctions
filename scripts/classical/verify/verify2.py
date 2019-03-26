# use simulated data to verify the coverage obtained using our method
# the same as verify.py except I allow the number of species detected at each timestep to vary

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypergeom, binom, uniform
import pickle as pickle
#import datetime
from scipy.special import logit


import sys
sys.path.insert(0,'../../..') # allows us to import undetected extinctions package

from undetected_extinctions.undetected_extinctions import SE_changed_E, inverse_midp


# simulates one instance of a SEUX outcome
def simulate2(U0, S0, mu, detn_array, n=None, d=None):
    '''
    simulate2(U0, S0, mu, detn_array, n=None, d=None)

    Simulates one possible SEUX outcome. Has extinction probability and no. spp detected each timestep constant in time.

    Parameters
    ----------
    U0: integer
        Initial number undetected species.
    S0: integer
        Initial number detected species.
    mu: float between 0 and 1
        Constant extinction probability, used to determine n
    detn: array of integers
        The no. of species detected each timestep (will be randomly chosen)
    n: list of integers, optional
        Number of species who survive each timestep (i.e. the hypergeometric sample size)
    d: list of integers, optional
        Number of species detected each timestep 

    Returns
    -------

    S, E, U, X, n, d: np arrays
        The no. of detected extant, detected extinct, undetected extant, undetected extinct, survivors, and detections at each timestep
    '''


    # obtain the number of species who survive each timestep
    # ---

    if n is None:


        # the number of species who survive each timestep is fixed for the experiment over which we obtain coverage
        n0 = U0+S0

        ni = n0
        n = list()

        di = np.random.choice( detn_array )
        d = list() # a place to store the detn that were randomly chosen

        tmax = U0 // min(detn_array) # the longest possible
        t = 0
        while (ni > 0) and (t <= tmax):

            n.append(ni)
            ni = binom(ni, 1-mu).rvs()

            d.append(di)
            di = np.random.choice( detn_array )

            t += 1

        n = np.array(n)
        d = np.array(d)


    # simulate the population subject to this n and d
    # ---

    undetected_remain = True # initialise as saying that there are still undetected species to be found

    S = [S0]
    E = [0]
    U = [U0]
    X = [0]
    t = 1

    while undetected_remain:

        # survival process
        Ut = hypergeom( S[t-1] + U[t-1], U[t-1], n[t] ).rvs()
        St = n[t] - Ut

        Et = E[t-1] + S[t-1] - St
        Xt = X[t-1] + U[t-1] - Ut

        # detection process
        detn = d[t-1]
        if detn >= Ut:

            # detect whatever remains
            St = St+Ut
            Ut = 0
            undetected_remain = False

        else:

            St = St+detn
            Ut = Ut-detn
        

        # append
        S.append(St); U.append(Ut); E.append(Et); X.append(Xt)

        t += 1

    # turn into numpy arrays
    S = np.array(S)
    E = np.array(E)
    U = np.array(U)
    X = np.array(X)

    return S, E, U, X, n, d


if __name__ == "__main__":


    # user set parameter values
    # ---

    nreps = 1000  # number of repetitions
    nsamples = 500 # number of samples from which to construct the CIs

    pcileV = np.array([90, 80, 70, 60, 50, 40, 30, 20, 10]) # list of percentiles to do

    U0 = 120    # initial number of undetected species
    S0 = 20     # initial number of detected species

    mu = 0.05 # probability of going extinct each timestep
    detn_array = [1,2,3] # no. species detected each year

    type_of_ci = 'midp' # just to make it consistent with verify.py

    # store parameters all together
    params = {
            'nsamples': nsamples,
            'nreps': nreps,
            'pcileV': pcileV,
            'U0': U0,
            'S0': S0,
            'mu': mu,
            'detn_array': detn_array,
            'type_of_ci': type_of_ci
            }


    # repeated simulations
    # ---

    # treatment of percentile depends on if we're doing one or two-sided
    edge_pV = (1-pcileV/100)/2

    cnt_withinV = np.zeros( edge_pV.shape, dtype=int )
    U0_meanV = list()
    n = None
    dfixed = None

    for nrep in range(nreps):

        print('doing rep ' + str(nrep) )


        # get the simulation (we can reuse the n we obtain)

        S_orig, E_orig, U_orig, X_orig, n, dfixed = simulate2(U0, S0, mu, detn_array, n, dfixed)

        S, E = SE_changed_E(S_orig, E_orig)     # collapse timesteps in which there were no detected extinctions
        d = S[1:] - S[:-1] + E[1:] - E[:-1]     # discoveries at each timestep
        psi = S[:-1] - (E[1:]-E[:-1])           # survivors at each timestep
        extns = E[1:] - E[:-1]                  # extinctions at each timestep
        T_idx = len(S)                          # number of timesteps
        tV = list( reversed(range(1,T_idx)) )   # list of timesteps to work backwards through

        # the minimum possible values of U_t are the number detected from t onwards
        min_poss_UV = [ sum(d[t-1:]) for t in range(1,T_idx) ]


        # repeatedly sample to construct the CI

        U0V = list()

        for nsample in range(nsamples):

            U = [0]*T_idx # a place to store this replicate

            for t in tV:

                alpha = uniform.rvs()
                S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]

                min_poss_U0 = max( ( min_poss_UV[t-1], U1+d0 ) )
                U[t-1] = inverse_midp(alpha, min_poss_U0, S0, S1, U1, d0, None, None)


            U0V.append( U[0] ) # store result


        # append mean

        U0_meanV.append( np.mean(U0V) )


        # count how many within intervals

        # construct the CI using percentiles
        CIV = [ ( np.percentile(U0V, (100-pcile)/2), np.percentile(U0V, 100-(100-pcile)/2) ) for pcile in pcileV ]
        cnt_within = np.array([ 1 if U0 >= U0_lo and U0 <= U0_hi else 0 for U0_lo, U0_hi in CIV ])
        cnt_withinV = cnt_withinV+cnt_within


    # plot coverage for each confidence level
    # ---

    if False:

        coverage = 100*cnt_withinV/nreps

        plt.plot( pcileV, pcileV, ls='dotted', color='black')
        plt.scatter( pcileV, coverage, color='black')
        plt.xlabel('nominal coverage desired')
        plt.ylabel('actual coverage obtained')
        plt.grid(True)
        plt.show()


    # save summary results to a pickle file
    # ---

    if True:

        fWriteName = 'verify_' + type_of_ci + '_U0_'  + str(U0) + '_nreps_' + str(nreps) + '_nsamples_' + str(nsamples) + '.pkl'
        f = open('../../../results/classical/verify2/' + fWriteName, 'wb')

        # a string explaining the pickle file
        ss  = 'Created by verify.py\n'
        ss += 'Contains the following:\n'
        ss += '0. ss, str: this string you are reading now.\n'
        ss += '1. params, dict: parameters used for the run.\n'
        ss += '2. n, np array: the sample size at each step.\n'
        ss += '3. cnt_withinV, np array of ints: the count of how many times the true U0 was within the bounds corresponding to pcileV.\n'
        ss += '4. U0_meanV, list of ints: the mean obtained for each rep (rounded).\n'
        ss += '5. U0_medianV, list of ints: the median obtained for each rep (rounded).\n'
        pickle.dump( ss, f ) # 0.
        pickle.dump( params, f )
        pickle.dump( n, f )
        pickle.dump( cnt_withinV, f )
        pickle.dump( U0_meanV, f )
        f.close()
