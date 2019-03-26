# use simulated data to verify the coverage obtained using our method

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypergeom, binom, uniform
import pickle as pickle
#import datetime
from scipy.special import logit


import sys
sys.path.insert(0,'../../..') # allows us to import undetected extinctions package

from undetected_extinctions.undetected_extinctions import SE_changed_E, inverse_midp

def inverse_tak1(alpha, S0, S1, U1, d0):
    '''
    U0_bnd = tak1(alpha, S0, S1, U1, d0)

    Find the bound for the number of undetected extant species at the previous timestep (U0)
    given the data at a given timestep (S1, U1, d0 etc.) and confidence level using the
    hypergeometric.

    (my modification of Tak's method, for comparing with the performance of midp)

    NOTE: using this function can cause the main function to blow out wildly for modest
    values of U0 (U0=300)

    alpha:
        float, confidence level (e.g. obtained by randomly sampling ~ U(0,1))
    S0:
        integer, the number of detected extant species at the previous timestep
    S1:
        integer, the number of detected extant species at the current timestep
    U1:
        integer, the number of undetected extant species at the current timestep
    d0:
        integer, the number of species detected during the previous timestep
    '''

    eps = 1e-10 # a numerical tolerance

    # define the function that we will be using
    # ---

    tak1_fnc = lambda U0, S0, S1, U1, d0: hypergeom.cdf( U1+d0, S0+U0, U0, S1+U1 )


    # obtain a sample value of U0 at our confidence level alpha
    # ---

    min_poss_U0 = U1+d0

    # first, check if we're in the situation where we accept the minimum possible value of U0
    # the second condition on alpha concerns the numerical issue below
    if ( tak1_fnc(min_poss_U0+1, S0, S1, U1, d0) < alpha ) or (1-alpha < eps):

        U0_bnd = min_poss_U0 # U0 is set to the minimum value, end search

    else:

        U0_next = min_poss_U0+1

        # there is a numerical issue for low min_poss_U0 towards end of time series
        # where tak1_fnc returns values near 1 that don't behave correctly
        # need to step out of this region
        while 1-tak1_fnc(U0_next, S0, S1, U1, d0) < eps:
            U0_next += 1

        # the function we'll use here on out, logit because slope changes less and it's unbounded
        logit_alpha = logit(alpha)
        f = lambda U0: logit( tak1_fnc(U0, S0, S1, U1, d0) ) - logit_alpha

        # initialise for my loop
        f_next = f(U0_next)
        notfound_U0_bnd = True

        # loop using the Newton's method to find sample value
        while notfound_U0_bnd:

            # update position

            U0_prev = U0_next; f_prev = f_next

            # do a Newton's step

            df = f(U0_prev+1) - f_prev              # get slope at new position
            U0_next = int( U0_prev - f_prev / df )  # next position
            f_next = f(U0_next)                     # value at next position

            # check if I've found it

            f_next_befor = f(U0_next-1); f_next_after = f(U0_next+1);

            if np.sign(f_next_befor) != np.sign(f_next):

                U0_bnd = U0_next-1 # always the smaller because this function decreases w increasing U0
                notfound_U0_bnd = False
        
            if np.sign(f_next) != np.sign(f_next_after):

                U0_bnd = U0_next
                notfound_U0_bnd = False

    return U0_bnd

# simulates one instance of a SEUX outcome
def simulate1(U0, S0, mu, detn, n=None):
    '''
    simulate1(U0, S0, mu, detn, n=None)

    Simulates one possible SEUX outcome. Has extinction probability and no. spp detected each timestep constant in time.

    Parameters
    ----------
    U0: integer
        Initial number undetected species.
    S0: integer
        Initial number detected species.
    mu: float between 0 and 1
        Constant extinction probability, used to determine n
    detn: integer
        The constant no. of species detected each timestep
    n: list of integers, optional
        Number of species who survive each timestep (i.e. the hypergeometric sample size)

    Returns
    -------

    S, E, U, X, n: np arrays
        The no. of detected extant, detected extinct, undetected extant, undetected extinct, and survivors, at each timestep
    '''


    # obtain the number of species who survive each timestep
    # ---

    if n is None:

        # the number of species who survive each timestep is fixed for the experiment over which we obtain coverage
        n0 = U0+S0

        ni = n0
        n = list()
        t = 0
        tmax = U0 // detn

        while (ni > 0) and (t <= tmax):

            n.append(ni)
            ni = binom(ni, 1-mu).rvs()
            t += 1

        n = np.array(n)


    # simulate the population subject to this n
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

    return S, E, U, X, n


if __name__ == "__main__":


    # user set parameter values
    # ---

    nreps = 1000  # number of repetitions
    nsamples = 500 # number of samples from which to construct the CIs

    pcileV = np.array([90, 80, 70, 60, 50, 40, 30, 20, 10]) # list of percentiles to do

    U0 = 150    # initial number of undetected species
    S0 = 20     # initial number of detected species

    mu = 0.1 # probability of going extinct each timestep
    detn = 2 # no. species detected each year

    type_of_ci = 'midp'
    # type_of_ci = 'tak1'

    # store parameters all together
    params = {
            'nsamples': nsamples,
            'nreps': nreps,
            'pcileV': pcileV,
            'U0': U0,
            'S0': S0,
            'mu': mu,
            'detn': detn,
            'type_of_ci': type_of_ci
            }


    # repeated simulations
    # ---

    # treatment of percentile depends on if we're doing one or two-sided
    edge_pV = (1-pcileV/100)/2

    cnt_withinV = np.zeros( edge_pV.shape, dtype=int )
    U0_meanV = list()
    n = None

    for nrep in range(nreps):

        print('doing rep ' + str(nrep) )


        # get the simulation (we can reuse the n we obtain)

        S_orig, E_orig, U_orig, X_orig, n = simulate1(U0, S0, mu, detn, n)

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

                if type_of_ci == 'midp':

                    min_poss_U0 = max( ( min_poss_UV[t-1], U1+d0 ) )
                    U[t-1] = inverse_midp(alpha, min_poss_U0, S0, S1, U1, d0, None, None)

                if type_of_ci == 'tak1':

                    U[t-1] = inverse_tak1(alpha, S0, S1, U1, d0)

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
        f = open('../../../results/classical/verify/' + fWriteName, 'wb')

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
