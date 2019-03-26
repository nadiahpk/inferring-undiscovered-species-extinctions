# Use simulated data based on the birds study to check the coverage of the method.
# There may be an issue with how the 'impossible' value is treated in the inverse_midp function,
# so the function below (find_U0_bnd) is more careful about ensuring that we do not venture two steps
# into the impossible region.
# It won't matter for the plant data (where timesteps without detected extinction events are collapsed 
# away), but it will matter if someone wants to apply the method to a dataset with long periods of 
# successive timesteps with no extinctions. 

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypergeom, binom, uniform
import pickle as pickle
#import datetime
from scipy.special import logit


import sys
sys.path.insert(0,'../../..') # allows us to import undetected extinctions package

from undetected_extinctions.undetected_extinctions import inverse_midp

def find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag=False):
    '''
    U0_bnd = find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag)

    Find the bound for the number of undetected extant species at the previous timestep (U0)
    given the data at a given timestep (S1, U1, d0 etc.) and confidence level using the
    mid-P method.

    The model used is the central hypergeometric distribution. 

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
    impossibleFlag:
        logical, whether or not our data is in the "impossible" region
    '''

    min_poss_U0 = U1 + d0 # the minimum possible value of U0 in reality

    # define the mid-P function that we will be using
    # ---

    midP_fnc = lambda U0, S0, S1, U1, d0: 0.5 * ( hypergeom.cdf( U1+d0, S0+U0, U0, S1+U1 ) + hypergeom.cdf( U1+d0-1, S0+U0, U0, S1+U1 ) )


    # obtain a sample value of U0 at our confidence level alpha
    # ---

    # first, check if we're in the situation where we wouldn't accept the minimum possible value of U0

    if midP_fnc(min_poss_U0, S0, S1, U1, d0) < alpha:

        if impossibleFlag:

            # don't take two steps in to the impossible region
            U0_bnd = min_poss_U0

        else:

            # U0 is set to the 'impossible' value, one less than minimum value, end search
            U0_bnd = min_poss_U0-1 
            impossibleFlag = True

    else:

        # We are trying to find the greatest value of U0, called U0_bnd, such that
        # midP_fnc(U0_bnd, ...) > alpha. 
        # In other words, we need to find the value such that the mid-P function switches:
        # midP_fnc(U0_bnd, ...) > alpha  and  midP_fnc(U0_bnd+1, ...) < alpha.
        # The search has two phases:
        # 1. find an interval [U0_lo, U0_hi] within which U0_bnd lies
        # 2. binary search between our upper and lower search bounds, U0_lo and U0_hi, to find U0_bnd


        # 1. find an interval [U0_lo, U0_hi] within which U0_bnd lies
        # ----

        # To quickly find the interval, we start at the minimum possible value for U0_bnd,
        # and work our way upwards until we find a U0_hi such that midP_fnc(U0_hi, ...) < alpha.
        # The step-size from U0_lo to the next U0 value is doubled each time, to account
        # for the fact that there is no upper limit to the U0 value

        # initialise for the loop
        U0_lo = min_poss_U0
        step_size = 1 # this will double at each step
        U0_hi_found = False

        while not U0_hi_found:

            # prospective upper
            U0_hi = U0_lo + step_size

            if midP_fnc(U0_hi, S0, S1, U1, d0) < alpha:

                U0_hi_found = True  # we've found an upper bound to search within

            else:

                U0_lo = U0_hi   # we didn't find an upper so this is new lower bound for search
                step_size *= 2  # double the step size for next time


        # 2. binary search between our upper and lower search bounds, U0_lo and U0_hi
        # ----

        # Once an interval [U0_lo, U0_hi] has been identified, we use a binary search
        # to find the value of U0_lo such that 
        # midP_fnc(U0_lo, ...) > alpha  and  midP_fnc(U0_lo+1, ...) < alpha.
        # This value of U0_lo is U0_bnd.
        # The binary search works by iteratively taking the mid-point U0_mid between U0_lo and U0_hi,
        # checking if it is a new low or high bound on the interval, and checking the condition above
        # to terminate the search.

        # have we already found bound?
        U0_bnd_found = U0_hi - U0_lo == 1 # if they're one apart, then we've found the switch point

        while not U0_bnd_found:

            U0_mid = int( (U0_lo+U0_hi) / 2 ) # floored midpoint
            alpha_mid = midP_fnc(U0_mid, S0, S1, U1, d0) # the alpha value at this mid-point

            if alpha_mid == alpha: # it's the actual bound (unlikely to happen)

                U0_bnd_found = True
                U0_lo = U0_mid # the bound is stored in U0_lo

            else:

                if alpha_mid < alpha: # the mid-point is a new upper bound

                    U0_hi = U0_mid

                else: # elif alpha_mid > alpha: # the mid-point is a new lower bound

                    U0_lo = U0_mid

                # if the interval bounds are one apart, then we've found U0_bnd
                U0_bnd_found = U0_hi - U0_lo == 1 

        # When we exit the while loop above then the bound has been found.
        # It is stored in U0_lo (i.e. the greatest value of U0 for which midP_fnc >= alpha)
        U0_bnd = U0_lo

        if U0_bnd > min_poss_U0:
            # we've moved out of the impossible region
            impossibleFlag = False

    if U0_bnd < 0:
        U0_bnd = 0 # don't allow negative numbers of undetected species



    return U0_bnd, impossibleFlag


# simulates one instance of a SEUX outcome, where parameter values are similar-ish to the birds study
def simulate_like_birds():
    '''
    simulate_like_birds()

    Simulates one possible SEUX outcome. 

    Returns
    -------

    S, E, U, X, n, d: np arrays
        The no. of detected extant, detected extinct, undetected extant, undetected extinct, survivors, and detections at each timestep
    '''

    # parameter values, chosen to be similar to birds study (Chishom et al. 2016)
    # ---

    U0 = 48;
    S0 = 167;

    # detections
    d = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 2, 3, 2, 2, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1])

    # survivals
    n = np.array([213, 213, 213, 212, 212, 212, 211, 211, 210, 209, 209, 208, 208, 207, 207, 205, 205, 205, 205, 205, 205, 205, 205, 203, 202, 202, 202, 202, 201, 201, 201, 201, 201, 200, 200, 200, 200, 196, 196, 196, 196, 196, 196, 196, 195, 195, 195, 195, 195, 193, 191, 191, 191, 190, 189, 189, 187, 187, 187, 187, 186, 186, 186, 186, 186, 186, 186, 185, 185, 185, 182, 182, 182, 182, 182, 182, 182, 182, 153, 153, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 146, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 141, 141, 141, 140, 139, 139, 139, 139, 138, 138, 138, 138, 138, 138, 137, 136, 136, 136])

    len_n = len(n)

    # The above have a fixed length of 133 timesteps, however in the random simulation,
    # we may need to go for more than 133 timesteps in order to reach U=0.
    # Therefore, I'll make the assumption that -- after we run out of bird data -- we'll
    # have all-but-one survives and 1 detection per timestep


    # simulate the population
    # ---

    undetected_remain = True # initialise as saying that there are still undetected species to be found

    S = [S0]
    E = [0]
    U = [U0]
    X = [0]
    t = 1

    while undetected_remain:

        if t < len_n:
            n_t = n[t-1]
            detn = d[t-1]
        else:
            n_t = S[t-1] + U[t-1] - 1
            detn = 1

        # survival process
        Ut = hypergeom( S[t-1] + U[t-1], U[t-1], n_t ).rvs()
        St = n_t - Ut

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

    return S, E, U, X


if __name__ == "__main__":


    # user set parameter values
    # ---

    nreps = 1000  # number of repetitions
    nsamples = 500 # number of samples from which to construct the CIs

    pcileV = np.array([90, 80, 70, 60, 50, 40, 30, 20, 10]) # list of percentiles to do

    # type_of_ci = 'midp' # just to make it consistent with verify.py
    type_of_ci = 'fbnd' # new idea about it, prevents two steps in

    # store parameters all together
    params = {
            'nsamples': nsamples,
            'nreps': nreps,
            'pcileV': pcileV,
            'U0': 48, # fixed for the birds simulation
            'S0': 167,
            'mu': None,
            'detn': None,
            'type_of_ci': type_of_ci
            }


    # repeated simulations
    # ---

    # treatment of percentile depends on if we're doing one or two-sided
    edge_pV = (1-pcileV/100)/2

    cnt_withinV = np.zeros( edge_pV.shape, dtype=int )
    U0_meanV = list()

    for nrep in range(nreps):

        print('doing rep ' + str(nrep) )


        # get the simulation 

        S, E, U, X = simulate_like_birds()
        U0 = U[0]

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

            impossibleFlag = False

            for t in tV: # work our way backwards

                alpha = uniform.rvs()
                S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]

                if type_of_ci == 'midp':

                    min_poss_U0 = max( ( min_poss_UV[t-1], U1+d0 ) )
                    U[t-1] = inverse_midp(alpha, min_poss_U0, S0, S1, U1, d0, None, None)

                if type_of_ci == 'fbnd':

                    U[t-1], impossibleFlag = find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag)

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
        f = open('../../../results/classical/verify4/' + fWriteName, 'wb')

        # a string explaining the pickle file
        ss  = 'Created by verify.py\n'
        ss += 'Contains the following:\n'
        ss += '0. ss, str: this string you are reading now.\n'
        ss += '1. params, dict: parameters used for the run.\n'
        ss += '2. None: to match previous.\n'
        ss += '3. cnt_withinV, np array of ints: the count of how many times the true U0 was within the bounds corresponding to pcileV.\n'
        ss += '4. U0_meanV, list of ints: the mean obtained for each rep (rounded).\n'
        ss += '5. U0_medianV, list of ints: the median obtained for each rep (rounded).\n'
        pickle.dump( ss, f ) # 0.
        pickle.dump( params, f )
        pickle.dump( None, f )
        pickle.dump( cnt_withinV, f )
        pickle.dump( U0_meanV, f )
        f.close()
