# functions for doing verification of the coverage on a range of possible mu and nu values

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
#import pickle as pickle

import sys
sys.path.insert(0,'../../../undetected_extinctions')    # so I can import the undetected extinctions package
from undetected_extinctions import find_U0_bnd          # find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag, omega, biasedurn)
from undetected_extinctions import SE_changed_E         # collapses timesteps with no detected extinctions

# simulates one instance of a SEUX outcome, but mu_fnc and nu_fnc are now a function of t
def simulate_t(U0, S0, mu_fnc, nu_fnc, T):
    '''
    S, E, U, X = simulate(U0, S0, mu_range, detn_array, T)

    Simulates one possible SEUX outcome. 

    Inputs
    ------

    U0: integer
        Initial number undetected species.
    S0: integer
        Initial number detected species.
    mu_fnc: function
        Accepts t as an input, returns a value of mu, the extinction probability
    nu_fnc: function
        Accepts t as an input, returns a value of nu, the detection probability
    T: integer
        The number of years to simulate

    Returns
    -------

    S, E, U, X: np arrays
        The no. of detected extant, detected extinct, undetected extant, undetected extinct, at each timestep
    '''

    # simulate the population subject to this mu range and detections array
    # ---

    S = [S0]
    E = [0]
    U = [U0]
    X = [0]

    for t in range(1, T+1):

        # how many, of both types, survive?
        mu = mu_fnc(t)
        n = stats.binom( U[t-1]+S[t-1], 1-mu ).rvs()

        # survival process
        Ut = stats.hypergeom( S[t-1] + U[t-1], U[t-1], n ).rvs()  # actually phi
        St = n - Ut                                         # actually psi

        Et = E[t-1] + S[t-1] - St
        Xt = X[t-1] + U[t-1] - Ut

        # detection process
        nu = nu_fnc(t)
        detn = stats.binom( Ut, nu ).rvs()

        St = St+detn
        Ut = Ut-detn

        # append
        S.append(St); U.append(Ut); E.append(Et); X.append(Xt)

    # turn into numpy arrays
    S = np.array(S)
    E = np.array(E)
    U = np.array(U)
    X = np.array(X)

    return S, E, U, X

# simulates one instance of a SEUX outcome
def simulate(U0, S0, mu_fnc, nu_fnc, T):
    '''
    S, E, U, X = simulate(U0, S0, mu_range, detn_array, T)

    Simulates one possible SEUX outcome. 

    Inputs
    ------

    U0: integer
        Initial number undetected species.
    S0: integer
        Initial number detected species.
    mu_fnc: function
        Accepts no inputs, returns a value of mu, the extinction probability
    nu_fnc: function
        Accepts no inputs, returns a value of nu, the detection probability
    T: integer
        The number of years to simulate

    Returns
    -------

    S, E, U, X: np arrays
        The no. of detected extant, detected extinct, undetected extant, undetected extinct, at each timestep
    '''

    # simulate the population subject to this mu range and detections array
    # ---

    S = [S0]
    E = [0]
    U = [U0]
    X = [0]

    for t in range(1, T+1):

        # how many, of both types, survive?
        mu = mu_fnc()
        n = stats.binom( U[t-1]+S[t-1], 1-mu ).rvs()

        # survival process
        Ut = stats.hypergeom( S[t-1] + U[t-1], U[t-1], n ).rvs()  # actually phi
        St = n - Ut                                         # actually psi

        Et = E[t-1] + S[t-1] - St
        Xt = X[t-1] + U[t-1] - Ut

        # detection process
        nu = nu_fnc()
        detn = stats.binom( Ut, nu ).rvs()

        St = St+detn
        Ut = Ut-detn

        # append
        S.append(St); U.append(Ut); E.append(Et); X.append(Xt)

    # turn into numpy arrays
    S = np.array(S)
    E = np.array(E)
    U = np.array(U)
    X = np.array(X)

    return S, E, U, X


def get_coverage_estimates(nsims, nsamples, pcileV, U0, S0, T, mu_fnc, nu_fnc, collapse=False):
    '''
    cnt_withinV, U0_meanV = get_coverage_estimates(params, mu_fnc, nu_fnc)

    Returns counts of true U0 within the confidence intervals, and U0 estimates (50 percentiles)

    Inputs
    ------

    nsims: integer
        Number of simulations to perform
    nsamples: integer 
        Number of samples from which to construct the CIs per simulation
    pcileV: array of integers
        List of percentiles for which to count the coverage
    U0: integer
        Initial number of undetected species
    S0: integer
        Initial number of detected species
    T: integer
        Number of timesteps to simulate
    mu_fnc: function
        Accepts no inputs, returns a value of mu, the extinction probability
    nu_fnc: function
        Accepts no inputs, returns a value of nu, the detection probability

    Returns
    -------
    cnt_withinV: array of integers
        A count of how many simulations had true U0 within the CI; corresponds to pcileV
    U0_meanV: array of floats
        The U0 estimates (50 percentile) corresponding to each of the nsims simulations performed
    '''


    # treatment of percentile depends on if we're doing one or two-sided
    edge_pV = (1-pcileV/100)/2

    # a place to store our return values
    cnt_withinV = np.zeros( edge_pV.shape, dtype=int )
    U0_meanV = list()


    # repeated simulations
    # ---

    for nsim in range(nsims):

        print('doing rep ' + str(nsim) )

        # get the simulation
        S_orig, E_orig, U_orig, X_orig = simulate_t(U0, S0, mu_fnc, nu_fnc, T)

        if collapse:
            # collapse timesteps in which no detected extinctions
            S, E = SE_changed_E(S_orig, E_orig)
        else:
            S = S_orig; E = E_orig

        d = S[1:] - S[:-1] + E[1:] - E[:-1]     # discoveries at each timestep
        psi = S[:-1] - (E[1:]-E[:-1])           # survivors at each timestep
        extns = E[1:] - E[:-1]                  # extinctions at each timestep
        T_idx = len(S)                          # number of timesteps
        tV = list( reversed(range(1,T_idx)) )   # list of timesteps to work backwards through

        # the minimum possible values of U_t are the number detected from t onwards + U_T
        U_T = U_orig[-1]
        min_poss_UV = [ U_T + sum(d[t-1:]) for t in range(1,T_idx) ]


        # repeatedly sample to construct the CI

        U0V = list()

        for nsample in range(nsamples):

            U = [0]*T_idx # a place to store this replicate
            U[-1] = U_T
            impossibleFlag = False

            for t in tV:

                S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]

                # if type_of_ci == 'midp':
                alpha = stats.uniform.rvs()
                min_poss_U0 = max( ( min_poss_UV[t-1], U1+d0 ) )
                U[t-1], impossibleFlag = find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag, None, None)

            U0V.append( U[0] ) # store result


        # append mean

        U0_meanV.append( np.mean(U0V) )


        # count how many within intervals

        # construct the CI using percentiles
        CIV = [ ( np.percentile(U0V, (100-pcile)/2), np.percentile(U0V, 100-(100-pcile)/2) ) for pcile in pcileV ]
        cnt_within = np.array([ 1 if U0 >= U0_lo and U0 <= U0_hi else 0 for U0_lo, U0_hi in CIV ])
        cnt_withinV = cnt_withinV+cnt_within


    return cnt_withinV, U0_meanV

# get an example of a simulation outcome and some sampled bounds
def get_example(U0, S0, mu_fnc, nu_fnc, T, negs):

    # get one simulation
    # ---

    S_orig, E_orig, U_orig, X_orig = simulate_t(U0, S0, mu_fnc, nu_fnc, T)

    S = S_orig; E = E_orig

    d = S[1:] - S[:-1] + E[1:] - E[:-1]     # discoveries at each timestep
    psi = S[:-1] - (E[1:]-E[:-1])           # survivors at each timestep
    extns = E[1:] - E[:-1]                  # extinctions at each timestep
    T_idx = len(S)                          # number of timesteps
    tV = list( reversed(range(1,T_idx)) )   # list of timesteps to work backwards through

    # the minimum possible values of U_t are the number detected from t onwards + U_T
    U_T = U_orig[-1]
    min_poss_UV = [ U_T + sum(d[t-1:]) for t in range(1,T_idx) ]


    # do negs number of examples of sampling bounds back in time
    # ---

    UV = list() # a place to store each of our examples

    for egs in range(negs):

        U = [0]*T_idx # a place to store this replicate
        U[-1] = U_T
        impossibleFlag = False

        for t in tV:

            S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]

            # if type_of_ci == 'midp':
            alpha = stats.uniform.rvs()
            min_poss_U0 = max( ( min_poss_UV[t-1], U1+d0 ) )
            U[t-1], impossibleFlag = find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag, None, None)

        # store this example
        UV.append(U)

    return S_orig, E_orig, U_orig, X_orig, UV



'''
# plot an example of a simulation outcome and some sampled bounds
def plot_example(U0, S0, mu_fnc, nu_fnc, T, negs, fName):

    # get one simulation
    # ---

    S_orig, E_orig, U_orig, X_orig = simulate_t(U0, S0, mu_fnc, nu_fnc, T)

    S = S_orig; E = E_orig

    d = S[1:] - S[:-1] + E[1:] - E[:-1]     # discoveries at each timestep
    psi = S[:-1] - (E[1:]-E[:-1])           # survivors at each timestep
    extns = E[1:] - E[:-1]                  # extinctions at each timestep
    T_idx = len(S)                          # number of timesteps
    tV = list( reversed(range(1,T_idx)) )   # list of timesteps to work backwards through

    # the minimum possible values of U_t are the number detected from t onwards + U_T
    U_T = U_orig[-1]
    min_poss_UV = [ U_T + sum(d[t-1:]) for t in range(1,T_idx) ]


    # do negs number of examples of sampling bounds back in time
    # ---

    UV = list() # a place to store each of our examples

    for egs in range(negs):

        U = [0]*T_idx # a place to store this replicate
        U[-1] = U_T
        impossibleFlag = False

        for t in tV:

            S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]

            # if type_of_ci == 'midp':
            alpha = stats.uniform.rvs()
            min_poss_U0 = max( ( min_poss_UV[t-1], U1+d0 ) )
            U[t-1], impossibleFlag = find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag, None, None)

        # store this example
        UV.append(U)


    # plot them
    # ---


    # plot the simulation as "true" values
    plt.plot(S_orig,  'green',  lw=1, label = r'$S_t$')
    plt.plot(E_orig,  'red',    lw=1, label = r'$E_t$')
    plt.plot(X_orig,  'blue',   lw=1, label = r'$X_t$')


    # plot each of the examples
    for i, U in enumerate(UV):

        if i == 0:
            plt.plot(U, 'orange', lw=0.5, alpha=0.5, label = r'$U_t^{[i]}$')
        else:
            plt.plot(U, 'orange', lw=0.5, alpha=0.5)

    # this one last so it's on top
    plt.plot(U_orig,  'black', lw=2, label = r'$U_t$')

    plt.xlabel('year')
    plt.ylabel('number of species')

    plt.legend(loc='upper right', fontsize='small')
    plt.tight_layout()
    plt.savefig(fName)
    plt.close()
'''
