import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

from small_verify_fncs import simulate, get_coverage_estimates

import sys
sys.path.insert(0,'../../../undetected_extinctions') # so I can import the undetected extinctions package
from undetected_extinctions import find_U0_bnd

#S, E, U, X = simulate(50, 10, lambda: 0.1, lambda: 0.1, 10)

# user set params
params = {
        'name_run': 'test',
        'nsamples': 10,
        'nsims': 10,
        'pcileV': np.array([90, 80, 70, 60, 50, 40, 30, 20, 10]),
        'U0': 50,
        'S0': 10,
        'T': 10
        }

mu_fnc = lambda: 0.1; nu_fnc = lambda: 0.1

nsamples = params['nsamples']   # number of samples from which to construct the CIs
nsims = params['nsims']         # number of repetitions
pcileV = params['pcileV']       # list of percentiles to check coverage for
U0 = params['U0']               # initial number of undetected species
S0 = params['S0']               # initial number of detected species
T = params['T']                 # length of timesteps to simulate (ends near U_T = 0)

# NOTE this is the main bit
# cnt_withinV, U0_meanV = get_coverage_estimates(nsims, nsamples, pcileV, U0, S0, T, mu_fnc, nu_fnc)

# NOTE this is in plot_tomove.py, I need to make a proper function of it

# NOTE this is a function that runs a quick example for us to plot
if True:

    # parameters
    # ---

    negs = 30 # how many example trajectories do we want to do?
    

    # get one simulation
    # ---

    S_orig, E_orig, U_orig, X_orig = simulate(U0, S0, mu_fnc, nu_fnc, T)

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
    #plt.legend(loc='best')
    plt.show()


'''

    if False:

        # plot coverage for each confidence level
        # ---


        coverage = 100*cnt_withinV/nsims

        plt.plot( pcileV, pcileV, ls='dotted', color='black')
        plt.scatter( pcileV, coverage, color='black')
        plt.xlabel('nominal coverage desired')
        plt.ylabel('actual coverage obtained')
        plt.grid(True)
        plt.show()

    else:


        # save summary results to a pickle file
        # ---

        fWriteName = 'verify_' + name_run + '_U0_'  + str(U0) + '_nsims' + str(nsims) + '_nsamples_' + str(nsamples) + '.pkl'
        f = open('../../../results/classical/small_verify/' + fWriteName, 'wb')

        # a string explaining the pickle file
        ss  = 'Created by verify.py\n'
        ss += 'Contains the following:\n'
        ss += ' - ss, str: this string you are reading now.\n'
        ss += ' - params, dict: parameters used for the run.\n'
        ss += ' - cnt_withinV, np array of ints: the count of how many times the true U0 was within the bounds corresponding to pcileV.\n'
        ss += ' - U0_meanV, list of ints: the mean obtained for each rep (rounded).\n'
        pickle.dump( ss, f ) # 0.
        pickle.dump( params, f )
        pickle.dump( cnt_withinV, f )
        pickle.dump( U0_meanV, f )
        f.close()
'''

