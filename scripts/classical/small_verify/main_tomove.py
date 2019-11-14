if __name__ == "__main__":


    # user set parameter values
    # ---

    nreps = 1000  # number of repetitions
    nsamples = 1000 # number of samples from which to construct the CIs

    pcileV = np.array([90, 80, 70, 60, 50, 40, 30, 20, 10]) # list of percentiles to do

    U0 = 120    # initial number of undetected species
    S0 = 20     # initial number of detected species

    mu_range = (0.01, 0.07) # probability of going extinct each timestep
    detn_array = [1,2,3] # no. species detected each year
    T = 30 # length of timesteps to simulate (ends near U_T = 0)

    type_of_ci = 'midp'

    # store parameters all together
    params = {
            'nsamples': nsamples,
            'nreps': nreps,
            'pcileV': pcileV,
            'U0': U0,
            'S0': S0,
            'mu_range': mu_range,
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
    d = None

    for nrep in range(nreps):

        print('doing rep ' + str(nrep) )


        # get the simulation (we can reuse the n we obtain)

        S_orig, E_orig, U_orig, X_orig = simulate3(U0, S0, mu_range, detn_array, T)

        # S, E = SE_changed_E(S_orig, E_orig)     # NOTE: don't collapse timesteps in which there were no detected extinctions
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

                alpha = uniform.rvs()
                S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]

                if type_of_ci == 'midp':

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


    if False:

        # plot coverage for each confidence level
        # ---


        coverage = 100*cnt_withinV/nreps

        plt.plot( pcileV, pcileV, ls='dotted', color='black')
        plt.scatter( pcileV, coverage, color='black')
        plt.xlabel('nominal coverage desired')
        plt.ylabel('actual coverage obtained')
        plt.grid(True)
        plt.show()

    else:


        # save summary results to a pickle file
        # ---

        fWriteName = 'verify_' + type_of_ci + '_U0_'  + str(U0) + '_nreps_' + str(nreps) + '_nsamples_' + str(nsamples) + '.pkl'
        f = open('../../../results/classical/small_verify/' + fWriteName, 'wb')

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

