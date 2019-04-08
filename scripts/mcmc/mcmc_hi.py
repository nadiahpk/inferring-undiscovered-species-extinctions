# run an mcmc chain, or continue a previous one

import sys
sys.path.insert(0,'../..') # allows us to import undetected extinctions package

import csv
import numpy as np
from scipy.stats import hypergeom, uniform
import random as random
import pickle

from undetected_extinctions.undetected_extinctions import frst_last_changed_E, get_SE


# parameters
# ---

U_T = 0 # assume that at the final timestep there are no undetected species remaining
nreps = 150000
start_type = 'hi' # or can be: 'hi', to high-start chains; 'lo', for low-start chains
previous_suffix = 4 # None, or set to a number to start at a previous point

'''
# for the first run, suffix = 0
previous_suffix = None # None, or set to a number to start at a previous point
nreps = 30000
'''

# where databases are
# ---

fname_frstlast = '../../data/processed/first_last_detns_final.csv'
dir_results = '../../results/mcmc/'


# get record of first and last observations
# ---

# read in data
csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
frst_last = [ ( int(row[1]), int(row[2]) ) for row in csv_f ]

# modify record, so that we only take intervals where the number of detected extinct species E changes
years_mod, frst_last_mod = frst_last_changed_E(frst_last)

# calculate the S and E for the timeseries
_, S, E = get_SE(frst_last_mod, years_mod)


# get important information about the record
# ---

d = S[1:] - S[:-1] + E[1:] - E[:-1]     # discoveries at each timestep
psi = S[:-1] - (E[1:]-E[:-1])           # survivors at each timestep
extns = E[1:] - E[:-1]                  # extinctions at each timestep

# the minimum possible no of undetected species at each timestep
U_min = np.array(list(reversed(np.cumsum(list(reversed( list(d) + [U_T] ))))))


# define the sampling functions I'll be reusing
# ---

# P(phi_t | U_t)
fnc_Ptt = lambda phi_t, t: hypergeom.pmf( phi_t, U[t]+S[t], U[t], phi_t + psi[t] )

# P(phi_{t+1} | phi_t) a likelihood
fnc_P1t = lambda phi_t, t: hypergeom.pmf( phi[t+1], phi_t-d[t]+S[t+1], phi_t-d[t], phi[t+1]+psi[t+1])


# initialise
# ---

if previous_suffix is None: # we are starting a new chain

    suffix = 0

    if start_type == 'lo':

        # starting this chain with the minimum possible value of U
        U = U_min 

    else:

        # starting this chain at a high value of U
        U = [ int(UU*1.5) for UU in U_min ]

else:

    # open the file where our previous results are stored
    fname_cont = dir_results + 'mcmc_' + start_type + '_' + str(previous_suffix) + '.pkl'
    f = open(fname_cont, 'rb')
    UV =pickle.load(f)
    f.close()

    U = list(UV[-1,:])
    suffix = previous_suffix + 1

phi = U[1:]+d
T = len(years_mod) - 1


# do nreps scans
# ---

UV = list()

for rep in range(nreps):

    if rep % 1000 == 0:
        print(str(rep))

    # update U_0 using the backwards sampling
    # ---

    # call u_0 = U_0 − phi_0 and sample u_0

    norm_L = S[0] * ( phi[0]+psi[0] ) / ( psi[0]*(psi[0]-1) )

    r = uniform.rvs() # choose our random number
    cum_normed = 0 # initialise
    u = 0
    while cum_normed < r:

        # step forward
        L = hypergeom.pmf( phi[0], S[0] + phi[0] + u, phi[0] + u, S[1] + U[1] )
        cum_normed += L/norm_L
        u += 1

    u -= 1
    U[0] = phi[0]+u


    # update each subsequent phi_0... and U_1... using the conditionals
    # ---

    for t in range(T-1):

        # get my old phi_t and calculate its likelihood
        phi_old = phi[t]
        L_old = fnc_Ptt(phi_old,t) * fnc_P1t(phi_old,t)

        # randomly choose a new phi_t from proposal distn and get likelihood
        phi_new = random.randint( phi[t+1]+d[t], U[t] )
        L_new = fnc_Ptt(phi_new,t) * fnc_P1t(phi_new,t)

        # accept if acceptance probability greater or equal to random number
        if L_new / L_old >= random.random():

            # update storage in prep for next gibbs update
            phi[t] = phi_new
            U[t+1] = phi_new - d[t]

        # otherwise we just keep the old value


    # store our U
    UV.append( [ UU for UU in U ] )

UV = np.array(UV)


# save to file
# ---

f = open(dir_results + 'mcmc_' + start_type + '_' + str(suffix) + '.pkl', 'wb')
pickle.dump( UV, f ) # 0.
f.close()
