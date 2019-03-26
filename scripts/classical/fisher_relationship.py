# get the relationship between the omega and the total number of species infered

import sys
sys.path.insert(0,'../../undetected_extinctions') # so I can import the undetected extinctions package

import csv
import numpy as np
from scipy.stats import uniform
import pickle

from rpy2.robjects.packages import importr # so I can import R package BiasedUrn for the Fisher variant
import rpy2.robjects.numpy2ri

from undetected_extinctions import frst_last_changed_E, get_SE, find_U0_bnd


# user parameters
# ---

percentile = 95
# no samples; omega values to explore; file name to plot to
nreps = 10000
omegaV = np.linspace(.15, .95, 17)

U_T = 0         # assume that at the final timestep there are no undetected species remaining

rpy2.robjects.numpy2ri.activate()
biasedurn = importr('BiasedUrn')

# where databases are etc.
# ---

fname_frstlast = '../../data/processed/first_last_detns_final.csv'
fname_results = '../../results/classical/fisher_relationship/fisher_relationship.csv' # where we append our data to


# get record of first and last observations, and other info
# ---

# read in data
csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)

# each row of frst_last_full corresponds to a particular species
frst_last = [ ( int(row[1]), int(row[2]) ) for row in csv_f ]

# modify record, so that we only take intervals where the number of detected extinct species E changes
years_mod, frst_last_mod = frst_last_changed_E(frst_last)

# calculate the S and E for the timeseries
_, S, E = get_SE(frst_last_mod, years_mod)

# get important information about the record
d = S[1:] - S[:-1] + E[1:] - E[:-1]     # discoveries at each timestep
psi = S[:-1] - (E[1:]-E[:-1])           # survivors at each timestep
extns = E[1:] - E[:-1]                  # extinctions at each timestep
T_idx = len(S)                          # number of timesteps

tV = list( reversed(range(1,T_idx)) ) # list of timesteps to work backwards through


# for each omega in the list, estimate N, and append to results file
# ---


for omega in omegaV:

    print(omega)

    # sampling to obtain estimates of total no. of species N
    # ---

    NV = list()
    for nrep in range(nreps):

        # work our way backwards through the time series, randomly sampling confidence leve

        if nrep % 10 == 0:
            print('omega = ', str(omega), '; nrep = ', str(nrep) )

        U = [0]*T_idx # a place to store this replicate
        U[-1] = U_T
        impossibleFlag = False

        for t in tV: # work our way backwards

            alpha = uniform.rvs()
            S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]
            U[t-1], impossibleFlag = find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag, omega, biasedurn)

        # calculate summary info and store
        N = S[0] + E[0] + U[0]              # assumes X(0) = 0
        NV.append(N)


    # calculate summary statistics for this value of omega
    # ---

    N_mean = np.mean(NV)
    N_lo = np.percentile(NV, (100-percentile)/2)
    N_hi = np.percentile(NV, 100 - (100-percentile)/2)


    # append to file
    # ---

    f = open(fname_results, 'a') # NOTE 'a', appending, so it can be done piecemeal

    # header order is: omega,N_mean,N_lo,N_hi,nreps
    f.write( ','.join( list(map(str, [omega, N_mean, N_lo, N_hi, nreps] )) ) )
    f.write('\n')

    f.close()

