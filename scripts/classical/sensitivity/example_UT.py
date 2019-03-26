# Run a full example of SEUX outcome for the scenario when U_T is actually high


import sys
sys.path.insert(0,'../../../undetected_extinctions') # so I can import the undetected extinctions package

import csv
import numpy as np
from scipy.stats import uniform
import matplotlib.pyplot as plt
import pickle

from undetected_extinctions import frst_last_changed_E, get_SE, find_U0_bnd


# user parameters
# ---

U_T = 500 # choose the largest

tf = 2015
percentile = 95 # CI that I want
omega = None    # central hypergeometric variant
biasedurn = None
nreps = 10000 # how many samples to take to estimate statistics


# where databases are etc.
# ---

fname_frstlast = '../../../data/processed/first_last_detns_final.csv'
dir_results = '../../../results/classical/sensitivity/sensitivity_UT/'


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
# ---

d = S[1:] - S[:-1] + E[1:] - E[:-1]     # discoveries at each timestep
psi = S[:-1] - (E[1:]-E[:-1])           # survivors at each timestep
extns = E[1:] - E[:-1]                  # extinctions at each timestep
T_idx = len(S)                          # number of timesteps

tV = list( reversed(range(1,T_idx)) ) # list of timesteps to work backwards through


# repeatedly sample
# ---

UM = list() # a place to store our replicates
XM = list() # a place to store our replicates

for nrep in range(nreps):

    # print replicate number every 10th

    if nrep % 10 == 0:
        print(nrep)

    # work our way backwards through the time series, randomly sampling confidence leve

    U = [0]*T_idx # a place to store this replicate
    U[-1] = U_T
    impossibleFlag = False

    for t in tV: # work our way backwards

        alpha = uniform.rvs()
        S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]
        U[t-1], impossibleFlag = find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag, omega, biasedurn)


    # calculate X_t from U_t

    N = S[0]+E[0]+U[0] # assumes X(0) = 0
    X = [ N - EE - SS - UU for EE, SS, UU in zip(E,S,U) ]

    # append results

    UM.append(U)
    XM.append(X)


# calculate statistics on sample
# ---

UM = np.array(UM)
XM = np.array(XM)

X_mean = np.mean(XM, axis=0)
U_mean = np.mean(UM, axis=0)

X_lo = np.percentile(XM, (100-percentile)/2,       axis=0)
X_hi = np.percentile(XM, 100 - (100-percentile)/2, axis=0)

U_lo = np.percentile(UM, (100-percentile)/2,       axis=0)
U_hi = np.percentile(UM, 100 - (100-percentile)/2, axis=0)


# store statistics to csv file
# ---

f = open(dir_results + 'example_UT.csv', 'w')
f.write('year,S,E,U_mean,X_mean,U_lo,U_hi,X_lo,X_hi\n')

for row in zip(years_mod,S, E, U_mean, X_mean, U_lo, U_hi, X_lo, X_hi):
    row_string = list(map( lambda v: str(v), row ))
    f.write(','.join(row_string))
    f.write('\n')

f.close()
