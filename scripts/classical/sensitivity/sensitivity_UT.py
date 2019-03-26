# find the sensitivity of the estimates to U_T > 0
# create a plot of proportion of species included in the record versus inferred proportion extinctions


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

# list of U_T to take
UTV = [0, 100, 200, 300, 500] # NOTE

tf = 2015
percentile = 95 # CI that I want
omega = None    # central hypergeometric variant
biasedurn = None
samples_per_UT = 1000 # how many samples to take to estimate statistics


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
d = S[1:] - S[:-1] + E[1:] - E[:-1]     # discoveries at each timestep
psi = S[:-1] - (E[1:]-E[:-1])           # survivors at each timestep
extns = E[1:] - E[:-1]                  # extinctions at each timestep
T_idx = len(S)                          # number of timesteps

tV = list( reversed(range(1,T_idx)) ) # list of timesteps to work backwards through

extn_rateM = list()
for U_T in UTV:

    print(U_T)

    extn_rateV = list()
    for sample in range(samples_per_UT):

        # work our way backwards through the time series, randomly sampling confidence leve

        U = [0]*T_idx # a place to store this replicate
        U[-1] = U_T
        impossibleFlag = False

        for t in tV: # work our way backwards

            alpha = uniform.rvs()
            S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]
            U[t-1], impossibleFlag = find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag, omega, biasedurn)

        # calculate summary info and store
        N = S[0] + E[0] + U[0]              # assumes X(0) = 0
        extn_rate = (N-S[-1]-U[-1]) / N     # calculate extinction rate
        extn_rateV.append(extn_rate)        # append results

    extn_rateM.append(extn_rateV)

extn_rateM = np.array(extn_rateM)


# store final result
# ---

f = open(dir_results + 'sensitivity_UT.pkl', 'wb')
pickle.dump( extn_rateM, f )
f.close()


# plot
# ---

mean = np.mean(extn_rateM, axis=1)
lo = np.percentile(extn_rateM, (100-percentile)/2, axis=1)
hi = np.percentile(extn_rateM, 100 - (100-percentile)/2, axis=1)

plt.plot(UTV, mean, color="black")
plt.fill_between(UTV, lo, hi, facecolor='black', alpha=0.5)
plt.xlabel('assumed undetected extant at end of record $U_{2015}$')
plt.ylabel('estimated extinction rate')
plt.ylim( (0.2, 0.8) )
plt.grid(True)
# plt.show()
plt.tight_layout()
plt.savefig(dir_results + 'sensitivity_UT.pdf')
plt.close()
