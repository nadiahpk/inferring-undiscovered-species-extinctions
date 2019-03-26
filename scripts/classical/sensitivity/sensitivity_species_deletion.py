# find the sensitivity of the estimates to random species deletions
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

# list of proportions of species to take
proportions = np.arange(0.3,1,0.1)

tf = 2015
percentile = 95 # CI that I want
omega = None    # central hypergeometric variant
biasedurn = None

# how many random subsets of the full species list we should take per proportion
subsets_per_proportion = 1000

# for each subset, how many samples we should take to estimate the mean
samples_per_subset = 30

# where databases are etc.
# ---

fname_frstlast = '../../../data/processed/first_last_detns_final.csv'
dir_results = '../../../results/classical/sensitivity/sensitivity_species_deletion/'


# get record of first and last observations, and other info
# ---

# read in data
csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)

# each row of frst_last_full corresponds to a particular species
frst_last_full = [ ( int(row[1]), int(row[2]) ) for row in csv_f ]

# need to know how many extant in full record
frst_full, last_full = zip(* frst_last_full )
no_extant_full = last_full.count(tf) # number known extant


# 
no_spp_full = len(frst_last_full)
subset_sizes = [ int(round(no_spp_full*propn)) for propn in proportions ]

idxs = list(range(no_spp_full))

extn_rateM = list()
for subset_size in subset_sizes:

    print(subset_size)
    extn_rate_samples = list()

    for sample in range(samples_per_subset):

        # take the random subset of the full species list, and treat the same as we did for the full list
        # ---

        # get random subset
        subset_idxs = np.random.choice( idxs, subset_size, replace=False )
        frst_last = [ frst_last_full[i] for i in subset_idxs ]

        # how many undetected extant species (i.e. species known to be extant that are missing from the subset
        frst, last = zip(* frst_last )
        no_extant_subset = last.count(tf)
        U_T = no_extant_full - no_extant_subset

        # modify record, so that we only take intervals where the number of detected extinct species E changes
        years_mod, frst_last_mod = frst_last_changed_E(frst_last)

        # calculate the S and E for the timeseries
        _, S, E = get_SE(frst_last_mod, years_mod)


        # get important information about the record
        d = S[1:] - S[:-1] + E[1:] - E[:-1]     # discoveries at each timestep
        psi = S[:-1] - (E[1:]-E[:-1])           # survivors at each timestep
        extns = E[1:] - E[:-1]                  # extinctions at each timestep
        T_idx = len(S)                          # number of timesteps


        # obtain a central estimate of extinction rate through repeated sampling of confidence level
        # ---

        extn_rateV = list()

        tV = list( reversed(range(1,T_idx)) ) # list of timesteps to work backwards through

        for nrep in range(samples_per_subset):

            # work our way backwards through the time series, randomly sampling confidence leve

            U = [0]*T_idx # a place to store this replicate
            U[-1] = U_T
            impossibleFlag = False

            for t in tV: # work our way backwards

                alpha = uniform.rvs()
                S0 = S[t-1]; S1 = S[t]; U1 = U[t]; d0 = d[t-1]
                U[t-1], impossibleFlag = find_U0_bnd(alpha, S0, S1, U1, d0, impossibleFlag, omega, biasedurn)


            # calculate summary info
            N = S[0] + E[0] + U[0]              # assumes X(0) = 0
            extn_rate = (N-S[-1]-U[-1]) / N     # calculate extinction rate

            # append results
            extn_rateV.append(extn_rate)

        extn_rate_samples.append( np.mean(extn_rate) )

    # append
    extn_rateM.append( extn_rate_samples )

    # store extn_rate_samples just in case
    # ---

    f = open(dir_results + 'sensitivity_spp_deletion_' + str(subset_size) + '.pkl', 'wb')
    pickle.dump( extn_rate_samples, f )
    f.close()


extn_rateM = np.array(extn_rateM)

# store final result
# ---

f = open(dir_results + 'sensitivity_spp_deletion_all.pkl', 'wb')
# a string explaining the pickle file
pickle.dump( extn_rateM, f )
f.close()


# get the full result to append its point to the graph
# ---

fname_classical = '../../../results/classical/classical_basic_result.csv'
csv_f = csv.reader( open(fname_classical) )
header = next(csv_f)
c_res = [ [ float(ri) for ri in row ] for row in csv_f ]
years_mod, S, E, U, X, _, _, _, _, = zip(* c_res )

N = S[0] + E[0] + U[0]              # assumes X(0) = 0
extn_rate = (N-S[-1]-U[-1]) / N     # calculate extinction rate

# plot
# ---

mean = np.mean(extn_rateM, axis=1)
lo = np.percentile(extn_rateM, (100-percentile)/2, axis=1)
hi = np.percentile(extn_rateM, 100 - (100-percentile)/2, axis=1)

mean2 = np.append( mean, [extn_rate] )
lo2 = np.append( lo, [extn_rate] )
hi2 = np.append( hi, [extn_rate] )
proportions2 = np.append( proportions, [1] )

plt.plot(proportions2, mean2, color="black")
plt.fill_between(proportions2, lo2, hi2, facecolor='black', alpha=0.5)
plt.xlabel('proportion of species included')
plt.ylabel('estimated extinction rate')
plt.ylim( (0.2, 0.8) )
plt.grid(True)
# plt.show()
plt.tight_layout()
plt.savefig(dir_results + 'sensitivity_species_deletion.pdf')
plt.close()
