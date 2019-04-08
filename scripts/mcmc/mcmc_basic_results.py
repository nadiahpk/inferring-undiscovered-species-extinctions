# obtain basic results from the two chains, the 95-percentiles

import sys
sys.path.insert(0,'../../undetected_extinctions') # allows us to import undetected extinctions package

import pickle
import numpy as np
import matplotlib.pyplot as plt
import csv

from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate() # so I can pass variables to it

from undetected_extinctions import frst_last_changed_E, get_SE


# parameters
# ---

suffixes = [0, 1, 2, 3, 4, 5] # the suffixes used for the file names
burnin = 15000 # discard the first 15000 iterations as burn in
percentile = 95 # CI that I want


# location of chains and info
# ---

dir_results = '../../results/mcmc/'
fname_years = '../../results/mcmc/years_mod.pkl'
fname_frstlast = '../../data/processed/first_last_detns_final.csv'


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


# obtain the chains, minus the burn-in, and prune them modestly
# ---

# chain 1
f = open(dir_results + 'mcmc_lo_0.pkl', 'rb');
UV1 = pickle.load( f );
f.close()
UV1 = UV1[burnin:,:] # remove burnin
for suffix in suffixes[1:]:
    f = open(dir_results + 'mcmc_lo_' + str(suffix) + '.pkl', 'rb');
    UVnext = pickle.load( f ); f.close()
    UV1 = np.append(UV1, UVnext, axis=0)

# chain 2
f = open(dir_results + 'mcmc_hi_0.pkl', 'rb');
UV2 = pickle.load( f );
f.close()
UV2 = UV2[burnin:,:] # remove burnin
for suffix in suffixes[1:]:
    f = open(dir_results + 'mcmc_hi_' + str(suffix) + '.pkl', 'rb');
    UVnext = pickle.load( f ); f.close()
    UV2 = np.append(UV2, UVnext, axis=0)

# combine UV
UV = np.append(UV1, UV2, axis=0)

# prune them to match what we used for the MCMC SE calculation
k = 100
UV = UV[::k,:]


# calculate X_t using U_t
# ---

N = S[0]+E[0]+UV[:,0] # assumes X(0) = 0

XV = list()
for i in range(len(N)):

    X = N[i]-E-S-UV[i,:]
    XV.append(X)

XV = np.array(XV)


# calculate statistics on sample
# ---

X_mean = np.mean(XV, axis=0)
U_mean = np.mean(UV, axis=0)

X_lo = np.percentile(XV, (100-percentile)/2,       axis=0)
X_hi = np.percentile(XV, 100 - (100-percentile)/2, axis=0)

U_lo = np.percentile(UV, (100-percentile)/2,       axis=0)
U_hi = np.percentile(UV, 100 - (100-percentile)/2, axis=0)


# store statistics to csv file
# ---

f = open(dir_results + 'mcmc_basic_result.csv', 'w')
f.write('year,S,E,U_mean,X_mean,U_lo,U_hi,X_lo,X_hi\n')

for row in zip(years_mod,S, E, U_mean, X_mean, U_lo, U_hi, X_lo, X_hi):
    row_string = list(map( lambda v: str(v), row ))
    f.write(','.join(row_string))
    f.write('\n')

f.close()
