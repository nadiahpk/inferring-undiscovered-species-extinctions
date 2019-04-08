# obtain summary results from the two chains, and the MCMC standard errors on the quantile estimates

import pickle
import numpy as np
import matplotlib.pyplot as plt

from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate() # so I can pass variables to it


# parameters
# ---

suffixes = [0, 1, 2, 3, 4, 5] # the suffixes used for the file names
burnin = 15000 # discard the first 15000 iterations as burn in NOTE -- in spare.pdf

# location of chains and info
# ---

dir_results = '../../results/mcmc/'
fname_years = '../../results/mcmc/years_mod.pkl'
#fname_chain1 = '../../results/mcmc/mcmc_lo_0.pkl'
#fname_chain2 = '../../results/mcmc/mcmc_hi_0.pkl'

mcmcse = importr('mcmcse') # import the R function

# obtain the chains and info, calculated ESS as we go
# ---

# get the years
f = open(fname_years, 'rb');
years_mod = pickle.load( f );
f.close()

# chain 1

f = open(dir_results + 'mcmc_lo_0.pkl', 'rb');
UV1 = pickle.load( f );
f.close()

# calculate ESS
ess1V = np.array( mcmcse.ess(UV1[:,:-1]) )

for suffix in suffixes[1:]:

    f = open(dir_results + 'mcmc_lo_' + str(suffix) + '.pkl', 'rb');
    UVnext = pickle.load( f ); f.close()

    essnextV = np.array( mcmcse.ess(UVnext[:,:-1]) )
    ess1V = ess1V + essnextV

    UV1 = np.append(UV1, UVnext, axis=0)

UV1 = UV1[burnin:,:] # remove burnin

# chain 2

f = open(dir_results + 'mcmc_hi_0.pkl', 'rb');
UV2 = pickle.load( f );
f.close()

# calculate ESS
ess2V = np.array( mcmcse.ess(UV2[:,:-1]) )

for suffix in suffixes[1:]:

    f = open(dir_results + 'mcmc_hi_' + str(suffix) + '.pkl', 'rb');
    UVnext = pickle.load( f ); f.close()

    essnextV = np.array( mcmcse.ess(UVnext[:,:-1]) )
    ess2V = ess2V + essnextV

    UV2 = np.append(UV2, UVnext, axis=0)

UV2 = UV2[burnin:,:] # remove burnin

# omit last variable U_T = 0

UV1 = UV1[:,:-1]
UV2 = UV2[:,:-1]

no_sams, no_vars = UV1.shape # number of samples, number of variables
no_sams = no_sams*2 # two chains

print('got up to here 1')

# find statistics on MCMC samples and plot
# ---

# effective sample size of chains
essV = ess1V+ess2V

print('got up to here 2')

# the minimum effective sample size needed for stable analysis
minESS = mcmcse.minESS(p = no_vars, alpha = .05, eps = .1)[0] # -- 1997
# mcmcse.minESS(p = 108, alpha = .05, eps = .05)[0] # -- 7987
# for theoretical details see: 
#   Vats, D., Flegal, J. M., and Jones, G. L. (2017a). Multivariate output analysis for Markov chain Monte Carlo. 
#   arXiv preprint arXiv:1512.07713.

# multiESS = mcmcse.multiESS( UV1 )[0] + mcmcse.multiESS( UV2 )[0] 
# -- something wrong with this, gave 91537 when ESS for variables separately ranged from only 807 to 14567 (first 150,000 samples)
# --- was I meant to divide it by the number of variables?

# plot to check that our sample sizes exceed the minimum needed
plt.scatter(years_mod[:-1], essV, color='black', label='MCMC results from two chains totalling ' + str(no_sams) + ' samples')
plt.xlabel('year $t$')
plt.ylabel(r'effective sample size for $U_t$')
plt.axhline(minESS, ls='dashed', color='black', label='minimum needed for 95 percentile at tol 0.1')
plt.grid(True)
plt.legend(loc='upper left')
#plt.show()
plt.ylim( (0,10000) )
plt.tight_layout()
plt.savefig(dir_results + 'mcmc_results_checkESS.pdf')
plt.close()

print('got up to here 3')

# credible intervals for selected U_t along with the MCMC Error for estimate and interval bounds
# ---

tV = [0, 32]

for t in tV:

    # bins set up so we have 1 for each integer
    min_data1 = min(UV1[:,t]); max_data1 = max(UV1[:,t]);
    min_data2 = min(UV2[:,t]); max_data2 = max(UV2[:,t]);

    bins1 = [ i-.5 for i in range(min_data1, max_data1+2) ]
    bins2 = [ i-.5 for i in range(min_data2, max_data2+2) ]

    # plot the comparison

    plt.figure(figsize=(10,5))

    # chain 1
    plt.hist( UV1[:,t], alpha=.5, bins=bins1, normed=True, color='red', label='chain 1') #, edgecolor='black' )

    # chain 2
    plt.hist( UV2[:,t], alpha=.5, bins=bins2, normed=True, color='blue', label='chain 2') #, edgecolor='black' )

    #plt.xlim( (2150,2650) ) # 
    plt.legend(loc='upper right')
    plt.xlabel('$U_{' + str(years_mod[t]) + '}$')
    plt.ylabel('normalised frequency')
    # plt.show()
    plt.tight_layout()
    plt.savefig(dir_results + 'mcmc_results_histogram_' + str(years_mod[t]) + '.pdf')
    plt.close()


    
# obtain the combined results and write to a csv
# ---

# combine UV
UV = np.append(UV1, UV2, axis=0)

# the function below runs out of memory if UV is too long, so shorten it by taking every k'th value
k = 100
UV_shorten = UV[::k,:]

# find estimate and credible interval for combined sample

vv = mcmcse.mcse_mat(x = UV_shorten, method="bm") # comes out as estimates then errors
est = np.array(vv[:no_vars])
est_se = np.array(vv[-no_vars:])

vv = mcmcse.mcse_q_mat(x=UV_shorten, q=0.025)
lo = np.array(vv[:no_vars])
lo_se = np.array(vv[-no_vars:])

vv = mcmcse.mcse_q_mat(x=UV_shorten, q=0.975)
hi = np.array(vv[:no_vars])
hi_se = np.array(vv[-no_vars:])

# write to csv

f = open(dir_results + 'mcmc_se_results.csv', 'w')

f.write('year,U_est,U_est SE,U_lo,U_lo SE,U_hi,U_hi SE\n')

for row in zip(years_mod, est, est_se, lo, lo_se, hi, hi_se):

    row_string = list(map( lambda v: str(v), row ))
    f.write(','.join(row_string))
    f.write('\n')

f.close()
