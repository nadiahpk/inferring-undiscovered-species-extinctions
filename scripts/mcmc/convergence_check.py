# check convergence of the two chains

import pickle
import numpy as np
import matplotlib.pyplot as plt
import pymc3 as pymc3


# parameters
# ---

suffixes = [0, 1] # results suffixes we'll use
sample_max = 40000 # how many iterations to do the plot for
end_bits = 2000 # start and end gap

suffixes = [0, 1] # NOTE test it works
sample_max = 100000
end_bits = 1000 # start and end gap

# location of chains and info
# ---

dir_results = '../../results/mcmc/'
fname_years = '../../results/mcmc/years_mod.pkl'
#fname_chain1 = '../../results/mcmc/mcmc_lo_0.pkl'
#fname_chain2 = '../../results/mcmc/mcmc_hi_0.pkl'


# obtain the chains and info
# ---

# get the years
f = open(fname_years, 'rb');
years_mod = pickle.load( f );
f.close()

# chain 1

f = open(dir_results + 'mcmc_lo_0.pkl', 'rb');
UV1_ = pickle.load( f );
f.close()

for suffix in suffixes[1:]:

    f = open(dir_results + 'mcmc_lo_' + str(suffix) + '.pkl', 'rb');
    UVnext = pickle.load( f ); f.close()
    UV1_ = np.append(UV1_, UVnext, axis=0)

# chain 2

f = open(dir_results + 'mcmc_hi_0.pkl', 'rb');
UV2_ = pickle.load( f );
f.close()

for suffix in suffixes[1:]:

    f = open(dir_results + 'mcmc_hi_' + str(suffix) + '.pkl', 'rb');
    UVnext = pickle.load( f ); f.close()
    UV2_ = np.append(UV2_, UVnext, axis=0)


# trim to the length of the plot
# ---

UV1 = UV1_[ :sample_max, : ]
UV2 = UV2_[ :sample_max, : ]

# find the last variable to satisfy $\hat{R} < 1.1$ and plot
# ---

GR_threshold = 1.1

# take log to correct for deviation from normal
log_UV1 = np.log(UV1); log_UV2 = np.log(UV2)

# sections to calculate GR on
lenUV = UV1.shape[0]
llV = list(range(-lenUV+end_bits,-end_bits))[0::1000] # from end_bits to end_bits from the end, steps of 1000
xV = [ lenUV+ll for ll in llV] # list of iteration numbers

# find the last to cross

last_x = 0 # initialise
last_i = 0
last_grV = list()
for i in range(len(years_mod)-1): # we don't do the last variable U_T bc it's always zero

    # calculate GR on sections

    grV = [ pymc3.gelman_rubin( np.vstack( (log_UV1[-ll:,i] , log_UV2[-ll:,i]) ) ) for ll in llV ]

    # if the Gelman-Rubin statistic calls below the threshold later than the latest we've found so far
    # then store

    threshold_x = next( x for x, gr in zip(xV, grV) if gr <= GR_threshold)

    if threshold_x > last_x:

        last_x = threshold_x
        last_i = i
        last_grV = grV
        


# plot the traces and the Gelman-Rubin statistic for the variable that was last to satisfy $\hat{R}$
# ---

f, axarr = plt.subplots(nrows=2, ncols=1, sharex=True)

# plot traces
axarr[0].plot(UV1[:,last_i], color='red', lw=0.5, label='chain 1')
axarr[0].plot(UV2[:,last_i], color='blue', lw=0.5, label='chain 2')
axarr[0].set_ylabel(r'$U_{' + str(years_mod[last_i]) + '}$ trace')
axarr[0].grid(True)
axarr[0].legend(loc='upper right')

# plot GR statistic
axarr[1].plot(xV, last_grV, lw=2, color='black', label=r'$U_{' + str(years_mod[last_i]) + '}$ trace')
axarr[1].text(last_x+2000, GR_threshold+0.01, 'iteration ' + str(last_x))
axarr[1].set_xlabel(r'iteration')
axarr[1].axhline(GR_threshold, ls='dotted', color='black', label='threshold')
axarr[1].axvline(last_x, ls='dotted', color='black')
axarr[1].set_ylabel(r'Gelman-Rubin statistic $\hat{R}$')
axarr[1].grid(True)
axarr[1].legend(loc='upper right')

#plt.show()
plt.tight_layout()
plt.savefig('../../results/mcmc/convergence_check.pdf')
plt.close()
