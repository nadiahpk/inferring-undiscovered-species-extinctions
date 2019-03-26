# a script for plotting the coverage of various methods (mid-P, clopper) of obtaining confidence intervals by sampling

import pickle as pickle
import matplotlib.pyplot as plt
import numpy as np

# user specify filename
# ---

# compare Tak and my methods
# fName = 'verify_tak1_U0_150_nreps_50_nsamples_100.pkl'
# fName = 'verify_midp_U0_150_nreps_50_nsamples_100.pkl'

# explore which part, reps or samples, we need to improve to get good fits
# fName = 'verify_midp_U0_150_nreps_50_nsamples_200.pkl' # this looks good, wobbles around the centre
# fName = 'verify_midp_U0_150_nreps_100_nsamples_100.pkl' # this wobbles less and also looks quite good

# compare again to Tak's method
# fName = 'verify_tak1_U0_150_nreps_100_nsamples_100.pkl' # looks terrible
# -- cf to 'verify_midp_U0_150_nreps_100_nsamples_100.pkl'

# compare my and Tak's methods for low U0
# fName = 'verify_tak1_U0_50_nreps_100_nsamples_100.pkl'
# fName = 'verify_midp_U0_50_nreps_100_nsamples_100.pkl' # mine's still good, but Tak's less bad than before

# compare my and Tak's methods for even lower U0
# fName = 'verify_midp_U0_20_nreps_100_nsamples_100.pkl'
# fName = 'verify_tak1_U0_20_nreps_100_nsamples_100.pkl'

# check how it goes when detn varies
# fName = 'verify_midp_U0_150_nreps_100_nsamples_100.pkl' # mu=0.1, d=2 constant
# fName = 'verify_midp_U0_120_nreps_100_nsamples_100.pkl' # mu=0.05, d=[1,3]

# --- no compression of timesteps, verify3

# fName = 'verify_midp_U0_20_nreps_100_nsamples_100.pkl' # seems fine
# fName = 'verify_midp_U0_120_nreps_100_nsamples_100.pkl' # mu=0.05, d=2 constant, quite often has zero detected extinctions; seems fine

# --- better sampling of the two figures to include in the appendix

# verify
# fName = 'verify_midp_U0_150_nreps_500_nsamples_200.pkl' # approx hour run time (9:11 -- > 11:19), I guess the dots are closer to the line?
# fName = 'verify_midp_U0_150_nreps_1000_nsamples_500.pkl' # approx 10 hours -- use for appendix

# verify2
# fName = 'verify_midp_U0_120_nreps_1000_nsamples_500.pkl' # 10:46 -

# -- similar to birds data
# fName = 'verify_midp_U0_48_nreps_50_nsamples_100.pkl' # -- can already see that this is not working well
# fName = 'verify_fbnd_U0_48_nreps_50_nsamples_100.pkl' # -- looks reasonable
fName = 'verify_fbnd_U0_48_nreps_1000_nsamples_500.pkl' # -- verify above, is a good fit

# read in the pickle file with the results
# ---

# directory where results are kept
# dirName = '../../../results/classical/verify/'    # constant detections
# dirName = '../../../results/classical/verify2/'   # choose detections from set of integers
# dirName = '../../../results/classical/verify3/'   # no collapsing of timesteps
dirName = '../../../results/classical/verify4/'   # data similar to birds study

f = open(dirName + fName, 'rb')
ss = pickle.load( f )
params = pickle.load( f )       # dict: parameters used for the run
n = pickle.load( f )            # np array: the sample size at each step
cnt_withinV = pickle.load( f )  # np array of ints: how many times true U0 was within the bounds corresponding to pcileV
U0_meanV = pickle.load( f )     # list of ints: the mean obtained for each rep (rounded)
f.close()


# get the parameter values and results out that we would like to plot
# ---

pcileV = params['pcileV']
nreps = params['nreps']
coverage = 100*cnt_withinV/nreps
U0 = params['U0']


# plot coverage for each confidence level
# ---

plt.plot( pcileV, pcileV, ls='dotted', color='black')
plt.scatter( pcileV, coverage, color='black')
plt.xlabel('nominal coverage desired (%)')
plt.ylabel('actual coverage obtained (%) over ' + str(nreps) + ' simulations')
plt.grid(True)
#plt.show()

plt.tight_layout()
plt.savefig(dirName + fName.split('.pkl')[0] + '.pdf')
plt.close()


# plot frequency histogram of the median U0 estimate
# ---

'''
min_U0 = min(U0_meanV); max_U0 = max(U0_meanV);
bins = np.array(range( min_U0, max_U0+2 )) - 0.5
plt.hist(U0_meanV, color='black', alpha=0.7, bins=bins)
'''

plt.axvline(U0, color='black', ls='dashed', label='true $U_0$')
plt.hist(U0_meanV, color='black', alpha=0.7)
plt.xlabel(r'inferred number of undetected extant species, $\bar{U}_0$')
plt.ylabel('frequency over ' + str(nreps) + ' simulations')
plt.grid(True)
plt.legend(loc='best')
#plt.show()

plt.tight_layout()
plt.savefig(dirName + fName.split('.pkl')[0] + '_histogram.pdf')
plt.close()
