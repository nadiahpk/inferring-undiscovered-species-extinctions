# plot a summary of the backwards stepping, the number of parameters versus the AIC found, and tightening

import pickle
import numpy as np
import matplotlib.pyplot as plt


# which to plot
# ---

# results from incrementally removing points

bkwds_frst_no_params = 120
bkwds_last_no_params = 386
bkwds_results_dir = '../../results/redetection_effort/backwards_stepping/' # where are results

# results from tightening the above

tight_results_dir = '../../results/redetection_effort/tighten/' # where are results
tight_frst_no_params = 158
tight_last_no_params = 224


# create a list of file name suffixes
# ---

bkwds_suffixes = np.arange(bkwds_frst_no_params, bkwds_last_no_params+2, 2) # no. of parameters is used as suffix
tight_suffixes = np.arange(tight_frst_no_params, tight_last_no_params+2, 2) # increments by two


# read in each results file and store info
# ---

# backwards-stepping incremental removal of points

# ready storage
bkwds_AICV = list()
bkwds_no_paramsV = list()

# loop over suffixes
for suffix in bkwds_suffixes:

    # read in info from that file

    fname = bkwds_results_dir + 'backwards_stepping_' + str(suffix) + '.pkl'
    f = open(fname, 'rb')
    _ = pickle.load( f )   # explanatory string
    cs = pickle.load( f )   # c_t points defining the spline
    ts = pickle.load( f )   # t points defining the spline
    _ = pickle.load( f ) # redetections info for each species used in the fitting
    AIC = pickle.load( f ) # redetections info for each species used in the fitting
    f.close()

    # store info
    no_params = len(cs) + len(ts)
    bkwds_AICV.append( AIC )
    bkwds_no_paramsV.append( no_params )

# tightening of backwards-stepping results

# ready storage
tight_AICV = list()
tight_no_paramsV = list()

# loop over suffixes
for suffix in tight_suffixes:
#for suffix in [172, 176, 178, 180, 184]:

    # read in info from that file

    fname = tight_results_dir + 'tighten_' + str(suffix) + '.pkl'
    f = open(fname, 'rb')
    _ = pickle.load( f )   # explanatory string
    cs = pickle.load( f )   # c_t points defining the spline
    ts = pickle.load( f )   # t points defining the spline
    _ = pickle.load( f ) # redetections info for each species used in the fitting
    AIC = pickle.load( f ) # redetections info for each species used in the fitting
    f.close()

    # store info
    no_params = len(cs) + len(ts)
    tight_AICV.append( AIC )
    tight_no_paramsV.append( no_params )


# plot
# ---

plt.scatter( bkwds_no_paramsV, bkwds_AICV, color='black', s=5, label='from incremental removal of points')
plt.scatter( tight_no_paramsV, tight_AICV, color='red', s=5, label='from maximimum likelihood fitting')

plt.xlabel( 'number of parameters\n' + r'($t$ and $c_t$ values defining redetections effort spline)')
plt.ylabel( 'AIC' )
# plt.title( 'from backwards-stepping algorithm' )
plt.legend( loc='best' )
plt.grid( True )
# plt.show()
plt.tight_layout()
plt.savefig('../../results/redetection_effort/summary_AIC_minimisation.pdf')
plt.close()

# -- plot says the minimum is no_params = 182
