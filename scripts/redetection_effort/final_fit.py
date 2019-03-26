# plot a comparison of the final fit to data

import csv
from scipy import interpolate
import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
import pickle

# some parameters for user
# ---

no_params = 182 # NOTE the no. parameters and suffix of file containing the function with the lowest AIC

# time range over which we'll do the redetection effort
t0 = 1822 # the Wallich collection
tf = 2015 # taken as our last date


# where databases are located
# ---

data_dir = '../../results/redetection_effort/tighten/' # where intermediate result kept


# read in best function
# ---

fname = data_dir + 'tighten_' + str(no_params) + '.pkl'
f = open(fname, 'rb')
_ = pickle.load( f )   # explanatory string
cs = pickle.load( f )   # c_t points defining the spline
ts = pickle.load( f )   # t points defining the spline
spps = pickle.load( f ) # redetections info for each species used in the fitting
f.close()


# compare the fit to the data used
# ---

# create a dictionary that converts the year into an index
tV = list(range(t0,tf+1))
t2i = { t: i for i, t in enumerate(tV) }

# fit a linear spline function to c_t to obtain the function c(t)
f = interpolate.interp1d( ts, cs, kind='linear')
cV_fit = f(tV)

# take into account redetn probs to improve last estimate

# calculate the intrinsic redetection probability of each species given the c(t)
# r_i \approx \frac{ \sum_{\tau} I_R(i,\tau) }{ \sum_{\tau \in \mathcal{T}_i} c(\tau) }
redetn_probs = { name: min( 1-1e-6, len(spp['redetns']) / sum( cV_fit[ t2i[max(t0,spp['frst'])] + 1 : t2i[min(tf,spp['last'])] ] ) ) for name, spp in spps.items() }

# no detections by timestep, the numerator \sum_{i \in \mathcal{S}_t} I_R(i,t)
all_redetns = reduce( lambda x,y: x+y, [ spp['redetns'] for spp in spps.values() ] )
no_redetns_bytime = [ all_redetns.count(t) for t in tV ]

# sum of redetection probabilities by timestep, the denominator \sum_{i \in \mathcal{S}_t} r_i(t)
sum_redetn_probs = [ sum( redetn_probs[name] for name in redetn_probs if t > spps[name]['frst'] and t < spps[name]['last'] ) for t in tV ]

# c_t is the ratio, c(t) \approx \frac{ \sum_{i \in \mathcal{S}_t} I_R(i,t) }{ \sum_{i \in \mathcal{S}_t} r_i(t) }
cV_data = [ n/s if s > 0 else np.nan for n, s in zip(no_redetns_bytime, sum_redetn_probs) ]


# plot comparison of fit to data
# ---

plt.figure(figsize=(8*.7,6*.7))
plt.scatter(tV, cV_data, marker='s', color='red', alpha=0.5, label='data')
plt.scatter(ts, cs, color='blue', alpha=0.5)
plt.plot(tV, cV_fit, color='blue', label='fitted $c(t)$')
plt.xlabel('time', fontsize='large')
plt.ylabel('redetection effort', fontsize='large')
plt.grid(True)
plt.legend(loc='best')
#plt.show()
plt.tight_layout()
plt.savefig('../../results/redetection_effort/final_fit.pdf')
plt.close()
