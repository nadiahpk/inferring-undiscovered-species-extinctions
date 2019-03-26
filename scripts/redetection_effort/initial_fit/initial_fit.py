# Obtain an initial redetection effort function, with a point at every year in the timeseries, by iterative solving

import csv
from scipy import interpolate
import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
import pickle

# some parameters
# ---

# time range over which we'll do the redetection effort
t0 = 1822 # the Wallich collection
tf = 2015 # taken as our last date

eps = 1e-6 # a small number used in place of 0 probability of redetection


# where databases are located
# ---

fname_redetns = '../../../results/redetection_effort/redetections_records.pkl'
fname_chosens = '../../../results/redetection_effort/chosen_spp.csv'


# read in databases, filter redetections records for chosen species
# ---

# read in list of species chosen
csv_f = csv.reader(open(fname_chosens))
header = next(csv_f)
spp_chosen = [ row[0] for row in csv_f ]

# get the redetections
f = open(fname_redetns,'rb')
spp_redetns = pickle.load(f) # { spp_name: { 'frst': yr first detected, 'last': yr last detected, 'redetns': list yrs redetected}
f.close()

# filter for our chosen species
spps = { spp_name: D for spp_name, D in spp_redetns.items() if spp_name in spp_chosen }


# create a dictionary that converts the year into an index
# ---

tV = list(range(t0,tf+1))
t2i = { t: i for i, t in enumerate(tV) }


# create the starting x vector, which defines points in the linear spline
# ---

# count how many species were available to be redetected at each t
no_can_redet = [ sum( 1 for spp in spps.values() if t > spp['frst'] and t < spp['last'] ) for t in tV ]

# count how many species were actually redetected at each t
no_did_redet = [ sum( 1 for spp in spps.values() if t in spp['redetns'] ) for t in tV ]

# the first estimate of c_t assumes all species have the same intrinsic redetection probability
# so c_t = number of species redetected at time t / number of species available to be redetected at time t
cV_data = np.array(no_did_redet) / np.array(no_can_redet)


# iteratively solve c_t
# ---

for i in range(50): # trial-and-error finds that 50 iterations is long enough to get a reasonable fit

    # c_t range set to between eps and 1-eps
    cs = [ min(max(eps,c),1-eps) for c in cV_data ]

    # fit a linear spline function to c_t to obtain the function c(t)
    f = interpolate.interp1d( tV, cs, kind='linear')
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


plt.scatter(tV, cV_data, color='red', marker='s', alpha=0.5, label='data')
plt.scatter(tV, cV_fit, color='blue', alpha=0.5)
plt.plot(tV, cV_fit, color='blue', label='fitted $c(t)$')
plt.xlabel('time')
plt.ylabel('redetection effort')
plt.grid(True)
plt.legend(loc='best')
#plt.show()
plt.tight_layout()
plt.savefig('../../../results/redetection_effort/initial_fit/initial_fit.pdf')
plt.close()


# save to pickle file
# ---

fName = '../../../results/redetection_effort/initial_fit/initial_fit.pkl'
f = open(fName, 'wb')
# a string explaining the pickle file
ss  = 'Created by initial_fit.py.\n'
ss += 'Contains the following:\n'
ss += '0. ss, string: this string you are reading now.\n'
ss += '1. cs, list of floats: redetection efforts corresponding to tV.\n'
ss += '2. tV: list ints: years.\n'
ss += '3. spps, dictionary: keys are names, and values are dictionary with frst, last, and redetns.\n'
pickle.dump( ss, f )
pickle.dump( cs, f )
pickle.dump( tV, f )
pickle.dump( spps, f )
f.close()
