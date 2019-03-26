# take a backwards_removes fit and tighten the fit

# from logLikelihoods import logLikelihoodSpline

import csv
from scipy import interpolate
import numpy as np
from functools import reduce
import operator # so I can mul(tiply) things
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import pickle

# calculates loglikelihood 
def logLikelihoodSpline(x, t0, tf, spps, t2i):

    tV = list(range(t0,tf+1))
    x = list(x) # comes in as np array

    # I have to split x into real meanings of times and collection efforts
    c_end = x[-2:]
    t_end = [t0, tf]
    halfway = int( (len(x)-2)/2 )
    ts = np.array( x[:halfway] + t_end )
    cs = np.array( x[halfway:-2] + c_end )

    # create the linear spline
    f = interpolate.interp1d( ts, cs, kind='linear')
    cV_fit = f(tV)

    # get the detn prob estimate for each species
    redetn_probs = { name: min( 1-1e-6, len(spp['redetns']) / sum( cV_fit[ t2i[max(t0,spp['frst'])] + 1 : t2i[min(tf,spp['last'])] ] ) ) for name, spp in spps.items() }

    # now calculate the likelihood of the data given the c and detns_prob_ests

    rV = list() # store probabilities at each timestep
    for t, c in zip(tV, cV_fit):

        # the probability of each timestep of data being how it was if c is true
        pp = list()
        for name, spp in spps.items():

            redetn_prob = redetn_probs[name]

            if t > spp['frst'] and t < spp['last']: # only counts those that were detected and extant

                if t in spp['redetns']:

                    pp.append( redetn_prob*c )

                else:

                    pp.append( 1-redetn_prob*c )

        r = reduce(operator.mul, pp, 1)

        rV.append(r)

    LL = -sum(np.log(rV))

    if False:
        print(x)
        print(LL)
        print('-----')

    return LL

# --------------------------------'

# some parameters
# ---

# time range over which we'll do the redetection effort
t0 = 1822 # the Wallich collection
tf = 2015 # taken as our last date
tV = list(range(t0,tf+1))
# dictionary making time to index of tV
t2i = { t: i for i, t in enumerate(tV) }
eps = 1e-6 # so redetection effort doesn't hit 0 or 1


# file naming
suffix = '1'
max_no_intermediate_pts = 81

# get list of expert-extant species' names
# ---

fInName = '../../../plants_databases/merged/first_last_detns.csv'
fIn = open(fInName) # open csv
csv_f = csv.reader(fIn) # get row iterator

# [(0, 'standard name'), (1, 'first detection'), (2, 'last detection'), (3, 'no. detns'), (4, 'expert extant?')]
header = next(csv_f)

extants = list()
for row in csv_f:

    if int(row[2]) >= 2015 or row[4] == 'yes': # extant

        extants.append(row[0])

fIn.close()


# get a list of species we'll use for fitting the redetections effort
# ---

fIn = open('choose_spp/all/choose_spp.csv') # open csv
csv_f = csv.reader(fIn) # get row iterator
header = next(csv_f) # [(0, 'standard name'), (1, 'lifespan')]
names = [ row[0] for row in csv_f ]
fIn.close()

# create a redetections dictionary
# ---

# get the species dictionaries with detections data
f = open('../../../plants_databases/merged/detns.pkl','rb')
spp_detns = pickle.load(f)
f.close()

spps = dict()
for name in names:

    spps[name] = dict()
    detns = spp_detns[name]['definite']
    spps[name]['frst'] = min(detns)

    if name in extants:

        spps[name]['last'] = tf+1
        # note counting the 2015 detection if occurs
        spps[name]['redetns'] = list(filter( lambda x: x <= tf and x >= t0, detns[1:] )) 

    else:

        spps[name]['last'] = max(detns)
        spps[name]['redetns'] = list(filter( lambda x: x <= tf and x >= t0, detns[1:-1] )) # excludes last detection


for suffix in range(75,105):

    # get the backwards_removes we want to tighten
    # ---

    fName = 'backwards_removes/backwards_removes_' + str(suffix) + '.pkl'
    f = open(fName, 'rb')
    ss = pickle.load( f )
    x = pickle.load( f )
    AIC = pickle.load( f )
    f.close()

    no_intermediate_pts = int( (len(x)-2)/2 )

    # I have to split x into real meanings of times and collection efforts
    no_intermediate_pts = int( (len(x)-2)/2 )
    ts = x[:no_intermediate_pts]
    cs = x[no_intermediate_pts:] # note includes c_t0 and c_tf at end

    # tighter bounds
    bnds = list()
    for idx in range(no_intermediate_pts+2):

        ci = cs[idx]
        bnd = ( max(eps, ci-0.05), min(1-eps, ci+0.05) )
        bnds.append(bnd)

    res = minimize( lambda c: logLikelihoodSpline( ts + list(c), t0, tf, spps, t2i ), cs, method='L-BFGS-B', bounds=bnds )

    cs = res.x
    new_x = ts + list(cs)

    LL = logLikelihoodSpline(new_x, t0, tf, spps, t2i)
    new_AIC = 2*len(new_x) + 2*LL

    print('suffix: ' + str(suffix) + ' old AIC: ' + str(round(AIC)) + ' new AIC: ' + str(round(new_AIC)) )

    # write to file
    # ---

    fName = 'backwards_removes/tighten_' + str(suffix) + '.pkl'
    f = open(fName, 'wb')
    # a string explaining the pickle file
    ss  = 'Created by tighten_removes.py.\n'
    ss += '0. ss, string: this string you are reading now.\n'
    ss += '1. x, list of floats and ints: Can be used as a new start value. \n'
    ss += '2. AIC, float: AIC for that x.\n'
    pickle.dump( ss, f ) # 0.
    pickle.dump( new_x, f )
    pickle.dump( new_AIC, f )
    f.close()

'''
if False: # check first point by plotting

        
    # I have to split x into real meanings of times and collection efforts
    c_end = x[-2:]
    t_end = [t0, tf]
    halfway = int( (len(x)-2)/2 )
    ts = np.array( x[:halfway] + t_end )
    cs = np.array( x[halfway:-2] + c_end )

    # create the linear spline
    f = interpolate.interp1d( ts, cs, kind='linear')
    cV_fit = f(tV)

    plt.scatter(ts, cs, marker='x')
    plt.plot(tV, cV_fit, label='AIC = ' + str(AIC) )


    # plot data on last go only

    # get the detn prob estimate for each species
    redetn_probs = { name: min( 1-1e-6, len(spp['redetns']) / sum( cV_fit[ t2i[max(t0,spp['frst'])] + 1 : t2i[min(tf,spp['last'])] ] ) ) for name, spp in spps.items() }

    # no detections by timestep
    all_redetns = reduce( lambda x,y: x+y, [ spp['redetns'] for spp in spps.values() ] )
    no_redetns_bytime = [ all_redetns.count(t) for t in tV ]

    sum_redetn_probs = [ sum( redetn_probs[name] for name in redetn_probs if t > spps[name]['frst'] and t < spps[name]['last'] ) for t in tV ]
    cV_data = [ n/s if s > 0 else np.nan for n, s in zip(no_redetns_bytime, sum_redetn_probs) ]

    plt.scatter(tV, cV_data, alpha=0.3)


    plt.legend(loc='best')
    plt.grid(True)
    plt.show()
'''
