# calculate the intrinsic redetection probability for all species that we can using the final-fit redetection effort function

import csv
import pickle
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt


# some fixed parameters
# ---

no_params = 182 # NOTE the no. parameters and suffix of file containing the function with the lowest AIC

t0 = 1822 # the Wallich collection
tf = 2015 # taken as our end date


# where the databases are located
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'
fname_detns = '../../data/processed/detections_records.pkl'
fname_redetn_effort = '../../results/redetection_effort/tighten/tighten_' + str(no_params) + '.pkl'


# turn definite-detections years into a redetection record
# ---

# get extant species based on the info we had before we did the Solow p-value check
csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
# extant if last detected since 1985, experts say extant, or common
extants = [ row[0] for row in csv_f if int(row[2]) > 1985 or row[4] == 'yes' or row[5] == 'yes' ]

# get the detections
f = open(fname_detns,'rb')
spp_detns = pickle.load(f)
f.close()

# turn into redetections
spp_redetns = dict()
for name, records in spp_detns.items():

    detns = records['definite'] # only include detections with a definite date

    if detns: # if this is not empty

        frst_detn = min(detns)

        if name in extants:

            last_detn = tf+1 # note counting the 2015 detection if occurs
            redetns = list(filter( lambda x: x <= tf and x >= t0, detns[1:] ))

        else:

            last_detn = max(detns)
            redetns = list(filter( lambda x: x <= tf and x >= t0, detns[1:-1] )) # excludes last detection

        if redetns: # if this is not empty we can include

            spp_redetns[name] = dict()
            spp_redetns[name]['frst'] = frst_detn
            spp_redetns[name]['redetns'] = redetns
            spp_redetns[name]['last'] = last_detn


# create function for redetection effort
# ---

f = open(fname_redetn_effort, 'rb')
_ = pickle.load( f )   # explanatory string
cs = pickle.load( f )   # c_t points defining the spline
ts = pickle.load( f )   # t points defining the spline
f.close()

# create a dictionary that converts the year into an index
tV = list(range(t0,tf+1))
t2i = { t: i for i, t in enumerate(tV) }

# fit a linear spline function to c_t to obtain the function c(t)
f = interpolate.interp1d( ts, cs, kind='linear')
cV_fit = f(tV)


# calculate each species' intrinsic redetection probability
# ---

# calculate the intrinsic redetection probability of each species given the c(t)
# r_i \approx \frac{ \sum_{\tau} I_R(i,\tau) }{ \sum_{\tau \in \mathcal{T}_i} c(\tau) }
redetn_probs = {
        name: min( 1-1e-6, len(spp['redetns']) / sum( cV_fit[ t2i[max(t0,spp['frst'])] + 1 : t2i[min(tf,spp['last'])] ] ) )
        for name, spp in spp_redetns.items() }


# plot and write to csv
# ---

# histogram of intrinsic redetection probabilities

# split into values from definitely extant and maybe extinct species
extant_redetnV = [ redetn_prob for name, redetn_prob in redetn_probs.items() if name in extants ]
extinct_redetnV = [ redetn_prob for name, redetn_prob in redetn_probs.items() if name not in extants ]

delta = 0.05
bins = np.arange(-delta/2, 1+delta, delta)

plt.figure(figsize=(.7*8, .7*6))
plt.hist([ extinct_redetnV, extant_redetnV], bins=bins,
        label=['presumed extinct', 'presumed extant'],
        color=['black','gray'], alpha=0.7, stacked=False)
plt.legend(loc='best')
plt.ylabel('number of species', fontsize='large')
plt.xlabel(r'species intrinsic redetection probability $r_i$', fontsize='large')
plt.xlim( (-delta, 1+delta) )
#plt.show()
plt.tight_layout()
plt.savefig('../../results/redetection_effort/intrinsic_redetection_probabilities_histogram.pdf')
plt.close()

# how intrinsic redetection probabilities vary with year of discovery

# create a dictionary D = { year_of_discovery: [ list of intrinsic redetection probs ] }
D = dict()
for name, redetn_prob in redetn_probs.items():

    frst = min( spp_detns[name]['definite'] + spp_detns[name]['inferred'] )
    if frst not in D:
        D[frst] = list()

    D[frst].append( redetn_prob )

# now average them and make a sorted list
frst_mean_redetnD = { frst: np.mean(redetnV) for frst, redetnV in D.items() }
frst_redetnV = sorted( frst_mean_redetnD.items(), key=lambda v: v[0] )
frst_detnV, redetnV = zip(* frst_redetnV )

# plot scatter and convolution
plt.scatter(frst_detnV, redetnV, alpha=0.7, label='year-averaged')

# convolve for running average
NV = [10, 20]
lsV = ['solid','dashed']

for N, ls in zip(NV, lsV):

    '''
    # one option is to pad the edges
    redetnV_padded = np.pad(redetnV, (N//2, N-1-N//2), mode='edge')
    redetn_convV = np.convolve(redetnV_padded, np.ones((N,))/N, mode='valid')
    plt.plot(frst_detnV, redetn_convV) # for valid method
    '''

    # but I prefer to chop off the plot, it's more obvious what's going on
    redetn_convV = np.convolve(redetnV, np.ones((N,))/N, mode='valid')
    plt.plot(frst_detnV[N//2:-(N-1-N//2)], redetn_convV, lw=3, ls=ls, label='window width ' + str(N) + ' running mean')


plt.xlabel('year of species\' first detection')
plt.ylabel('intrinsic redetection probability averaged over species')
plt.ylim( (0,1.2) )
plt.grid(True)
plt.legend(loc='upper center', ncol=2)
#plt.show()
plt.tight_layout()
plt.savefig('../../results/redetection_effort/intrinsic_redetection_probabilities_vs_frstdetn.pdf')
plt.close()

# csv file

f = open('../../results/redetection_effort/intrinsic_redetection_probabilities.csv','w')
f.write('standard name,intrinsic redetection probability\n')

for name in sorted(redetn_probs.keys()):

    redetn_prob = redetn_probs[name]
    f.write( name + ',' + str(redetn_prob) + '\n')

f.close()

