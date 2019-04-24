# how does the coverage of our space vary as we change the minimum lifespan included?

import pickle
import matplotlib.pyplot as plt


# some parameters
# ---

t0 = 1822 # the Wallich collection
tf = 2015 # taken as our end date


# where the databases are located
# ---

fname_redetns = '../../results/redetection_effort/redetections_records.pkl'


# get the redetections
# ---

f = open(fname_redetns,'rb')
spp_redetns = pickle.load(f) # { spp_name: { 'frst': yr first detected, 'last': yr last detected, 'redetns': list yrs redetected}
f.close()


# count the 'coverage' in terms of number of species given a cut-off min_lifespan

min_lifespan = 30
coverage = { t: 0 for t in range(t0, tf+1) }

cnt_spp = 0
for spp, D in spp_redetns.items():

    frst = D['frst']
    last = D['last']

    lifespan = last-frst
    if lifespan >= min_lifespan:
        cnt_spp += 1

        for t in range( max(t0,frst+1), min(last,tf+1) ):
            coverage[t] += 1

tV,cV = zip(* sorted( coverage.items() , key = lambda v: v[0] ) )

plt.plot(tV, cV, lw=2, color='black')
plt.axvline(t0, color='black', ls='dotted')
plt.axvline(tf, color='black', ls='dotted')
plt.xlabel('year', fontsize="x-large")
plt.ylim( (0,1400) )
plt.text(tV[0]+5, cV[0]+10, '(' + str(tV[0]) + ',' + str(cV[0]) + ')')
plt.ylabel('number of species with\nrecords covering year', fontsize="x-large")
plt.grid(True)
#plt.show()
plt.tight_layout()
plt.savefig('../../results/redetection_effort/coverage_redetections.pdf')
plt.close()

'''
# frequency plot of lifespans
# ---

lifespans = [ v['last']-v['frst'] for v in spp_redetns.values() ]
bins = [ i-0.5 for i in range( max(lifespans)+2 ) ]
plt.hist(lifespans, bins=bins)
plt.show()

# frequency plot of number of redetections
# ---

cnt_redetnsV = [ len(v['redetns']) for v in spp_redetns.values() ]
bins = [ i-0.5 for i in range( min(cnt_redetnsV), max(cnt_redetnsV)+2 ) ]
plt.hist(cnt_redetnsV, bins=bins)
plt.show()

'''

