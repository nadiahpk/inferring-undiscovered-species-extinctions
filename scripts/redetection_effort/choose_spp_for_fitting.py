# choose which species' redetection records will be used for determining the redetection effort in time
# exclude species that:
#   1. have last detection between 1985 and 2015
#   2. have no redetections between t0 and tf
#   3. have records shorter than min_lifespan = 30 years

from functools import reduce
import csv
import pickle
import matplotlib.pyplot as plt


# some parameters
# ---

min_lifespan = 30

t0 = 1822
tf = 2015


# where databases are located
# ---

fname_redetns = '../../results/redetection_effort/redetections_records.pkl'



# get the redetections
# ---

f = open(fname_redetns,'rb')
spp_redetns = pickle.load(f) # { spp_name: { 'frst': yr first detected, 'last': yr last detected, 'redetns': list yrs redetected}
f.close()


# create the list of species chosen
# ---

spp_chosen = list()
for spp_name, D in spp_redetns.items():

    # pull out info about redetections

    frst = D['frst']
    last = D['last']
    redetns = D['redetns']

    if not ( last >= 1985 and last <= 2015 ): # exclude species that have last detection between 1985 and 2015

        lifespan = last-frst

        if lifespan >= min_lifespan: # exclude species that have records shorter than min_lifespan = 30 years

            if list(filter(lambda x: x >= t0 and x <= tf, redetns)): # exclude species that have no redetections between t0 and tf

                spp_chosen.append( spp_name )


# plot coverage
# ---

coverage = { t: 0 for t in range(t0, tf+1) }

for spp_name in spp_chosen:

    frst = spp_redetns[spp_name]['frst']
    last = spp_redetns[spp_name]['last']

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
plt.savefig('../../results/redetection_effort/chosen_spp_coverage_redetections.pdf')
plt.close()


# write list to csv
# ---

f = open('../../results/redetection_effort/chosen_spp.csv','w')
f.write('standard names\n')

for name in sorted(spp_chosen):

    f.write( name + '\n')

f.close()
