# can we get a sense of how the date of discovery affects the probability of going extinct?
# plot the proportion of currently-extinct vs the year of discovery

import csv
import numpy as np
import matplotlib.pyplot as plt

# params
# ---

t0 = 1796; tf = 1984


# where the databases are located
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'

# read in year of discovery, whether presumed extant or extinct, put in list in dictionary
# ---

csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)

D = dict()
for row in csv_f:

    frst = int(row[1])
    last = int(row[2])

    # is it presumed extinct or extant?
    
    #  expert_extant   or common          or last detected in last 30 years
    if row[4] == 'yes' or row[5] == 'yes' or last >= 1985:
        extinct = 0
    else:
        extinct = 1

    # append to our dictionary
    if frst in D:
        D[frst].append(extinct)
    else:
        D[frst] = [extinct]


# calculate the average proportion of species extinct by year of discovery
# ---

frst, propn_extinct = zip(* [ (frst, np.mean(D[frst])) for frst in sorted(D.keys()) if frst < 1985 ] )

fig = plt.figure(figsize=(.7*8, .7*6)) # plotting as we go
ax = fig.add_subplot(111)

ax.scatter(frst, propn_extinct, s=5, color='black', label='year-averaged')

# the above is obviously hard to read because sometimes there aren't that many species detected in a year
# use a rolling window

NV = [10, 30]
labelV = ['10 yr window', '30 yr window']

for N, label in zip(NV, labelV):

    N2 = int(N/2)
    yr_mids = [ t0+N2, tf-N2, N]
    yr_midV = list(range(yr_mids[0], yr_mids[1]+1))

    avgs = list()
    for yr in yr_midV:

        biglist = list()
        for yri in range( yr-N2, yr+N2+1):
            if yri in D:
                biglist += D[yri]

        # take average
        avgs.append(np.mean(biglist))

    ax.plot(yr_midV, avgs, label=label)

ax.legend(loc='best')
ax.set_xlabel('year of first detection (discovery)')
ax.set_ylabel('proportion of species presumed extinct')
#plt.show()
plt.tight_layout()
plt.savefig('../../results/describe_data/propn_extinct_v_discovery.pdf')
plt.close()
