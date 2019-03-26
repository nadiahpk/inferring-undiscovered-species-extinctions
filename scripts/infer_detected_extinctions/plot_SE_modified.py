# plot the number of detected extant (S) and extinct (E) species in time using the raw data,
# our expert and heuristic info, and the Solow's p-value for species to which the heuristic applies

import csv
import matplotlib.pyplot as plt
import numpy as np


# user parameters
# ---

t0 = 1822
tf = 2015
p_cutoffs = [ 0.01, 0.05, 0.1 ] # three cut-offs for the p-value explored
lsV = ['dotted', 'dashed', 'solid'] # line style for each


# where the databases are located
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'
fname_solowps = '../../results/infer_detected_extinctions/solows_pvalue.csv'


# read in the intermediate data
# ---

# read in species names and their Solow p-value

csv_f = csv.reader(open(fname_solowps))
header = next(csv_f)
spp_ps = { row[0]: float(row[1]) for row in csv_f }

# read in information about species first and last detection, and if experts think they are extant
csv_f = csv.reader(open(fname_frstlast)) # open csv
header = next(csv_f)
#              name  , frst_detn  , last_detn  , expert, common
spp_flec = [ ( row[0], int(row[1]), int(row[2]), row[4], row[5] ) for row in csv_f ]


# for each critical p-value, for each species, 
# - give it a new date of last detection based on expert knowledge, the heuristic, and Solow's p-value
# - plot the resulting figure
# ---

# for each p-value cut-off

SV = list()
EV = list() # place to store

for p_cutoff in p_cutoffs:

    # for each species, 

    fl = list() # place to store results

    for name, frst_detn, last_detn, expert, common in spp_flec:

        new_last_detn = last_detn

        if last_detn > 2015 or expert or common: # if known to be extant

            new_last_detn = 2015

        else:

            if last_detn > tf-30: # if unsure but heuristic applies to this species

                if name not in spp_ps: # we couldn't calculate the Solow p-value for this species

                    new_last_detn = 2015 # so we default to the null hypothesis that it is still extant

                else:

                    p_solow = spp_ps[name]
                    
                    if p_solow >= p_cutoff: # then we say it is still extant

                        new_last_detn = 2015


        # store the species' new info
        fl.append( (frst_detn, new_last_detn) )


    # turn first and last detections into running S and E

    frstObs, lastObs = zip(* fl )
    frstObs = np.array(frstObs)
    lastObs = np.array(lastObs)

    y0 = min( frstObs ) # year of first observation
    yf = max( lastObs ) # year of last observation
    years = range(y0, yf+1) # y0, y0+1, ..., yf inclusive

    S = np.array( [ sum( (frstObs <= t) & (lastObs >= t) ) for t in years ] )
    E = np.array( [ sum( lastObs < t ) for t in years ] )

    SV.append(S)
    EV.append(E)


# plot them and save figure
# ---

labelV = [ '$P < ' + str(p_cutoff) + '$' for p_cutoff in p_cutoffs ]

#plt.figure(figsize=(10,5))
for S, E, label, ls in zip( SV, EV, labelV, lsV):
    plt.plot(years, S, 'green',  lw=2, label = r'$S_t$ with ' + label, ls=ls, alpha=0.8)
    plt.plot(years, E, 'red',    lw=2, label = r'$E_t$ with ' + label, ls=ls, alpha=0.8)

plt.axvline(1985, ls='dotted', color='black') # indicate where the heuristic located
plt.xlabel('year')
plt.ylabel('number of species')
plt.legend(loc='best',ncol=3,fontsize='small')
plt.xlim( (1795,2020) )
plt.ylim( (-50,2075) )
plt.grid(True)
#plt.show()
plt.tight_layout()
plt.savefig('../../results/infer_detected_extinctions/plot_SE_modified.pdf')
plt.close()
