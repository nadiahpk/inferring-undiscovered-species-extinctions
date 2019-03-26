# plot the number of detected extant (S) and extinct (E) species in time using the raw data
# three plots, using: 
#   (1) raw first and last detections; 
#   (2) first and last detections + expert info; 
#   (3) first and last detections + expert info + 30 years heuristic.

import csv
import matplotlib.pyplot as plt
import numpy as np


# where the databases are located
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'


# read in the data
# ---

csv_f = csv.reader(open(fname_frstlast)) # open csv
header = next(csv_f)

#              name  , frst_detn  , last_detn  , expert, common
spp_flec = [ ( row[0], int(row[1]), int(row[2]), row[4], row[5] ) for row in csv_f ]


# convert data into first and last detections three ways
# ---

fl_raw = list() # raw first last
fl_exp = list() # raw first last + expert info
fl_heu = list() # raw first last + expert info + heuristic

for name, frst_detn, last_detn, expert, common in spp_flec:

    new_last_detn = last_detn

    # using detection record only

    if last_detn > 2015:
        new_last_detn = 2015

    fl_raw.append( (frst_detn, new_last_detn) )

    # using detection record + expert info

    if expert or common:
        new_last_detn = 2015

    fl_exp.append( (frst_detn, new_last_detn) )

    # using detection record + expert info + heuristic

    if last_detn > 1985:
        new_last_detn = 2015

    fl_heu.append( (frst_detn, new_last_detn) )


# turn first and last detections into running S and E
# ---

flV = [fl_raw, fl_exp, fl_heu]

SV = list()
EV = list()
for frst_last in flV:

    frstObs, lastObs = zip(* frst_last )
    frstObs = np.array(frstObs)
    lastObs = np.array(lastObs)

    # find the S and E
    # ---

    y0 = min( frstObs ) # year of first observation
    yf = max( lastObs ) # year of last observation

    years = range(y0, yf+1) # y0, y0+1, ..., yf inclusive

    # produce S and E from y0-1 to yf

    S = np.array( [ sum( (frstObs <= t) & (lastObs >= t) ) for t in years ] )
    E = np.array( [ sum( lastObs < t ) for t in years ] )

    SV.append(S)
    EV.append(E)


# plot them
# ---

labelV = ['detection record', 'expert knowledge', '30 year heuristic']
lsV = ['dotted', 'dashed', 'solid']

#plt.figure(figsize=(10,5))
for S, E, label, ls in zip( SV, EV, labelV, lsV):
    plt.plot(years, S, 'green',  lw=3, label = r'$S_t$ ' + label, ls=ls, alpha=0.8)
    plt.plot(years, E, 'red',    lw=3, label = r'$E_t$ ' + label, ls=ls, alpha=0.8)

plt.axvline(1985, ls='dotted', color='black') # indicate where the heuristic located
plt.xlabel('year')
plt.ylabel('number of species')
plt.legend(loc='best',ncol=3,fontsize='small')
plt.xlim( (1795,2020) )
plt.ylim( (-50,2075) )
plt.grid(True)
#plt.show()
plt.tight_layout()
plt.savefig('../../results/describe_data/plot_SE_raw.pdf')
plt.close()
