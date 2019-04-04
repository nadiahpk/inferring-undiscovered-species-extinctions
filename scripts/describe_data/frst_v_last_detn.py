# plot a scatter diagram of first versus last detection date for each species

import csv
import numpy as np
import matplotlib.pyplot as plt


# where the databases are located
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'


# read in frst, last detection tuples
# ---

csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
frst_lasts = [ ( int(row[1]), int(row[2]) ) for row in csv_f ]

frsts, lasts = zip(* frst_lasts )

# plot
# ---

fig = plt.figure(figsize=(.7*8, .7*6))
ax = fig.add_subplot(111)
ax.scatter(frsts, lasts, s=2, alpha=0.5)
ax.set_xlabel('year of first detection (discovery)')
ax.set_ylabel('year of last detection')
ax.set_aspect('equal', 'box')
#plt.show()
plt.tight_layout()
plt.savefig('../../results/describe_data/frst_v_last_detn.pdf')
plt.close()
