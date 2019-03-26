# use first and last detections, t0-tf range, expert knowledge, heuristic, and Solow p-value
# to make our final record of first and last detections for every species

import csv


# parameters
# ---

t0 = 1822 # first year to include
tf = 2015 # last year to include


# where databases are
# ---

fname_frstlast = '../../data/processed/first_last_detns_psolow.csv'


# turn first-last detections database into a list of first and last observations
# ---

# read in first and last detections and other info
csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
fl_info = [ row for row in csv_f ]

fl = dict()
for row in fl_info:

    spp_name = row[0]
    frst = int(row[1])
    last = int(row[2])

    # update the first detection if needed
    if frst < t0:

        frst = t0

    # update the last detection as neeed

    if row[4] == 'yes' or row[5] == 'yes': # experts say extant or Chong says common

        last = 2015

    else:

        if last > tf-30: # if extant according to the heuristic

            if row[6] != 'yes': # do not reject heuristic according to the Solow P-value

                last = 2015

    fl[spp_name] = (frst, last)


# write to csv file
# ---

f = open('../../data/processed/first_last_detns_final.csv','w')

f.write('standard name,first detection,last detection\n')

for name in sorted(fl.keys()):

    frst = str( fl[name][0] )
    last = str( fl[name][1] )
    f.write( name + ',' + frst + ',' + last + '\n')


f.close()
