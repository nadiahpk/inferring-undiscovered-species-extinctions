# add another column to our first-last detections file marking if the heuristic, 
# that the species is still extant bc it was last seen within the past 30 years,
# was rejected by the Solow p-value

import csv



# locations of intermediate data
# ---

# first and last detections and initial expert determinations
fname_frstlast = '../../data/processed/first_last_detns.csv' 

# contains both the Solow P-value, and experts' second opinions on whether
# the species designated by the P-value as extinct were in fact extinct
fname_solow_expert = '../../data/inferred_extinct_second_opinion/seen_recently_inferred_extinct_second_expert_opinion.csv'


# read in data
# ---

# first and last detections and other info
csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
fl = { row[0]: row[1:] for row in csv_f } # { standard_name: [ frst, last, no_detns, expert_extant?, chong_common?] }

# experts' second opinions on those designated extinct by Solow P-value
csv_f = csv.reader(open(fname_solow_expert))
header = next(csv_f)
sl = [ row for row in csv_f ] # [ family, spp_name, frst, last, Solow_P, life_form, cons_status, cultivated, expert_2nd_opinion, comment


# for each designated extinct by Solow, append its info to fl
# ---

# make room on each row
for spp_name in fl:

    fl[spp_name].append('') # empty space where I'll put a 'yes' if the Solow test designated it extinct

# for each designated extinct by Solow, (1) mark 'yes' for designated extinct by Solow, and (2) update the expert info.
for row in sl:

    # 1. mark 'yes' because designated extinct by Solow

    spp_name = row[1].lower()
    fl[spp_name][5] = 'yes'

    # 2. update the expert info

    expert_2nd_opinion = row[8].strip().lower()

    if expert_2nd_opinion == 'extant':

        fl[spp_name][3] = 'yes'


# save result to file
# ---

f = open('../../data/processed/first_last_detns_psolow.csv','w')

f.write('standard name,first detection,last detection,no. detns,expert extant?,chong common?,Solow reject heuristic?\n')

for spp_name in sorted(fl):

    row = [spp_name] + fl[spp_name]

    f.write( ','.join(row) )
    f.write( '\n' )

f.close()
