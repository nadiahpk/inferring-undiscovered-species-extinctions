# add another column to our first-last detections file marking if the heuristic, 
# that the species is still extant bc it was last seen within the past 30 years,
# was rejected by the Solow p-value

import csv



# locations of intermediate data
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'
fname_rejects = '../../results/infer_detected_extinctions/seen_recently_inferred_extinct.csv'


# read in data
# ---

# first and last detections and other info
csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
fl = [ row for row in csv_f ]

# list of species rejected as extant according to Solow's P < 0.1
csv_f = csv.reader(open(fname_rejects))
header = next(csv_f)
rejects = [ row[0] for row in csv_f ]


# append the Solow p-value information as needed
for row in fl:

    spp_name = row[0]

    if spp_name in rejects:

        row.append('yes')

    else:

        row.append('')


f = open('../../data/processed/first_last_detns_psolow.csv','w')

f.write('standard name,first detection,last detection,no. detns,expert extant?,chong common?,Solow reject heuristic?\n')

for row in fl:

    f.write( ','.join(row) )
    f.write( '\n' )

f.close()
