# create the redetections records, which is detections apart from first and last detection

import csv
import pickle

# some fixed parameters
# ---

t0 = 1822 # the Wallich collection
tf = 2015 # taken as our end date


# where the databases are located
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'
fname_detns = '../../data/processed/detections_records.pkl'


# identify species extant after 2015 (these I can count detection at 2015 for)
# ---

csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
extants = [ row[0] for row in csv_f if int(row[2]) > 2015 or row[4] ] # so either seen after 2015 or experts' say it is extant


# turn definite detections years into a redetection record
# ---

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


# output files of redetections records
# ---

# pickle file

f = open('../../results/redetection_effort/redetections_records.pkl','wb')
pickle.dump( spp_redetns, f )
f.close()

