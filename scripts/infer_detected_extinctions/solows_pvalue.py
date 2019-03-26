# find the p-values for extinction and the estimated extinction date
# using the Solow method corrected for redetection effort

from scipy import interpolate
import pickle


# user parameters
# ---

t0 = 1822
tf = 2015
tmin = t0 # the earliest record included in the Solow calculation


# locations of intermediate data
# ---

fname_effortyears = '../../results/redetection_effort/redetection_effort_years.pkl'
fname_detns = '../../data/processed/detections_records.pkl'


# create functions that map from year to effort-year and back again
# ---

# read in year and effort-year vectors
f = open(fname_effortyears, 'rb')
_ = pickle.load( f )   # explanatory string
tV = pickle.load( f )  # years
ccV = pickle.load( f ) # effort-years
f.close()

# limit of the range of our effort-years
cc0 = ccV[0]; ccf = ccV[-1]

# functions to convert between, linear splines that obtain values within the range
t2cc = interpolate.interp1d( tV, ccV, kind='linear' )
cc2t = interpolate.interp1d( ccV, tV, kind='linear' )


# find Solow's p-value for each species where possible
# ---

# get species detections records
f = open(fname_detns, 'rb')
spp_detns = pickle.load( f )
f.close()

# loop through species, finding p-value

spp_ps = dict() # place to store

for name, D in spp_detns.items():

    detns = D['definite'] # only using definite detections

    if detns: # if non-empty

        if max(detns) >= tf:

            spp_ps[name] = 1 # p-value set to 1 if seen on or after our last date

        else:

            ts = list(filter( lambda t: t >= tmin and t <= tf , detns )) # filter only detections within timeseries range
            n = len(ts)-1

            if n >= 1: # need more than one detection to apply Solow method

                # convert detection record in years to record in effort-years
                ccs = t2cc(ts) 

                # apply the Solow equation in terms of effort-years
                p_solow = ( (ccs[-1]-ccs[0]) / (ccf-ccs[0]) )**n

                # save our result to the dictionary
                spp_ps[name] = p_solow


# save p-value results to a csv file
# ---

f = open('../../results/infer_detected_extinctions/solows_pvalue.csv', 'w')
f.write('standard name,Solow p-value\n')

for name in sorted(spp_ps.keys()):

    p_solow = spp_ps[name]
    f.write( name + ',' + str(p_solow) + '\n')

f.close()
