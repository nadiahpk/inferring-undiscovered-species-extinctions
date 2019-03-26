# create a list of species who were last detected in the past 30 years (a heuristic) but are inferred,
# by the Solow p-value, to be extinct

import csv


# user parameters
# ---

tf = 2015
p_cutoffs = [ 0.01, 0.05, 0.1 ] # three cut-offs for the p-value explored


# locations of intermediate data
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'
fname_solowps = '../../results/infer_detected_extinctions/solows_pvalue.csv'


# for each species satisfying the heuristic, if Solow's p-value is below cut-off, make a note of it
# ---

# read in species names and info about them
csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
spp_infos = { row[0]: { 'frst': int(row[1]), 'last': int(row[2]), 'experts': row[4], 'chong_common': row[5] } for row in csv_f }

# read in species names and their Solow p-value
csv_f = csv.reader(open(fname_solowps))
header = next(csv_f)
spp_ps = { row[0]: row[1] for row in csv_f }

# the heuristic is that species who were last seen in the past 30 years are presumed extant
heur = tf - 30 # 1985 in our case

# loop through species

spp_infer = dict()
for name, info in spp_infos.items():

    last = info['last']

    if last > heur and last < tf: # seen in last 30 years but not known from detn record to be extant

        if name in spp_ps: # if we were able to apply the Solow method to it

            pvalue = float(spp_ps[name]) # Solow p-value

            extinct = list() # list of whether extinct according to various p-values
            extinctFlag = False
            for p_cutoff in p_cutoffs:

                if pvalue < p_cutoff:

                    extinct.append('yes')
                    extinctFlag = True

                else:

                    extinct.append('no')

            if extinctFlag: # was inferred extinct according to at least one of our p-value cutoffs

                # store info about this species
                frst = info['frst']
                experts = info['experts']
                common = info['chong_common']

                spp_infer[name] = {
                        'frst': frst,
                        'last': last,
                        'pvalue': pvalue,
                        'extincts': extinct,
                        'experts': experts,
                        'common': common }

# write results to a file
# ---

f = open('../../results/infer_detected_extinctions/seen_recently_inferred_extinct.csv', 'w')

# headers
f.write('standard name,first detection,last detection,expert extant?,chong common?,Solow p-value')
for p_cutoff in p_cutoffs:
    f.write(',extinct at P-crit ' + str(p_cutoff))
f.write('\n')

# results
for name in sorted(spp_infer.keys()):

    v = spp_infer[name]
    f.write( name + ',' + str(v['frst']) + ',' + str(v['last']) + ',' + v['experts'] + ',' + v['common'] + ',' + str(v['pvalue']) + ',' + ','.join(v['extincts']) + '\n' )

f.close()
