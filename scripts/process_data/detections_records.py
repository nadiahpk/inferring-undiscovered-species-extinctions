# write all detections to detns, and first-last to first_last_detns.csv

import csv
import pickle


# where the databases are located
# ---

fname_ourdb = '../../data/cleaned_plants_database/merged.csv'


# read in stitched_3, creating detns a dictionary of species standard names and detection years
# ---

fIn = open(fname_ourdb) # open csv
csv_f = csv.reader(fIn) # get row iterator
# [(0, 'gbif ID'), (1, 'catalogue number or barcode'), (2, 'record or collection number'), (3, 'herbarium'),
#  (4, 'collector'), (5, 'co-collectors'), (6, 'collection date'), (7, 'location'), (8, 'family'),
#  (9, 'species'), (10, 'det by'), (11, 'det date'), (12, 'standardise year'), (13, 'standardise name')]
header = next(csv_f)

detns = dict()
for row in csv_f:

    # prepare a place for this name if needed
    name = row[13]
    if name not in detns:
        detns[name] = { 'inferred': set(), 'definite': set() }

    # get the year
    year_s = row[12]

    if 'c' in year_s: # if it's uncertain, inferred from collector

        year_s = year_s.replace('c','').replace('.','')
        years = year_s.split('-')

        if len(years) == 1:
            year = int(years[0])
        else:
            min_yr = int(years[0]); max_yr = int(years[1])
            year = round(min_yr + 0.5*(max_yr-min_yr))

        # append
        detns[name]['inferred'].add(year)

    else: # if it's a definite date

        year = int(year_s)

        # append
        detns[name]['definite'].add(year)

fIn.close()

# sort them for convenience
detns_sorted = { spp_name: { 'inferred': sorted( v['inferred'] ), 'definite': sorted( v['definite'] ) } for spp_name, v in detns.items() }


# output files of detections records
# ---

# pickle file

f = open('../../data/processed/detections_records.pkl','wb')
pickle.dump( detns_sorted, f )
f.close()

# csv file

f = open('../../data/processed/detections_records.csv','w')
f.write('standard name\n')

for name in sorted(detns_sorted.keys()):

    f.write( name )

    v = detns_sorted[name]
    all_detns = sorted( set( v['inferred'] + v['definite'] ) )

    for detn in all_detns:

        f.write( ',' + str(detn) )

    f.write('\n')

f.close()

