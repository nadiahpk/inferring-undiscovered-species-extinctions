# write the first and last detections, whether experts said it was extant, and if Chong says it is common

import csv

import sys
sys.path.insert(0,'../..') # allows us to import undetected extinctions package

from undetected_extinctions.helpers.name_cleaning import clean_species_name


# where the databases are located
# ---

fname_ourdb = '../../data/cleaned_plants_database/merged.csv'
fname_experts_extant = '../../data/cleaned_plants_database/experts_extant.csv'
fname_common = '../../data/processed/common_species_from_chong.csv'


# read in experts' extants and commons
# ---

# experts' extants

fIn = open(fname_experts_extant) # open csv
csv_f = csv.reader(fIn) # get row iterator
header = next(csv_f)
extants = [ clean_species_name(row[1]) for row in csv_f ]
fIn.close()

# Chong's common species

fIn = open(fname_common) # open csv
csv_f = csv.reader(fIn) # get row iterator
# header = next(csv_f) # has no header
commons = [ row[0] for row in csv_f ]
fIn.close()


# read in our database and get each detection
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
            year = round(min_yr + 0.5*(max_yr-min_yr)) # take average

        # append
        detns[name]['inferred'].add(year)

    else: # if it's a definite date

        year = int(year_s)

        # append
        detns[name]['definite'].add(year)

fIn.close()

if False: # for double-checking

    common_not_ourdb = [ name for name in commons if name not in detns ]
    # -- empty

    extant_not_ourdb = [ name for name in extants if name not in detns ]
    # -- returns the following list, 
    # 'aglaia tenuicaulis',         # new discovery not added to herbarium databases yet
    # 'aglaia yzermannii',          # new discovery not added to herbarium databases yet
    # 'ampelocissus cinnamomea',    # new discovery not added to herbarium databases yet
    # 'ficus stricta',              # 2 records no barcode removed
    # 'hanguana neglecta',          # complicated reidentification not explained by experts, removed
    # 'hanguana nitens',            # new discovery not added to herbarium databases yet
    # 'hanguana rubinea',           # new discovery not added to herbarium databases yet
    # 'hanguana triangulata',       # new discovery not added to herbarium databases yet
    # 'lygodium circinnatum',       # sent to lygodium longifolium in redownloaded database
    # 'plectocomiopsis geminiflora' # 1 record no barcode removed
    # 'zingiber singapurense'       # new discovery not added to herbarium databases yet


# get first-last detections using all data (definite + inferred)
# ---

# get first-last
fl = dict()
for name, v in detns.items():

    all_detns = set(list(v['inferred']) + list(v['definite']))
    fl[name] = ( min(all_detns), max(all_detns), len(all_detns) )


# write to csv file
# ---

f = open('../../data/processed/first_last_detns.csv','w')

f.write('standard name,first detection,last detection,no. detns,expert extant?,chong common?\n')

for name in sorted(fl.keys()):

    frst = str( fl[name][0] )
    last = str( fl[name][1] )
    no_detns = str( fl[name][2] )
    f.write( name + ',' + frst + ',' + last + ',' + no_detns )

    if name in extants:
        f.write(',yes')
    else:
        f.write(',')

    if name in commons:
        f.write(',yes\n')
    else:
        f.write(',\n')

f.close()
