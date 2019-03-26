# obtain a list of common species from Chong et al. (2009) Checklist also found in our db

import csv

import sys
sys.path.insert(0,'../..') # allows us to import undetected extinctions package

from undetected_extinctions.helpers.name_cleaning import clean_species_name


# where the two databases to compare are located
# ---

fname_ourdb = '../../data/cleaned_plants_database/merged.csv'
fname_chong = '../../data/chong_2009_checklist/Chong_native_spp.csv'


# some names in Chong and their equivalent name in our database
# ---

chong2db = {
        'ardisia colorata': 'ardisia sanguinolenta',
        'campnosperma auriculata': 'campnosperma auriculatum',
        'costus speciosus': 'cheilocostus speciosus',
        'clerodendrum laevifolium': 'clerodendrum disparifolium',
        'diospyros lanceaefolia': 'diospyros lanceifolia',
        'embelia lampani': 'embelia lampanii',
        'ficus aurantiaca': 'ficus punctata',
        'gymnanthera nitida': 'gymnanthera oblonga',
        'hypolytrum nemorum': 'hypolytrum nemorum var. proliferum',
        'premna foetida': 'premna serratifolia',
        'tectaria singaporeana': 'tectaria singaporiana',
        'melanthera biflora': 'wollastonia biflora',
}



# get set of standardised names from each
# ---

# our database

fIn = open(fname_ourdb) # open csv
csv_f = csv.reader(fIn) # get row iterator
# [(0, 'gbif ID'), (1, 'catalogue number or barcode'), (2, 'record or collection number'), (3, 'herbarium'),
#  (4, 'collector'), (5, 'co-collectors'), (6, 'collection date'), (7, 'location'), (8, 'family'),
#  (9, 'species'), (10, 'det by'), (11, 'det date'), (12, 'standardise year'), (13, 'standardise name')]
header = next(csv_f)
ourdb = set([ row[13] for row in csv_f ])
fIn.close()

# the Chong checklist
fIn = open(fname_chong) # open csv
csv_f = csv.reader(fIn) # get row iterator
header = next(csv_f)
chong = set([ clean_species_name(row[1]) for row in csv_f if row[5] == 'common' ])
fIn.close()

# convert names where needed
chong_ = [ name if name not in chong2db else chong2db[name] for name in chong ]


# find in our db
# ---

common_in_our_db = [ name for name in chong_ if name in ourdb ]

# write to a csv
# ---

f = open('../../data/processed/common_species_from_chong.csv','w')
for name in sorted(common_in_our_db):
    f.write(name)
    f.write('\n')
f.close()
