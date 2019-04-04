# summary of collections by herbaria

import csv
import numpy as np
import matplotlib.pyplot as plt


# where the databases are located
# ---

fname_ourdb = '../../data/cleaned_plants_database/merged.csv'
fname_spp2fam = '../../data/cleaned_plants_database/species2family.csv'


# create dictionary to convert species name to family name
# ---

csv_f = csv.reader(open(fname_spp2fam))
header = next(csv_f)
spp2famD = { row[0]: row[1] for row in csv_f }


# read through each row of our database, and record info about each herbarium
# ---

herbaria = dict() # { herbarium: { 'no specimens': x, 'species': set of species, 'families': set of families } }

fIn = open(fname_ourdb) # open csv
csv_f = csv.reader(fIn) # get row iterator
# [(0, 'gbif ID'), (1, 'catalogue number or barcode'), (2, 'record or collection number'), (3, 'herbarium'),
#  (4, 'collector'), (5, 'co-collectors'), (6, 'collection date'), (7, 'location'), (8, 'family'),
#  (9, 'species'), (10, 'det by'), (11, 'det date'), (12, 'standardise year'), (13, 'standardise name')]
header = next(csv_f)

for row in csv_f:

    herbarium = row[3].strip()
    spp_name = row[13]
    fam_name = spp2famD[spp_name]

    # append if needed
    if herbarium not in herbaria:
        herbaria[herbarium] = { 'no specimens': 0, 'species': set(), 'families': set() }

    herbaria[herbarium]['no specimens'] += 1
    herbaria[herbarium]['species'].add(spp_name)
    herbaria[herbarium]['families'].add(fam_name)


# there are 63 of these so we'll need to cut it back
# ---

if False:

    # print in order of who has the most
    cnt_h = sorted( [ (herbaria[h]['no specimens'],h) for h in herbaria ], key=lambda v: v[0] )
    for cnt,h in cnt_h:
        print(str(cnt) + ' ' + h )

    '''
    Notes:

        KEW and Kew Herbarium appear separately
        Naturalis may also appear as NHN

    Top no specimens:
        253 NHMUK BOT
        253 E E Royal Botanic Garden Edinburgh Herbarium -- Royal Botanic Garden Edinburgh
        303 US Botany -- United States National Herbarium
        1106 KEW -- The Herbarium at the Royal Botanic Gardens Kew
        2371 Naturalis Botany -- Naturalis Biodiversity Center
        7574 SINU -- The Herbarium of the National University of Singapore 
        21867 SING -- The Singapore Herbarium

    So we'll split into int the top 4 and then "other"
    '''

# take the top 4 and put the rest in other
# ---

key2summary = {
        'SING': 'SING',
        'SINU': 'SINU',
        'Naturalis Botany': 'Naturalis Botany',
        'KEW': 'KEW',
        'Kew Herbarium': 'KEW',
        }

summaryD = { 
        'SING': { 'no specimens': 0, 'species': set(), 'families': set(), 'pretty name': 'Singapore Botanic Gardens Herbarium' },
        'SINU': { 'no specimens': 0, 'species': set(), 'families': set(), 'pretty name': 'Herbarium of the Lee Kong Chian Natural History Museum' },
        'Naturalis Botany': { 'no specimens': 0, 'species': set(), 'families': set(), 'pretty name': 'National Herbarium of the Netherlands' },
        'KEW': { 'no specimens': 0, 'species': set(), 'families': set(), 'pretty name': 'Herbarium at the Royal Botanic Gardens Kew' },
        'Other': { 'no specimens': 0, 'species': set(), 'families': set(), 'pretty name': 'Other' },
    }

for herbarium, D in herbaria.items():

    # get key for where to append
    if herbarium in key2summary:
        key = key2summary[herbarium]
    else:
        key = 'Other'

    # append
    summaryD[key]['no specimens'] += D['no specimens']
    summaryD[key]['species'].update( D['species'] )
    summaryD[key]['families'].update( D['families'] )

# add an entry for "all herbaria"
all_herb = { 'no specimens': 0, 'species': set(), 'families': set(), 'pretty name': 'All combined' }
for herbarium, D in summaryD.items():

    all_herb['no specimens'] += D['no specimens']
    all_herb['species'].update( D['species'] )
    all_herb['families'].update( D['families'] )

summaryD['all'] = all_herb


# now calculate summary statistics
summary = dict()
for h, D in summaryD.items():

    summary[h] = dict()
    summary[h]['no specimens'] = D['no specimens']
    summary[h]['no species'] = len(D['species'])
    summary[h]['no families'] = len(D['families'])
    summary[h]['pretty name'] = D['pretty name']


# write to csv
# ---

# order I want to present them in
order = [ 'SING', 'SINU', 'Naturalis Botany', 'KEW', 'Other', 'all' ]

f = open('../../results/describe_data/summary_collections.csv','w')
f.write('collection,no. specimens,no. species, no. families\n')

for h in order:

    f.write( summary[h]['pretty name'] )
    f.write( ',' + str(summary[h]['no specimens']) )
    f.write( ',' + str(summary[h]['no species']) )
    f.write( ',' + str(summary[h]['no families']) )
    f.write( '\n')

f.close()
