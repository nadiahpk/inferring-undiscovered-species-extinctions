# obtain a basic description of all the species seen in the past 30 years that the Solow p-value designated as extinct
# particularly interested in species designated extinct that are know (from Chong et al. (2009)) to be common

import csv

import sys
sys.path.insert(0,'../..') # allows us to import undetected extinctions package

from undetected_extinctions.helpers.name_cleaning import clean_species_name

# where the databases are located
# ---

fname_spps = '../../results/infer_detected_extinctions/seen_recently_inferred_extinct.csv'
fname_chong = '../../data/chong_2009_checklist/Chong_native_spp.csv'
fname_spp2fam = '../../data/cleaned_plants_database/species2family.csv'
fname_db2chong = '../../data/chong_2009_checklist/database2chong.csv'


# read in the intermediate data
# ---

# create dictionary to map from database name to name in Chong

csv_f = csv.reader(open(fname_db2chong))
header = next(csv_f)
db2chong = { row[0]: row[1] for row in csv_f }

# create dictionary to map species name to family

csv_f = csv.reader(open(fname_spp2fam))
header = next(csv_f)
spp2fam = { row[0]: row[1] for row in csv_f }

# read in species designated extinct their Solow p-value

csv_f = csv.reader(open(fname_spps))
header = next(csv_f)
spp_extincts = { row[0]: { 'frst': row[1], 'last': row[2], 'pvalue': float(row[5]) } for row in csv_f }

# read in chong
csv_f = csv.reader(open(fname_chong))
header = next(csv_f)
spp_chong = { clean_species_name(row[1]): { 'life_form': row[4].strip(), 'status': row[5].strip(), 'cultivated': row[6].strip() } for row in csv_f }



# make rows for a latex table
# ---

table = list() # family, species, first, last, p-value, lifeform, status, cultivated

for spp_name, D in spp_extincts.items():

    fam_name = spp2fam[spp_name]

    table_row = [ fam_name, '\emph{' + spp_name.capitalize() + '}', D['frst'], D['last'], '{:0.4f}'.format(D['pvalue']) ]

    spp_name2 = ' '.join(spp_name.split()[0:2])
    if spp_name in spp_chong:

        C = spp_chong[spp_name2]
        table_row += [ C['life_form'], C['status'], C['cultivated'] ]

    else:

        C = spp_chong[db2chong[spp_name]]
        #table_row += [ '', '', '' ]
        table_row += [ C['life_form'], C['status'], C['cultivated'] ]

    table.append( table_row )


# write table
# ---

f = open('../../results/infer_detected_extinctions/describe_seen_recently_inferred_extinct.tex', 'w')

table.sort(key = lambda r: r[4] ) # want to go in sorted by p-value

for row in table:

    f.write( '\t & '.join(row) )
    f.write( '\\\\ \n' )

f.close()
