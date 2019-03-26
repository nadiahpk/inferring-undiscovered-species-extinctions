# cumulative number of specimens

import csv
import numpy as np
import matplotlib.pyplot as plt


# where the databases are located
# ---

fname_ourdb = '../../data/cleaned_plants_database/merged.csv'


# count number of specimens in each year
# ---

fIn = open(fname_ourdb) # open csv
csv_f = csv.reader(fIn) # get row iterator
# [(0, 'gbif ID'), (1, 'catalogue number or barcode'), (2, 'record or collection number'), (3, 'herbarium'),
#  (4, 'collector'), (5, 'co-collectors'), (6, 'collection date'), (7, 'location'), (8, 'family'),
#  (9, 'species'), (10, 'det by'), (11, 'det date'), (12, 'standardise year'), (13, 'standardise name')]
header = next(csv_f)

specimens = dict() # { year: count }
for row in csv_f:

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

    else: # if it's a definite date

        year = int(year_s)


    # prepare a place for this year if needed

    if year not in specimens:
        specimens[year] = 0


    # add to count

    specimens[year] += 1


fIn.close()

# calculate the cumulative
# ---

yearV, countV = zip(* sorted( specimens.items(), key=lambda v: v[0] ) )
cumul_countV = np.cumsum(countV)


# plot
# ---

plt.figure(figsize=(4*0.7,3*0.7))
plt.plot(yearV, cumul_countV, color='black', lw=2)
plt.xlim( (1795,2020) )
plt.ylim( (0,35000) )
plt.ylabel( 'specimens' )
plt.xlabel( 'year' )
plt.grid(True)
#plt.show()
plt.tight_layout()
plt.savefig('../../results/describe_data/cumul_specimens.pdf')
plt.close()
