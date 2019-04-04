# cumulative number of species, of families

import csv
import numpy as np
import matplotlib.pyplot as plt


# where the databases are located
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'
fname_spp2fam = '../../data/cleaned_plants_database/species2family.csv'


# read in species names and year of first detection
# ---

csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
spp_yrs = [ (row[0], int(row[1])) for row in csv_f ]


# read in the dictionary from species name to family name
# ---

csv_f = csv.reader(open(fname_spp2fam))
header = next(csv_f)
spp2famD = { row[0]: row[1] for row in csv_f }


# create list of families seen
# ---

famS = set([ spp2famD[spp_yr[0]] for spp_yr in spp_yrs ])


# turn into family and year of first detection
# ---

# create a dictionary of { fam: years seen }
fam_yrD = { fam: 3000 for fam in famS } # initialise

for spp, yr in spp_yrs:

    fam = spp2famD[spp]

    if yr < fam_yrD[fam]:

        fam_yrD[fam] = yr


# find
# ---

_, year_sppV = zip(* spp_yrs )
year_famV = list(fam_yrD.values())


yearV = list(range(1795,2020+1))
cnt_sppV = [ year_sppV.count(yr) for yr in yearV ]
cnt_famV = [ year_famV.count(yr) for yr in yearV ]

cumul_cnt_sppV = np.cumsum(cnt_sppV)
cumul_cnt_famV = np.cumsum(cnt_famV)


# plot years with new species
# ---

end_spp = yearV.index(max(year_sppV))
end_fam = yearV.index(max(year_famV))

plt.figure(figsize=(4*0.7,3*0.7))
plt.plot(yearV[:end_spp+1], cumul_cnt_sppV[:end_spp+1], color='black', lw=2)
plt.xlim( (1795,2020) )
plt.ylim( (0,2200) )
plt.ylabel( 'species' )
plt.xlabel( 'year' )
plt.grid(True)
#plt.show()
plt.tight_layout()
plt.savefig('../../results/describe_data/cumul_species.pdf')
plt.close()

plt.figure(figsize=(4*0.7,3*0.7))
plt.plot(yearV[:end_fam+1], cumul_cnt_famV[:end_fam+1], color='black', lw=2)
plt.xlim( (1795,2020) )
plt.ylim( (0,180) )
plt.ylabel( 'families' )
plt.xlabel( 'year' )
plt.grid(True)
#plt.show()
plt.tight_layout()
plt.savefig('../../results/describe_data/cumul_families.pdf')
plt.close()


# write to csv in case we need particular values
# ---

f = open('../../results/describe_data/cumul_species_families.csv','w')
f.write('year,cumulative no. species,cumulative no. families\n')
for year, cnt_spp, cnt_fam in zip(yearV, cumul_cnt_sppV, cumul_cnt_famV):
    f.write(str(year) + ',' + str(cnt_spp) + ',' + str(cnt_fam) + '\n')
f.close()
