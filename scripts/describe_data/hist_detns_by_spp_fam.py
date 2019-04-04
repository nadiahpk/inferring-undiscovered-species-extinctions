# plot a histogram of the no. of detections by species and by families

import csv
import numpy as np
import matplotlib.pyplot as plt


# where the databases are located
# ---

fname_frstlast = '../../data/processed/first_last_detns.csv'
fname_spp2fam = '../../data/cleaned_plants_database/species2family.csv'


# read in species names and no. of detections
# ---

csv_f = csv.reader(open(fname_frstlast))
header = next(csv_f)
spp_nodetns = [ (row[0], int(row[3])) for row in csv_f ]


# turn into family and no. of detections
# ---

# read in the dictionary from species name to family name
csv_f = csv.reader(open(fname_spp2fam))
header = next(csv_f)
spp2famD = { row[0]: row[1] for row in csv_f }


# create list of families seen
# ---

famS = set([ spp2famD[spp_nodetn[0]] for spp_nodetn in spp_nodetns ])

# create a dictionary of { fam: years seen }
fam_nodetns = { fam: 0 for fam in famS } # initialise

for spp, nodetns in spp_nodetns:

    fam = spp2famD[spp]
    fam_nodetns[fam] += nodetns


# create the list of no. of detections by spp and by family
# ---

_, nodetns_spp = zip(* spp_nodetns )
nodetns_fam = list(fam_nodetns.values())


# plot histograms
# ---

# by species

bins = np.arange( min(nodetns_spp)-0.5, max(nodetns_spp)+1.5, 1 )

plt.figure(figsize=(.7*8, .7*6))
plt.hist(nodetns_spp, bins=bins, color='black', alpha=0.7)
plt.ylabel('frequency', fontsize='large')
plt.xlabel('no. of detections of species', fontsize='large')
#plt.show()
plt.tight_layout()
plt.savefig('../../results/describe_data/hist_detns_by_spp.pdf')
plt.close()

# by families

# I've turned off the fine binning here bc histogram v. spread out
# bins = np.arange( min(nodetns_fam)-0.5, max(nodetns_fam)+1.5, 1 ) 

plt.figure(figsize=(.7*8, .7*6))
plt.hist(nodetns_fam, color='black', alpha=0.7)
plt.ylabel('frequency', fontsize='large')
plt.xlabel('no. of detections of family', fontsize='large')
#plt.show()
plt.tight_layout()
plt.savefig('../../results/describe_data/hist_detns_by_fam.pdf')
plt.close()
