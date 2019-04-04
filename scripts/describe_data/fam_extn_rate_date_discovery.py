# create a table of the extinction rate per family

import csv
import matplotlib.pyplot as plt
import numpy as np

# where databases are kept
# ---

fname_species2family = '../../data/cleaned_plants_database/species2family.csv' # header: standardise name,family
fname_extants = '../../data/processed/first_last_detns.csv' # header: standard name,first detection,last detection
fname_frstdetn = '../../data/processed/detections_records.csv' # header: standard name,first detection,second detection


# read in species names and families, and create dictionary from species name to family
# ---

csv_f = csv.reader(open( fname_species2family ))
header = next(csv_f)
spp2fam = { row[0]: row[1] for row in csv_f }


# read in species names and date of last detection (after accounting for expert info etc)
# and create dictionary from species name to status (extant/extinct)
# ---

csv_f = csv.reader(open( fname_extants ))
header = next(csv_f)

spp2status = dict()
for row in csv_f:

    spp_name = row[0]
    last_detn = int(row[2])
    expert_extant = row[4]
    common = row[5]

    if last_detn >= 1985 or expert_extant == 'yes' or common == 'yes':

        spp2status[spp_name] = 'extant'

    else:

        spp2status[spp_name] = 'extinct'

# read in species names and date of discovery
# ---

csv_f = csv.reader(open( fname_frstdetn ))
header = next(csv_f)
spp2disc = { row[0]: int(row[1]) for row in csv_f }


# create a dictionary of families with value list: [ no_extinct, no_total, date_disc ]
# ---

famD = dict()
for spp_name, status in spp2status.items():

    family = spp2fam[spp_name]

    if family not in famD:
        famD[family] = [0,0,2015]

    # add species to total
    famD[family][1] += 1

    # add to extinct count if extinct
    if spp2status[spp_name] == 'extinct':
        famD[family][0] += 1

    # earliest discovery
    disc = spp2disc[spp_name]
    if disc < famD[family][2]:
        famD[family][2] = disc


# write a list: [ date_of_discovery, proportion_extinct, total_no_spp ]
# ---

fam_disc_prop_tot = [ [ fam, v[2], v[0]/v[1], v[1] ] for fam, v in famD.items() ]   # all families
filt_fam_disc_prop_tot = [ v for v in fam_disc_prop_tot if v[2] > 0 ]               # only families with some extinctions


# create a scatter plot of the date of first dis
# ---

# pull out data for plotting
fams, discs, props, tots = zip(* fam_disc_prop_tot )
f_fams, f_discs, f_props, f_tots = zip(* filt_fam_disc_prop_tot )

# sort out marker size, proportional to spp richness of famiily
scaler = 4
s = np.array(tots) * scaler

# plt.grid(True)
plt.scatter(discs, props, color='black', edgecolor='none', alpha=0.7, s=s)

plt.xlabel('year of discovery of family')
plt.ylabel('proportion of species in family extinct')

# create the legend showing relative marker sizes
plt.xlim( (1790, 2010) )
plt.scatter(3000, 1, color='black', edgecolor='none', alpha=0.7, s=100*scaler, label='100 species')
plt.scatter(3000, 1, color='black', edgecolor='none', alpha=0.7, s=50*scaler, label='50 species')
plt.scatter(3000, 1, color='black', edgecolor='none', alpha=0.7, s=10*scaler, label='10 species')
plt.scatter(3000, 1, color='black', edgecolor='none', alpha=0.7, s=1*scaler, label='1 species')
leg = plt.legend(loc='upper center', title=r'family size:', ncol=4, bbox_to_anchor=(0.5, 1.30), borderpad=0.7)
leg._legend_box.align = "left"
#plt.ylim( (-0.05, 1.2) )

# plt.show()
plt.tight_layout(rect=[0,0,1,.99])
plt.savefig('../../results/describe_data/fam_extn_rate_date_discovery.pdf')
plt.close()

'''
# write to quadruple list, ordered by speciousness and alphabetically secondarily
# ---

fam_extincts_extants_total = [ [fam, v[0], v[1]-v[0], v[1]] for fam, v in famD.items() ]
fam_extincts_extants_total.sort( key=lambda v: (v[3], v[0]) )


# plot
# ---

fams, extincts, extants, totals = zip(* fam_extincts_extants_total )
inds = list(range(len(fams)))
bar_width = 0.5

plt.figure(figsize=(12,17))

p1 = plt.barh(inds, extincts, bar_width, color='red', alpha=0.8)
p2 = plt.barh(inds, extants, bar_width, left=extincts, color='blue', alpha=0.8)

plt.ylabel('families')
plt.xlabel('no. species per family')
plt.yticks(inds, fams, fontsize='x-small')
plt.ylim(-1,len(fams))
plt.legend((p1[0], p2[0]), ('no. extinct species', 'no. extant species'), loc='lower right')

# plt.show()
plt.grid(axis='x')



'''
