# create a table of the extinction rate per family

import csv
import matplotlib.pyplot as plt

# where databases are kept
# ---

fname_species2family = '../../data/cleaned_plants_database/species2family.csv' # header: standardise name,family
fname_extants = '../../data/processed/first_last_detns.csv' # header: standard name,first detection,last detection


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


# create a dictionary of families with value list: [ no_extinct, no_total ]
# ---

famD = dict()
for spp_name, status in spp2status.items():

    family = spp2fam[spp_name]

    if family not in famD:
        famD[family] = [0,0]

    # add species to total
    famD[family][1] += 1

    # add to extinct count if extinct
    if status == 'extinct':
        famD[family][0] += 1


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
plt.tight_layout()
plt.savefig('../../results/describe_data/extn_rate_per_family.pdf')
plt.close()



'''
# check
# ---

cnt_extinct = 0
cnt_extant = 0
for spp_name, status in spp2status.items():

    if status == 'extant':
        cnt_extant += 1
    else:
        cnt_extinct += 1
'''
