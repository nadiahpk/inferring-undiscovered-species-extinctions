import matplotlib.pyplot as plt
import pickle
import numpy as np
import csv


# name which file
# ---

fname = '../../../results/classical/sensitivity/sensitivity_species_deletion/sensitivity_spp_deletion_all.pkl'

xlabel = 'proportion of species included'
proportions = np.arange(0.3,1,0.1)
percentile = 95 # CI that I want


# load result
# ---

f = open(fname, 'rb')
extn_rateM = pickle.load( f )
f.close()

# and classical result

fname_classical = '../../../results/classical/classical_basic_result.csv'
csv_f = csv.reader( open(fname_classical) )
header = next(csv_f)
c_res = [ [ float(ri) for ri in row ] for row in csv_f ]
years_mod, S, E, U, X, _, _, _, _, = zip(* c_res )

N = S[0] + E[0] + U[0]              # assumes X(0) = 0
extn_rate = (N-S[-1]-U[-1]) / N     # calculate extinction rate


# plot
# ---

mean = np.mean(extn_rateM, axis=1)
lo = np.percentile(extn_rateM, (100-percentile)/2, axis=1)
hi = np.percentile(extn_rateM, 100 - (100-percentile)/2, axis=1)

mean2 = np.append( mean, [extn_rate] )
lo2 = np.append( lo, [extn_rate] )
hi2 = np.append( hi, [extn_rate] )
proportions2 = np.append( proportions, [1] )

plt.figure(figsize=(4*0.7,3*0.7))
plt.plot(proportions2, mean2, color="black")
plt.fill_between(proportions2, lo2, hi2, facecolor='black', alpha=0.5)
plt.xlabel('proportion of species included')
plt.ylabel('extinction rate')
plt.ylim( (0, 1) )
plt.grid(True)
# plt.show()
plt.tight_layout()
plt.savefig('../../../results/classical/sensitivity/sensitivity_species_deletion/sensitivity_species_deletion.pdf')
plt.close()

