import matplotlib.pyplot as plt
import pickle
import numpy as np


# name which file
# ---

fname = '../../../results/classical/sensitivity/sensitivity_UT/sensitivity_UT.pkl'
xlabel = 'assumed undetected extant at end of record $U_{2015}$'
UTV = [0, 100, 200, 300, 500] # NOTE
percentile = 95 # CI that I want


# load result
# ---

f = open(fname, 'rb')
extn_rateM = pickle.load( f )
f.close()


# plot
# ---

mean = np.mean(extn_rateM, axis=1)
lo = np.percentile(extn_rateM, (100-percentile)/2, axis=1)
hi = np.percentile(extn_rateM, 100 - (100-percentile)/2, axis=1)

plt.figure(figsize=(4*0.7,3*0.7))
plt.plot(UTV, mean, color="black")
plt.fill_between(UTV, lo, hi, facecolor='black', alpha=0.5)
plt.xlabel('assumed $U_{2015}$')
plt.ylabel('extinction rate')
plt.ylim( (0, 1) )
plt.grid(True)
# plt.show()
plt.tight_layout()
plt.savefig(fname.split('pkl')[0] + 'pdf')
plt.close()
