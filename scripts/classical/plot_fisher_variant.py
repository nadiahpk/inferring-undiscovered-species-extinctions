# plot the main results figure, the classical and Bayesian SEUX results combined

import csv
import matplotlib.pyplot as plt
import pickle
import numpy as np


# where are the results stored
# ---

# fname_classical = '../../results/classical/classical_fisher_184.csv'
# fname_classical = '../../results/classical/classical_fisher_180.csv'
# fname_classical = '../../results/classical/classical_fisher_160.csv'
fname_classical = '../../results/classical/classical_fisher_170.csv'
fname_classical = '../../results/classical/classical_fisher_172.csv'


# read in both
# ---

# classical results
csv_f = csv.reader( open(fname_classical) )
header = next(csv_f)
c_res = [ [ float(ri) for ri in row ] for row in csv_f ]

years_mod, S, E, c_U_means, c_X_means, c_U_los, c_U_his, c_X_los, c_X_his = zip(* c_res )


# plot them
# ---

# classical result

plt.plot(years_mod, S,       'green',  lw=2, label = r'detected extant, $S_t$')
plt.plot(years_mod, E,       'red',    lw=2, label = r'detected extinct, $E_t$')
plt.plot(years_mod, c_U_means, 'orange', lw=1, label = r'classical inferred undetected extant, $\bar{U}_t$')
plt.plot(years_mod, c_X_means, 'blue',   lw=1, label = r'classical inferred undetected extinct, $\bar{X}_t$')

plt.text( 2016, S[-1], str(int(round(S[-1]))), color='green' , verticalalignment='center', fontsize='x-small')
plt.text( 2016, E[-1], str(int(round(E[-1]))), color='red' , verticalalignment='bottom', fontsize='x-small')
plt.text( 1821, S[0], str(int(round(S[0]))), color='green' , verticalalignment='center', horizontalalignment='right', fontsize='x-small')
plt.text( 2016, c_X_means[-1], str(int(round(c_X_means[-1]))), color='blue' , verticalalignment='top', fontsize='x-small')
plt.text( 1821, c_U_means[0], str(int(round(c_U_means[0]))), color='orange' , verticalalignment='top', horizontalalignment='right', fontsize='x-small')

plt.fill_between(years_mod, c_U_los, c_U_his, facecolor='orange', alpha=0.3)
plt.fill_between(years_mod, c_X_los, c_X_his, facecolor='blue', alpha=0.3)

plt.legend(loc='upper right', ncol=2, fontsize='x-small')
plt.xlabel('year')
plt.ylabel('number of plant species')
plt.grid(True)
#plt.ylim( (-100, 2600) )
plt.xlim( (1807, 2030) )
plt.tight_layout()
plt.savefig(fname_classical.split('csv')[0] + 'pdf')
plt.close()
