# plot the main results figure, the classical and Bayesian SEUX results combined
# double check that the old inverse_midp function worked okay for our
# data, despite the mistake regarding impossible values

import csv
import matplotlib.pyplot as plt
import pickle
import numpy as np


# where are the results stored
# ---

fname_mcmc = '../../results/mcmc/mcmc_basic_result.csv'
fname_classical = '../../results/classical/classical_basic_result.csv'


# read in both
# ---

# mcmc results
csv_f = csv.reader( open(fname_mcmc) )
header = next(csv_f)
m_res = [ [ float(ri) for ri in row ] for row in csv_f ]

# classical results
csv_f = csv.reader( open(fname_classical) )
header = next(csv_f)
c_res = [ [ float(ri) for ri in row ] for row in csv_f ]

years_mod, S, E, c_U_means, c_X_means, c_U_los, c_U_his, c_X_los, c_X_his = zip(* c_res )
_, _, _, m_U_means, m_X_means, m_U_los, m_U_his, m_X_los, m_X_his = zip(* m_res )


# plot them
# ---

# classical result

plt.plot([], [], lw=0, label='Inferred from data:')

plt.plot(years_mod, S,       'green',  lw=2, label = r'detected extant, $S_t$')
plt.plot(years_mod, E,       'red',    lw=2, label = r'detected extinct, $E_t$')

plt.plot([], [], lw=0, label='Classical inference:')

plt.plot(years_mod, c_U_means, 'orange', lw=1, label = r'undetected extant, ${U}_t$')
plt.plot(years_mod, c_X_means, 'blue',   lw=1, label = r'undetected extinct, ${X}_t$')

plt.text( 2016, S[-1], str(int(round(S[-1]))), color='green' , verticalalignment='center', fontsize='x-small')
plt.text( 2016, E[-1], str(int(round(E[-1]))), color='red' , verticalalignment='bottom', fontsize='x-small')
plt.text( 1821, S[0], str(int(round(S[0]))), color='green' , verticalalignment='center', horizontalalignment='right', fontsize='x-small')
plt.text( 2016, c_X_means[-1], str(int(round(c_X_means[-1]))), color='blue' , verticalalignment='top', fontsize='x-small')
plt.text( 1821, c_U_means[0], str(int(round(c_U_means[0]))), color='orange' , verticalalignment='top', horizontalalignment='right', fontsize='x-small')

plt.fill_between(years_mod, c_U_los, c_U_his, facecolor='orange', alpha=0.3)
plt.fill_between(years_mod, c_X_los, c_X_his, facecolor='blue', alpha=0.3)

# mcmc bayesian result

plt.plot([], [], lw=0, label='Bayesian inference:')

plt.plot(years_mod, m_U_means, 'brown', lw=1, label = r'undetected extant, ${U}_t$')
plt.plot(years_mod, m_X_means, 'magenta',   lw=1, label = r'undetected extinct, ${X}_t$')

plt.text( 2016, m_X_means[-1], str(int(round(m_X_means[-1]))), color='magenta' , verticalalignment='center', fontsize='x-small')
plt.text( 1821, m_U_means[0], str(int(round(m_U_means[0]))), color='brown' , verticalalignment='bottom', horizontalalignment='right', fontsize='x-small')

plt.fill_between(years_mod, m_U_los, m_U_his, facecolor='brown', alpha=0.3)
plt.fill_between(years_mod, m_X_los, m_X_his, facecolor='magenta', alpha=0.3)

plt.xlabel('year')
plt.ylabel('number of plant species')
plt.grid(True)
plt.ylim( (-100, 2600) )
plt.xlim( (1807, 2030) )

# sort out legend
leg = plt.legend(loc='upper right', ncol=3, fontsize='x-small')
# leg = plt.legend(handles, labels, ncol=4)

for vpack in leg._legend_handle_box.get_children():
    for i, hpack in enumerate(vpack.get_children()):
        if i % 3 == 0:
            hpack.get_children()[0].set_width(-7)


#plt.show()
plt.tight_layout()
plt.savefig('../../results/main_figure/main_figure.pdf')
plt.close()
