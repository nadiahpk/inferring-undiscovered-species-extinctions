# plot the main results figure, the classical and Bayesian SEUX results combined

import csv
import matplotlib.pyplot as plt
import pickle
import numpy as np


# where are the results stored
# ---

dname_resV = [
        '../../../results/classical/sensitivity/sensitivity_UT/',
        '../../../results/classical/sensitivity/sensitivity_species_deletion/',
        '../../../results/classical/sensitivity/sensitivity_species_deletion/'
        ]

fname_resV = [
        'example_UT.csv',
        'example_species_deletion.csv',
        'example_species_deletion_UT0.csv'
        ]

for dname_res, fname_res in zip(dname_resV, fname_resV):

    # read in results
    # ---

    csv_f = csv.reader( open(dname_res + fname_res) )
    header = next(csv_f)
    res = [ [ float(ri) for ri in row ] for row in csv_f ]

    years_mod, S, E, U_means, X_means, U_los, U_his, X_los, X_his = zip(* res )


    # plot them
    # ---

    plt.plot(years_mod, S,       'green',  lw=2, label = r'detected extant, $S_t$')
    plt.plot(years_mod, E,       'red',    lw=2, label = r'detected extinct, $E_t$')
    plt.plot(years_mod, U_means, 'orange', lw=1, label = r'inferred undetected extant, $\bar{U}_t$')
    plt.plot(years_mod, X_means, 'blue',   lw=1, label = r'inferred undetected extinct, $\bar{X}_t$')

    plt.text( 2016, S[-1], str(int(round(S[-1]))), color='green' , verticalalignment='center', fontsize='x-small')
    plt.text( 2016, E[-1], str(int(round(E[-1]))), color='red' , verticalalignment='bottom', fontsize='x-small')
    plt.text( 1821, S[0], str(int(round(S[0]))), color='green' , verticalalignment='center', horizontalalignment='right', fontsize='x-small')
    plt.text( 2016, X_means[-1], str(int(round(X_means[-1]))), color='blue' , verticalalignment='top', fontsize='x-small')
    plt.text( 1821, U_means[0], str(int(round(U_means[0]))), color='orange' , verticalalignment='top', horizontalalignment='right', fontsize='x-small')

    plt.fill_between(years_mod, U_los, U_his, facecolor='orange', alpha=0.3)
    plt.fill_between(years_mod, X_los, X_his, facecolor='blue', alpha=0.3)

    #plt.legend(loc='upper right', ncol=2, fontsize='x-small')
    plt.xlabel('year')
    plt.ylabel('number of plant species')
    plt.grid(True)
    #plt.ylim( (-100, 2600) )
    plt.xlim( (1807, 2030) )
    #plt.show()
    plt.tight_layout()
    plt.savefig(dname_res + fname_res.split('.')[0] + '.pdf')
    plt.close()

# plot legend only
# ---

if False:

    import pylab

    fig = pylab.figure()
    figlegend = pylab.figure()
    ax = fig.add_subplot(111)
    lineS = ax.plot([0,1], [0,1], 'green',  lw=2)
    lineE = ax.plot([0,1], [0,1], 'red',    lw=2)
    lineU = ax.plot([0,1], [0,1], 'orange', lw=1)
    lineX = ax.plot([0,1], [0,1], 'blue',   lw=1)
    legend = figlegend.legend( 
                                (lineS[0], lineE[0], lineU[0], lineX[0]) , 
                                (
                                    r'detected extant, $S_t$', 
                                    r'detected extinct, $E_t$', 
                                    r'inferred undetected extant, $\bar{U}_t$', 
                                    r'inferred undetected extinct, $\bar{X}_t$'),
                                'center',
                                fontsize='large'
                                )

    #fig.show()
    #figlegend.show()
    figlegend.canvas.draw()
    #figlegend.savefig('../../../results/classical/sensitivity/examples_legend.pdf',
        #bbox_inches=legend.get_window_extent().transformed(figlegend.dpi_scale_trans.inverted()))
    figlegend.savefig('../../../results/classical/sensitivity/examples_legend.pdf')
