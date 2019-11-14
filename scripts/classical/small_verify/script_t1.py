
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pickle

from small_verify_fncs import get_coverage_estimates, get_example

import sys
sys.path.insert(0,'../../../undetected_extinctions') # so I can import the undetected extinctions package
from undetected_extinctions import find_U0_bnd


# parameters we'll keep constant across all variants
# ---

dName = '../../../results/classical/small_verify/' # directory to write results
pcileV = np.array([90, 80, 70, 60, 50, 40, 30, 20, 10]) # list of percentiles to check coverage for

negs = 100                      # number of example runs to plot

# switches to control behaviour
plotPmu = True
plotegs = True
getCIs = True
plotCIs = True

# create the different runs to perform
# ---

runs = {
        'const': { # completed
            'U0': 50,
            'S0': 10,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'discrete',
            'mu_fnc': lambda t: 0.1,
            'mu_pdf_fnc': lambda mu: 1 if mu == 0.1 else 0,
            'muV': [0.1],
            'nu_fnc_type': 'discrete',
            'nu_fnc': lambda t: 0.1,
            'nu_pdf_fnc': lambda nu: 1 if nu == 0.1 else 0,
            'nuV': [0.1],
            },
        'const2': { # NOTE -- a bit longer, but not a priority
            'U0': 120,
            'S0': 30,
            'nsims': 500,
            'nsamples': 100,
            'T': 5,
            'mu_fnc_type': 'discrete',
            'mu_fnc': lambda t: 0.1,
            'mu_pdf_fnc': lambda mu: 1 if mu == 0.1 else 0,
            'muV': [0.1],
            'nu_fnc_type': 'discrete',
            'nu_fnc': lambda t: 0.1,
            'nu_pdf_fnc': lambda nu: 1 if nu == 0.1 else 0,
            'nuV': [0.1],
            },
        'const3': { # completed
            'U0': 1200,
            'S0': 300,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'discrete',
            'mu_fnc': lambda t: 0.1,
            'mu_pdf_fnc': lambda mu: 1 if mu == 0.1 else 0,
            'muV': [0.1],
            'nu_fnc_type': 'discrete',
            'nu_fnc': lambda t: 0.1,
            'nu_pdf_fnc': lambda nu: 1 if nu == 0.1 else 0,
            'nuV': [0.1],
            },
        'const4': { # completed 
            'U0': 100,
            'S0': 600,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'discrete',
            'mu_fnc': lambda t: 0.15,
            'mu_pdf_fnc': lambda mu: 1 if mu == 0.15 else 0,
            'muV': [0.15],
            'nu_fnc_type': 'discrete',
            'nu_fnc': lambda t: 0.07,
            'nu_pdf_fnc': lambda nu: 1 if nu == 0.07 else 0,
            'nuV': [0.07],
            },
        'beta1': { # completed
            'U0': 50,
            'S0': 10,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'continuous',
            'mu_fnc': lambda t: stats.beta(4,36).rvs(),
            'mu_pdf_fnc': lambda mu: stats.beta(4,36).pdf(mu),
            'muV': np.linspace(0,1,100),
            'nu_fnc_type': 'continuous',
            'nu_fnc': lambda t: stats.beta(4,36).rvs(),
            'nu_pdf_fnc': lambda nu: stats.beta(4,36).pdf(nu),
            'nuV': np.linspace(0,1,100),
            },
        'beta2': { # completed
            'U0': 120,
            'S0': 30,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'continuous',
            'mu_fnc': lambda t: stats.beta(4,36).rvs(),
            'mu_pdf_fnc': lambda mu: stats.beta(4,36).pdf(mu),
            'muV': np.linspace(0,1,100),
            'nu_fnc_type': 'continuous',
            'nu_fnc': lambda t: stats.beta(4,36).rvs(),
            'nu_pdf_fnc': lambda nu: stats.beta(4,36).pdf(nu),
            'nuV': np.linspace(0,1,100),
            },
        'beta3': { # completed 
            'U0': 800,
            'S0': 300,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'continuous',
            'mu_fnc': lambda t: stats.beta(50,450).rvs() if stats.uniform.rvs() < 0.5 else stats.beta(100,400).rvs(),
            'mu_pdf_fnc': lambda mu: 0.5*stats.beta(50,450).pdf(mu) + 0.5*stats.beta(100,400).pdf(mu),
            'muV': np.linspace(0,1,100),
            'nu_fnc_type': 'continuous',
            'nu_fnc': lambda t: stats.beta(50,450).rvs() if stats.uniform.rvs() < 0.5 else stats.beta(100,400).rvs(),
            'nu_pdf_fnc': lambda nu: 0.5*stats.beta(50,450).pdf(nu) + 0.5*stats.beta(100,400).pdf(nu),
            'nuV': np.linspace(0,1,100),
            },
        'fnct': { # completed
            'U0': 120,
            'S0': 30,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: 0.3 if t == 3 else 0.05,
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: 0.1,
            'nuV': [0.1],
            },
        'fnct1': { # completed
            'U0': 1200,
            'S0': 300,
            'nsims': 500,
            'nsamples': 100,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: 0.3 if t == 3 else 0.05,
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: 0.1,
            'nuV': [0.1],
            },
        'fnct2': { # completed -- varying nu in time
            'U0': 800,
            'S0': 300,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: 0.1,
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: 0.3 if t == 3 else 0.05,
            'nuV': [0.1],
            },
        'fnct3': { # completed -- redo of fnct1 for longer
            'U0': 1200,
            'S0': 300,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: 0.3 if t == 3 else 0.05,
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: 0.1,
            'nuV': [0.1],
            },
        'fnct4': { # completed
            'U0': 800,
            'S0': 300,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: 0.1,
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: [0.05, 0.1, 0.2, 0.4, 1][t-1],
            'nuV': [0.1],
            },
        'fnct5': { # completed
            'U0': 800,
            'S0': 300,
            'nsims': 500,
            'nsamples': 100,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: [0.1, 0.1, 0.3, 0.2, .1][t-1],
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: [0.05, 0.4, 0.1, 0.1, .1][t-1],
            'nuV': [0.1],
            },
        'fnct5l': { # completed
            'U0': 800,
            'S0': 300,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: [0.1, 0.1, 0.3, 0.2, .1][t-1],
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: [0.05, 0.4, 0.1, 0.1, .1][t-1],
            'nuV': [0.1],
            },
        'fnct6': { # completed
            'U0': 500,
            'S0': 500,
            'nsims': 500,
            'nsamples': 100,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: [0.3, 0.2, 0.1, 0.1, .1][t-1],
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: [0.05, 0.05, 0.2, 0.4, .8][t-1],
            },
        'fnct6l': { # completed
            'U0': 500,
            'S0': 500,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: [0.3, 0.2, 0.1, 0.1, .1][t-1],
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: [0.05, 0.05, 0.2, 0.4, .8][t-1],
            },
        'fnct7': { # completed
            'U0': 500,
            'S0': 600,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: [0.4, 0.1, 0.4, 0.1, .1][t-1],
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: [0, 0.5, 0.1, 0.1, 1][t-1],
            },
        'fnct8': { # NOTE -- running now, a bit crazy idea
            'U0': 500,
            'S0': 600,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: [0, 0.5, 0.1, 0.1, 0.7][t-1],
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: [0.4, 0.1, 0.4, 0.1, .1][t-1],
            },
        'fnct9': { # NOTE - running now, wobble fnct7 differently
            'U0': 500,
            'S0': 600,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: [0.1, 0.4, 0.1, 0.4, .1][t-1],
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: [0.5, 0, 0.2, 0, 1][t-1],
            },
        'fnct10': { # NOTE - running now - lots of early extinctions, only detect later
            'U0': 500,
            'S0': 600,
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'of_t',
            'mu_fnc': lambda t: [0.6, 0.1, 0.05, 0.05, 0][t-1],
            'nu_fnc_type': 'of_t',
            'nu_fnc': lambda t: [0, 0, 0, 0.3, 1][t-1],
            },
        }

todo_list = ['fnct9']

# for runName, run in runs.items():
for runName in todo_list:

    run = runs[runName]

    U0 = run['U0']
    S0 = run['S0']

    # plot P(mu) and P(nu)
    # ---

    if plotPmu:

        # mu's pdf

        mu_fnc_type = run['mu_fnc_type']

        if mu_fnc_type == 'of_t':

            # should make a plot of t vs mu_t 

            T = run['T']
            mu_fnc = run['mu_fnc']

            tV = np.arange(1,T+1)
            muV = [ mu_fnc(t) for t in tV ]

            plt.figure(figsize=(8*0.3,6*0.3))
            plt.scatter(tV, muV, color='black')

            plt.xlabel(r'$t$', fontsize='x-large')
            plt.ylabel(r'$\mu_t$', fontsize='x-large')
            plt.xticks( tV )
            plt.yticks( [0,1] )
            plt.tight_layout()
            plt.savefig(dName + runName + '_tmu.pdf')
            plt.close()

        else:

            # should make a plot of mu_t vs P(mu_t)

            mu_pdf_fnc = run['mu_pdf_fnc']
            muV = run['muV']

            p_mu = [ mu_pdf_fnc(mu) for mu in muV ]

            plt.figure(figsize=(8*0.3,6*0.3))

            if mu_fnc_type == 'continuous':
                plt.plot(muV, p_mu, 'black', lw=3)
            else:
                plt.stem(muV, p_mu, 'black', markerfmt = 'ko', basefmt=" ")

            plt.xlabel(r'$\mu_t$', fontsize='x-large')
            plt.ylabel(r'$P(\mu_t)$', fontsize='x-large')
            plt.xticks( [0,1] )
            plt.yticks( [0,int(np.ceil(max(p_mu)))] )
            plt.tight_layout()
            plt.savefig(dName + runName + '_Pmu.pdf')
            plt.close()

        # nu's pdf

        nu_fnc_type = run['mu_fnc_type']

        if nu_fnc_type == 'of_t':

            # should make a plot of t vs nu_t 

            T = run['T']
            nu_fnc = run['nu_fnc']

            tV = np.arange(1,T+1)
            nuV = [ nu_fnc(t) for t in tV ]

            plt.figure(figsize=(8*0.3,6*0.3))
            plt.scatter(tV, nuV, color='black')

            plt.xlabel(r'$t$', fontsize='x-large')
            plt.ylabel(r'$\nu_t$', fontsize='x-large')
            plt.xticks( tV )
            plt.yticks( [0,1] )
            plt.tight_layout()
            plt.savefig(dName + runName + '_tnu.pdf')
            plt.close()

        else:

            nu_pdf_fnc = run['nu_pdf_fnc']
            nuV = run['nuV']

            p_nu = [ nu_pdf_fnc(nu) for nu in nuV ]

            plt.figure(figsize=(8*0.3,6*0.3))

            if nu_fnc_type == 'continuous':
                plt.plot(nuV, p_nu, 'black', lw=3)
            else:
                plt.stem(nuV, p_nu, 'black', markerfmt = 'ko', basefmt=" ")

            plt.xlabel(r'$\nu_t$', fontsize='x-large')
            plt.ylabel(r'$P(\nu_t)$', fontsize='x-large')
            plt.xticks( [0,1] )
            plt.yticks( [0,int(np.ceil(max(p_nu)))] )
            plt.tight_layout()
            plt.savefig(dName + runName + '_Pnu.pdf')
            plt.close()


    if plotegs:

        T = run['T']
        mu_fnc = run['mu_fnc']
        nu_fnc = run['nu_fnc']
        fName = dName + runName + '_eg.pdf' # where to save the plot

        S_orig, E_orig, U_orig, X_orig, UV = get_example(U0, S0, mu_fnc, nu_fnc, T, negs)

        plt.figure(figsize=(8*0.5,6*0.5))

        # plot the simulation as "true" values
        plt.plot(S_orig,  'green',  lw=1, label = r'$S_t$')
        plt.plot(E_orig,  'red',    lw=1, label = r'$E_t$')
        plt.plot(X_orig,  'blue',   lw=1, label = r'$X_t$')

        # plot each of the examples
        for i, U in enumerate(UV):

            if i == 0:
                plt.plot(U, 'orange', lw=0.5, alpha=0.5, label = r'$U_t^{[i]}$')
            else:
                plt.plot(U, 'orange', lw=0.5, alpha=0.5)

        # this one last so it's on top
        plt.plot(U_orig,  'black', lw=2, label = r'$U_t$')

        plt.xlabel('year')
        plt.ylabel('number of species')

        plt.legend(loc='upper right', fontsize='small')
        plt.tight_layout()
        plt.savefig(fName)
        plt.close()

    if getCIs:

        nsims = run['nsims']
        nsamples = run['nsamples']
        T = run['T']
        mu_fnc = run['mu_fnc']
        nu_fnc = run['nu_fnc']
        fName = dName + runName + '_coverage.pkl' # where to save the outcome

        cnt_withinV, U0_meanV = get_coverage_estimates(nsims, nsamples, pcileV, U0, S0, T, mu_fnc, nu_fnc)

        # save result to pickle file
        f = open(fName, 'wb')
        pickle.dump( pcileV, f )
        pickle.dump( nsims, f )
        pickle.dump( cnt_withinV, f )
        pickle.dump( U0_meanV, f )
        f.close()

    if plotCIs:

        fName = dName + runName + '_coverage.pkl' # where to save the outcome

        f = open(fName, 'rb')
        pcileV = pickle.load( f )
        nsims = pickle.load( f )
        cnt_withinV = pickle.load( f )
        U0_meanV = pickle.load( f )
        f.close()

        coverage = 100*cnt_withinV/nsims

        plt.figure(figsize=(8*0.5,6*0.5))

        plt.plot( pcileV, pcileV, ls='dotted', color='black')
        plt.scatter( pcileV, coverage, color='black')
        plt.xlabel('nominal coverage desired')
        plt.ylabel('actual coverage obtained')

        fName = dName + runName + '_coverage.pdf' # where to save the outcome

        plt.xlim( (0,100) )
        plt.ylim( (0,100) )
        plt.tight_layout()
        plt.savefig(fName)
        plt.close()


'''
If I have time
        'const2': {
            'nsims': 1000,
            'nsamples': 1000,
            'T': 5,
            'mu_fnc_type': 'discrete',
            'mu_fnc': lambda: 0.1 if stats.uniform.rvs() < 0.5 else 0.05,
            'mu_pdf_fnc': lambda mu: 0.5 if mu == 0.1 else ( 0.5 if mu == 0.05 else 0 ),
            'muV': [0.05, 0.1],
            'nu_fnc_type': 'discrete',
            'nu_fnc': lambda: 0.1 if stats.uniform.rvs() < 0.5 else 0.05,
            'nu_pdf_fnc': lambda nu: 0.5 if nu == 0.1 else ( 0.5 if nu == 0.05 else 0 ),
            'nuV': [0.05, 0.1],
            },
'''
