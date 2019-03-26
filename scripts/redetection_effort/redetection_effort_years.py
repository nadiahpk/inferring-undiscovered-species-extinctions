# create a plot of years versus redetection-effort-years, and save to a file for easy use

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import pickle
#import csv


# some parameters for user
# ---

no_params = 182 # NOTE the no. parameters and suffix of file containing the function with the lowest AIC

# time range over which we'll do the redetection effort
t0 = 1822 # the Wallich collection
tf = 2015 # taken as our last date


# where databases are located
# ---

data_dir = '../../results/redetection_effort/tighten/' # where intermediate result kept


# create redetection effort function
# ---

# read in info about best function
fname = data_dir + 'tighten_' + str(no_params) + '.pkl'
f = open(fname, 'rb')
_ = pickle.load( f )   # explanatory string
cs = pickle.load( f )   # c_t points defining the spline
ts = pickle.load( f )   # t points defining the spline
spps = pickle.load( f ) # redetections info for each species used in the fitting
f.close()

# this is a function that accepts a vector of ts (years) and returns the sampling effort in that year
t2c_fnc = interpolate.interp1d( ts, cs, kind='linear')


# obtain cumulative effort-years, save
# ---

tV = np.array(range( t0, tf+1 )) # vector of years
cV = t2c_fnc( tV )
ccV = np.cumsum( cV )

fname = '../../results/redetection_effort/redetection_effort_years.pkl'
f = open(fname, 'wb')
# a string explaining the pickle file
ss  = 'Created by redetection_effort_years.py.\n'
ss += 'Contains the following:\n'
ss += '0. ss, string: this string you are reading now.\n'
ss += '1. tV, np array of integers: years. \n'
ss += '2. ccV: np array of floats: corresponding cumulative effort years. \n'
pickle.dump( ss, f )
pickle.dump( tV, f )
pickle.dump( ccV, f )
f.close()

# plot
# ---

plt.plot(tV, ccV, color='black', lw=3)
plt.xlabel('years')
plt.ylabel('effort-years')
plt.grid(True)
#plt.show()
plt.tight_layout()
plt.savefig('../../results/redetection_effort/redetection_effort_years.pdf')
plt.close()
